#include "eager/bvcnf.hpp"

#include "features.hpp"
#include "logger.hpp"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <sys/stat.h>

namespace ttc
{
namespace eager
{

namespace
{

bool fileExists(const std::string& path)
{
  struct stat st;
  return ::stat(path.c_str(), &st) == 0;
}

// Locate an approxmc binary. Prefer one on PATH, then known install locations.
std::string locateApproxMc()
{
  const char* env = std::getenv("APPROXMC");
  if (env != nullptr && fileExists(env))
  {
    return env;
  }
  for (const std::string& cand : {std::string(std::getenv("HOME") ? std::getenv("HOME") : "")
                                      + "/bins/approxmc",
                                  std::string("/usr/local/bin/approxmc"),
                                  std::string("/usr/bin/approxmc")})
  {
    if (!cand.empty() && fileExists(cand))
    {
      return cand;
    }
  }
  // Fall back to PATH lookup via the bare name.
  return "approxmc";
}

}  // namespace

BvCnfResult writeBvCnf(cvc5::Solver& solver,
                       const std::vector<cvc5::Term>& assertions,
                       const std::vector<cvc5::Term>& projectionVars,
                       const std::string& path)
{
  auto& tm = ttc::getTermBuilder(solver);
  cvc5::Term formula;
  if (assertions.empty())
  {
    formula = tm.mkBoolean(true);
  }
  else if (assertions.size() == 1)
  {
    formula = assertions[0];
  }
  else
  {
    formula = tm.mkTerm(cvc5::Kind::AND, assertions);
  }

  ttc::BitblastedCnfData data =
      ttc::getBitblastedCnf(solver, formula, projectionVars);

  // Count clauses (0-separated) in the flat DIMACS vector.
  std::size_t numClauses = 0;
  for (uint32_t v : data.clauses)
  {
    if (v == 0)
    {
      ++numClauses;
    }
  }

  std::ofstream out(path);
  if (!out)
  {
    throw std::runtime_error("Unable to open CNF output file: " + path);
  }
  out << "p cnf " << data.numVars << ' ' << numClauses << '\n';

  // Independent support / sampling set for projected counting. CryptoMiniSat
  // and ApproxMC accept multiple "c ind ... 0" lines; keep each line short.
  constexpr std::size_t kPerLine = 10;
  for (std::size_t i = 0; i < data.samplingVars.size(); i += kPerLine)
  {
    out << "c ind";
    for (std::size_t j = i; j < i + kPerLine && j < data.samplingVars.size();
         ++j)
    {
      out << ' ' << data.samplingVars[j];
    }
    out << " 0\n";
  }

  for (uint32_t v : data.clauses)
  {
    if (v == 0)
    {
      out << "0\n";
    }
    else
    {
      out << static_cast<int32_t>(v) << ' ';
    }
  }
  out.flush();

  BvCnfResult res;
  res.path = path;
  res.numVars = data.numVars;
  res.numClauses = numClauses;
  res.numSamplingVars = data.samplingVars.size();
  return res;
}

std::optional<std::string> runApproxMc(const std::string& cnfPath,
                                       std::uint64_t seed,
                                       int verbosity)
{
  std::string binary = locateApproxMc();
  std::ostringstream cmd;
  cmd << '"' << binary << '"' << " --seed " << seed << " --verb "
      << (verbosity > 1 ? 1 : 0) << " " << '"' << cnfPath << '"' << " 2>&1";

  Log(2) << "Running ApproxMC: " << cmd.str() << std::endl;

  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.str().c_str(), "r"),
                                                pclose);
  if (!pipe)
  {
    Log(0) << "Error: failed to launch approxmc" << std::endl;
    return std::nullopt;
  }

  std::optional<std::string> count;
  std::array<char, 4096> buf{};
  std::string line;
  while (fgets(buf.data(), static_cast<int>(buf.size()), pipe.get()) != nullptr)
  {
    line = buf.data();
    if (verbosity > 1)
    {
      std::cout << line;
    }
    // Parse "s mc <count>".
    std::istringstream iss(line);
    std::string s, mc;
    if ((iss >> s >> mc) && s == "s" && mc == "mc")
    {
      std::string value;
      if (iss >> value)
      {
        count = value;
      }
    }
  }
  return count;
}

}  // namespace eager
}  // namespace ttc
