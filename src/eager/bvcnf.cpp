#include "eager/bvcnf.hpp"

#include "features.hpp"
#include "logger.hpp"

#include <approxmc/approxmc.h>
#include <arjun/arjun.h>

#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <gmpxx.h>
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

// Count clauses (0-separated) in the flat DIMACS literal vector.
std::size_t countClauses(const std::vector<uint32_t>& clauses)
{
  std::size_t n = 0;
  for (uint32_t v : clauses)
  {
    if (v == 0)
    {
      ++n;
    }
  }
  return n;
}

// Write a bit-blasted CNF (flat 0-separated literals + 1-based sampling set) to
// `path` in DIMACS form, emitting the sampling set as "c ind ... 0" lines.
void writeDimacs(const ttc::BitblastedCnfData& data, std::size_t numClauses,
                 const std::string& path)
{
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

// Conjoin the assertions into a single formula for bit-blasting.
cvc5::Term buildFormula(cvc5::Solver& solver,
                        const std::vector<cvc5::Term>& assertions)
{
  auto& tm = ttc::getTermBuilder(solver);
  if (assertions.empty())
  {
    return tm.mkBoolean(true);
  }
  if (assertions.size() == 1)
  {
    return assertions[0];
  }
  return tm.mkTerm(cvc5::Kind::AND, assertions);
}

}  // namespace

BvCnfResult writeBvCnf(cvc5::Solver& solver,
                       const std::vector<cvc5::Term>& assertions,
                       const std::vector<cvc5::Term>& projectionVars,
                       const std::string& path)
{
  cvc5::Term formula = buildFormula(solver, assertions);
  ttc::BitblastedCnfData data =
      ttc::getBitblastedCnf(solver, formula, projectionVars);

  std::size_t numClauses = countClauses(data.clauses);
  writeDimacs(data, numClauses, path);

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

std::optional<BvCountResult> countBvModels(
    cvc5::Solver& solver,
    const std::vector<cvc5::Term>& assertions,
    const std::vector<cvc5::Term>& projectionVars,
    std::uint64_t seed,
    double epsilon,
    double delta,
    int verbosity,
    const std::string* dumpCnfPath)
{
  cvc5::Term formula = buildFormula(solver, assertions);

  double bbStart = Log.elapsed();
  ttc::BitblastedCnfData data =
      ttc::getBitblastedCnf(solver, formula, projectionVars);
  double bbEnd = Log.elapsed();

  BvCountResult res;
  res.numVars = data.numVars;
  res.numClauses = countClauses(data.clauses);
  res.numSamplingVars = data.samplingVars.size();
  res.bitblastSeconds = bbEnd - bbStart;

  // Optionally persist the same bit-blasted CNF (the --cnf option).
  if (dumpCnfPath != nullptr)
  {
    writeDimacs(data, res.numClauses, *dumpCnfPath);
  }

  // Build an Arjun SimplifiedCNF directly from the flat DIMACS literals, mirror
  // the CSB (~/solvers/csb ApxMC) integration: minimize the projection support
  // with Arjun, then count the projected models with ApproxMC -- all in-process
  // over an exact-GMP field, with no intermediate file round-trip.
  std::unique_ptr<CMSat::FieldGen> fg =
      std::make_unique<ArjunNS::FGenMpz>();
  ApproxMC::AppMC appmc(fg);
  ArjunNS::SimplifiedCNF cnf(fg);

  cnf.new_vars(data.numVars);

  std::vector<CMSat::Lit> clause;
  for (uint32_t v : data.clauses)
  {
    if (v == 0)
    {
      cnf.add_clause(clause);
      clause.clear();
    }
    else
    {
      int32_t lit = static_cast<int32_t>(v);
      uint32_t var = static_cast<uint32_t>(std::abs(lit)) - 1;
      clause.emplace_back(var, lit < 0);
    }
  }
  // Tolerate a trailing clause not terminated by a 0 separator.
  if (!clause.empty())
  {
    cnf.add_clause(clause);
  }

  // Sampling vars come back 1-based from the bit-blaster; Arjun/ApproxMC use
  // 0-based variable indices.
  std::vector<uint32_t> samplVars;
  samplVars.reserve(data.samplingVars.size());
  for (uint32_t v : data.samplingVars)
  {
    samplVars.push_back(v - 1);
  }
  cnf.set_weighted(false);
  cnf.set_sampl_vars(samplVars);

  ArjunNS::Arjun arjun;
  arjun.set_verb(verbosity > 1 ? verbosity : 0);

  // Arjun preprocessing tuned for ApproxMC, matching the CSB ApxMC settings.
  ArjunNS::SimpConf sc;
  sc.appmc = true;
  sc.oracle_vivify = true;
  sc.oracle_vivify_get_learnts = true;
  sc.oracle_sparsify = false;
  sc.iter1 = 2;
  sc.iter2 = 0;

  ArjunNS::SimplifiedCNF simp = arjun.standalone_get_simplified_cnf(cnf, sc);

  res.simplifiedVars = simp.nvars;
  res.simplifiedSamplingVars = simp.sampl_vars.size();

  appmc.set_seed(static_cast<uint32_t>(seed));
  appmc.set_verbosity(verbosity > 1 ? 1 : 0);
  appmc.set_epsilon(epsilon);
  appmc.set_delta(delta);

  appmc.new_vars(simp.nvars);
  for (const auto& cl : simp.clauses)
  {
    appmc.add_clause(cl);
  }
  for (const auto& cl : simp.red_clauses)
  {
    appmc.add_clause(cl);
  }
  appmc.set_multiplier_weight(simp.multiplier_weight);
  appmc.set_sampl_vars(simp.sampl_vars);

  ApproxMC::SolCount solCount = appmc.count();

  // Absolute count = multiplier * cellSolCount * 2^hashCount, in exact GMP.
  mpz_class result;
  mpz_class cellSolCount(solCount.cellSolCount);
  mpz_mul_2exp(result.get_mpz_t(), cellSolCount.get_mpz_t(),
               solCount.hashCount);

  // Arjun folds unconstrained projection variables into multiplier_weight (e.g.
  // a free n-bit variable becomes a 2^n factor). With the unweighted FGenMpz
  // field this is an integer FMpz; the weighted path yields a rational FMpq.
  mpq_class multiplier(1);
  const CMSat::Field* ptr = appmc.get_multiplier_weight().get();
  if (ptr != nullptr)
  {
    if (const ArjunNS::FMpz* mult = dynamic_cast<const ArjunNS::FMpz*>(ptr))
    {
      multiplier = mpq_class(mult->val);
    }
    else if (const ArjunNS::FMpq* mult =
                 dynamic_cast<const ArjunNS::FMpq*>(ptr))
    {
      multiplier = mult->val;
    }
  }
  mpq_class finalCount = multiplier * result;
  finalCount.canonicalize();

  // The bit-vector model count is an integer; emit it as such.
  if (finalCount.get_den() == 1)
  {
    res.count = finalCount.get_num().get_str();
  }
  else
  {
    res.count = finalCount.get_str();
  }
  return res;
}

}  // namespace eager
}  // namespace ttc
