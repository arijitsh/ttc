#ifndef TTC_EAGER_BVCNF_HPP
#define TTC_EAGER_BVCNF_HPP

#include <cvc5/cvc5.h>

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace ttc
{
namespace eager
{

// Statistics describing the model-preserving bit-blasted CNF that was written.
struct BvCnfResult
{
  std::string path;             // file the CNF was written to
  std::uint32_t numVars = 0;    // number of CNF variables
  std::size_t numClauses = 0;   // number of clauses
  std::size_t numSamplingVars = 0;  // size of the projection / sampling set
};

// Bit-blast the (conjunction of the) given assertions into a model-preserving
// CNF and write it to `path` in DIMACS form. The bits of every bit-vector
// variable in `projectionVars` are emitted as the independent-support
// ("c ind ... 0") sampling set so that a projected model counter run on the
// file counts exactly the bit-vector models of the formula.
BvCnfResult writeBvCnf(cvc5::Solver& solver,
                       const std::vector<cvc5::Term>& assertions,
                       const std::vector<cvc5::Term>& projectionVars,
                       const std::string& path);

// Run ApproxMC on the given CNF file and return the reported model count
// (the value printed on the "s mc <count>" line). Returns std::nullopt if
// ApproxMC could not be run or produced no count.
std::optional<std::string> runApproxMc(const std::string& cnfPath,
                                       std::uint64_t seed,
                                       int verbosity);

// Statistics from the in-process eager bit-vector counting path.
struct BvCountResult
{
  std::string count;                // decimal model count (exact GMP integer)
  std::uint32_t numVars = 0;        // number of CNF variables (pre-Arjun)
  std::size_t numClauses = 0;       // number of clauses (pre-Arjun)
  std::size_t numSamplingVars = 0;  // size of the projection / sampling set
  std::uint32_t simplifiedVars = 0; // CNF variables after Arjun simplification
  std::size_t simplifiedSamplingVars = 0;  // sampling vars after Arjun
  double bitblastSeconds = 0.0;     // time spent bit-blasting to CNF
};

// Bit-blast the (conjunction of the) given assertions into a model-preserving
// CNF in memory, minimize the projection support with Arjun, and count the
// projected models with ApproxMC -- all in-process, with no intermediate CNF
// file. The bits of every bit-vector variable in `projectionVars` form the
// independent-support sampling set. When `dumpCnfPath` is non-null the same
// bit-blasted CNF is also written to that file in DIMACS form (the `--cnf`
// option), but the count itself is always produced from the in-memory CNF.
// Returns std::nullopt if the count could not be produced.
std::optional<BvCountResult> countBvModels(
    cvc5::Solver& solver,
    const std::vector<cvc5::Term>& assertions,
    const std::vector<cvc5::Term>& projectionVars,
    std::uint64_t seed,
    double epsilon,
    double delta,
    int verbosity,
    const std::string* dumpCnfPath);

}  // namespace eager
}  // namespace ttc

#endif  // TTC_EAGER_BVCNF_HPP
