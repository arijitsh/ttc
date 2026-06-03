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

}  // namespace eager
}  // namespace ttc

#endif  // TTC_EAGER_BVCNF_HPP
