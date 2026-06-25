#pragma once

#include <cvc5/cvc5.h>

#include <cstdint>
#include <iosfwd>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include "parser.hpp"

namespace ttc::tdd
{
// A Theory Decision Diagram (TDD): a reduced, ordered binary decision diagram
// whose decision nodes are the Boolean projection variables and whose leaves
// carry the residual theory (LRA) constraint reached along a path. A leaf whose
// constraint is theory-unsatisfiable is the FALSE terminal; every other leaf is
// a feasible region (it "contains the range of each real variable" as a
// conjunction of LRA atoms).
//
// The diagram is built bottom-up: each assertion is compiled into a small TDD
// and the assertions are combined with a classic apply(AND, .) operation that
// merges two smaller nodes into a larger one, hash-consing internal nodes and
// theory leaves as it goes. The (weighted) model count is the sum over every
// root->leaf path that does not reach FALSE of the product of the per-literal
// weights, with leaf value 1.
struct TddResult
{
  boost::multiprecision::cpp_int modelCount = 0;  // unweighted projected count
  long double weightedCount = 0.0L;               // == modelCount when unweighted
  bool hasWeights = false;
  std::size_t numNodes = 0;       // reduced internal (decision) nodes
  std::size_t numLeaves = 0;      // distinct feasible theory leaves
  std::size_t numVars = 0;        // Boolean decision variables
  std::uint64_t smtCalls = 0;     // theory feasibility checks issued to cvc5
};

class TheoryDD
{
 public:
  // The solver must already carry the parsed assertions (it is used, via
  // push/assert/checkSat/pop, only to test feasibility of theory leaves).
  TheoryDD(cvc5::Solver& solver,
           std::vector<cvc5::Term> boolVars,
           std::vector<cvc5::Term> realVars,
           std::unordered_map<cvc5::Term, TTCParser::LiteralWeight> weights,
           bool hasWeights);

  // Bottom-up compile the conjunction of the assertions into the diagram.
  void compile(const std::vector<cvc5::Term>& assertions);

  TddResult result() const;

  // Emit the reduced diagram in Graphviz dot format.
  void writeDot(std::ostream& os) const;

  std::uint64_t smtCalls() const { return d_smtCalls; }

 private:
  // Terminal ids are reserved: 0 = FALSE, 1 = TRUE.
  static constexpr int kFalse = 0;
  static constexpr int kTrue = 1;

  struct Node
  {
    bool leaf = false;          // true: theory leaf; false: decision node
    int var = -1;               // decision-variable index (decision nodes)
    int lo = kFalse;            // var = false child
    int hi = kFalse;            // var = true child
    cvc5::Term constraint;      // residual theory constraint (leaf nodes)
  };

  // --- diagram construction -------------------------------------------------
  int compileTerm(const cvc5::Term& term);
  int makeNode(int var, int lo, int hi);
  int makeLeaf(const cvc5::Term& constraint);
  int negate(int node);
  int apply(cvc5::Kind op, int f, int g);  // op in {AND, OR}

  // --- theory leaf helpers --------------------------------------------------
  bool leafSat(const cvc5::Term& constraint);
  std::string canonKey(const cvc5::Term& term) const;
  void collectAssoc(const cvc5::Term& term, cvc5::Kind op,
                    std::vector<std::string>& keys) const;

  // --- term inspection ------------------------------------------------------
  int projIndex(const cvc5::Term& term) const;  // -1 if not a projection var
  bool containsProj(const cvc5::Term& term) const;
  int minProjIndex(const cvc5::Term& term) const;  // INT_MAX if none

  // --- counting -------------------------------------------------------------
  int levelOf(int node) const;  // decision level; numVars for terminals/leaves
  long double factor(int varIdx) const;        // wpos + wneg
  long double wpos(int varIdx) const;
  long double wneg(int varIdx) const;

  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_boolVars;
  std::vector<cvc5::Term> d_realVars;
  std::unordered_map<cvc5::Term, TTCParser::LiteralWeight> d_weights;
  std::unordered_map<cvc5::Term, int> d_varIndex;
  bool d_hasWeights = false;

  std::vector<Node> d_nodes;  // index 0 = FALSE, 1 = TRUE
  int d_root = kFalse;

  std::map<std::tuple<int, int, int>, int> d_nodeTable;       // (var,lo,hi)->id
  std::unordered_map<std::string, int> d_leafTable;           // canonKey->id
  std::unordered_map<std::string, bool> d_satCache;           // canonKey->sat
  std::map<std::tuple<int, int, int>, int> d_applyMemo;       // (op,f,g)->id
  std::unordered_map<int, int> d_negMemo;                     // node->!node
  mutable std::unordered_map<std::string, int> d_projCache;   // term->minProj

  std::uint64_t d_smtCalls = 0;
};

}  // namespace ttc::tdd
