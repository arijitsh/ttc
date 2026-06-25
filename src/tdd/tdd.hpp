#pragma once

#include <cvc5/cvc5.h>

#include <cstdint>
#include <functional>
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
  std::uint64_t feasiblePaths = 0;  // complete feasible assignments (Shannon)
};

class TheoryDD
{
 public:
  // How the diagram is built:
  //  - Shannon  : theory-pruned Shannon expansion. The whole formula stays on
  //               the solver as a feasibility oracle; an infeasible Boolean
  //               prefix prunes its entire subtree immediately, and the
  //               terminals are just TRUE/FALSE (a reduced feasibility BDD).
  //               Scales to many decision variables -- the default.
  //  - Regions  : bottom-up apply(AND, .) of the assertions, carrying the
  //               residual LRA region in each leaf (so leaves expose the range
  //               of each real variable). Richer for --printdd but does not
  //               scale: the per-path regions never merge.
  //  - Planned  : bottom-up apply(AND, .) like Regions, but a planning phase
  //               first chooses the decision-variable order and the order in
  //               which assertions are conjoined so as to minimize the number
  //               of theory checks (i.e. distinct residual regions). Bottom
  //               atom nodes are built first, then combined one assertion at a
  //               time in the planned schedule.
  enum class Mode
  {
    Shannon,
    Regions,
    Planned
  };

  // One candidate variable order considered by the planning phase, scored by
  // the induced width of its elimination order on the variable-interaction
  // graph (a proxy for the number of theory checks the build will issue).
  struct PlanCandidate
  {
    std::string name;
    int width = 0;
    bool usable = true;
  };

  // The solver must already carry the parsed assertions (it is used, via
  // push/assert/checkSat/pop, only to test feasibility).
  TheoryDD(cvc5::Solver& solver,
           std::vector<cvc5::Term> boolVars,
           std::vector<cvc5::Term> realVars,
           std::unordered_map<cvc5::Term, TTCParser::LiteralWeight> weights,
           bool hasWeights,
           Mode mode = Mode::Shannon);

  // Live progress, reported periodically during compile() so callers can render
  // a progress table in the style of the other engines.
  struct Progress
  {
    int level = 0;             // current decision depth being explored
    std::size_t nodes = 0;     // diagram nodes allocated so far
    std::uint64_t checks = 0;  // theory feasibility checks issued so far
    std::uint64_t paths = 0;   // complete feasible assignments found so far
  };

  // Run only the planning phase (Planned mode), so callers can report the chosen
  // order before the (possibly long) build. compile() reuses it if already run.
  void runPlanning(const std::vector<cvc5::Term>& assertions);

  // Bottom-up compile the conjunction of the assertions into the diagram.
  void compile(const std::vector<cvc5::Term>& assertions);

  TddResult result() const;

  // Emit the reduced diagram in Graphviz dot format.
  void writeDot(std::ostream& os) const;

  std::uint64_t smtCalls() const { return d_smtCalls; }

  // Planning-phase report (Planned mode): the candidate orders considered, the
  // index of the chosen one, and the resulting decision-variable order. Empty
  // until compile() runs in Planned mode.
  const std::vector<PlanCandidate>& planCandidates() const
  {
    return d_planCandidates;
  }
  int planChosen() const { return d_planChosen; }
  std::vector<std::string> decisionOrderNames() const;

  // Register a callback fired repeatedly during compilation. The callback is
  // expected to rate-limit its own output (e.g. by elapsed time).
  void setProgressCallback(std::function<void(const Progress&)> cb)
  {
    d_progress = std::move(cb);
  }

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

  // --- planning phase (Planned mode) ----------------------------------------
  // Build the variable-interaction graph, score candidate elimination orders by
  // induced width, pick the cheapest, and derive the decision-variable order
  // (reorders d_boolVars / d_varIndex) and the assertion apply schedule.
  void plan(const std::vector<cvc5::Term>& assertions);
  // Collect the tracked variables (Boolean projvars + real vars) of a term as
  // interaction-graph vertex ids.
  void collectVars(const cvc5::Term& term, std::vector<int>& out) const;
  // Greedy elimination order over the interaction graph (min-fill if byFill,
  // else min-degree); returns the order and sets widthOut to its induced width.
  std::vector<int> greedyOrder(const std::vector<std::vector<int>>& adj,
                               bool byFill, int& widthOut) const;
  // Induced width of an arbitrary elimination order on the interaction graph.
  int inducedWidth(const std::vector<std::vector<int>>& adj,
                   const std::vector<int>& order) const;

  // --- diagram construction -------------------------------------------------
  int compilePlanned(const std::vector<cvc5::Term>& assertions);
  int buildShannon(int level);             // theory-pruned Shannon expansion
  int compileTerm(const cvc5::Term& term);  // Regions mode (apply leaves)
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
  Mode d_mode = Mode::Shannon;

  std::vector<Node> d_nodes;  // index 0 = FALSE, 1 = TRUE
  int d_root = kFalse;

  std::map<std::tuple<int, int, int>, int> d_nodeTable;       // (var,lo,hi)->id
  std::unordered_map<std::string, int> d_leafTable;           // canonKey->id
  std::unordered_map<std::string, bool> d_satCache;           // canonKey->sat
  std::map<std::tuple<int, int, int>, int> d_applyMemo;       // (op,f,g)->id
  std::unordered_map<int, int> d_negMemo;                     // node->!node
  mutable std::unordered_map<std::string, int> d_projCache;   // term->minProj

  std::uint64_t d_smtCalls = 0;
  std::uint64_t d_paths = 0;  // complete feasible assignments found (Shannon)
  std::uint64_t d_leafCount = 0;  // distinct feasible leaves materialized
  std::function<void(const Progress&)> d_progress;

  // Planning-phase results (Planned mode).
  std::unordered_map<cvc5::Term, int> d_vertexOf;  // var term -> graph vertex id
  std::vector<PlanCandidate> d_planCandidates;
  int d_planChosen = -1;
  std::vector<std::size_t> d_applySchedule;  // assertion indices, apply order
  bool d_planned = false;
};

}  // namespace ttc::tdd
