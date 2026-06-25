#include "tdd/tdd.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <ostream>
#include <stdexcept>

#include "features.hpp"

namespace ttc::tdd
{
using boost::multiprecision::cpp_int;

TheoryDD::TheoryDD(cvc5::Solver& solver,
                   std::vector<cvc5::Term> boolVars,
                   std::vector<cvc5::Term> realVars,
                   std::unordered_map<cvc5::Term, TTCParser::LiteralWeight> weights,
                   bool hasWeights)
    : d_solver(solver),
      d_boolVars(std::move(boolVars)),
      d_realVars(std::move(realVars)),
      d_weights(std::move(weights)),
      d_hasWeights(hasWeights)
{
  for (std::size_t i = 0; i < d_boolVars.size(); ++i)
  {
    d_varIndex.emplace(d_boolVars[i], static_cast<int>(i));
  }
  auto& tm = ttc::getTermBuilder(d_solver);
  // Reserve the two terminals. FALSE carries the empty (unsatisfiable) region;
  // TRUE carries the trivially-true region.
  d_nodes.resize(2);
  d_nodes[kFalse] = Node{true, -1, kFalse, kFalse, tm.mkBoolean(false)};
  d_nodes[kTrue] = Node{true, -1, kFalse, kFalse, tm.mkBoolean(true)};
}

// ---------------------------------------------------------------------------
// Term inspection
// ---------------------------------------------------------------------------
int TheoryDD::projIndex(const cvc5::Term& term) const
{
  auto it = d_varIndex.find(term);
  return it == d_varIndex.end() ? -1 : it->second;
}

bool TheoryDD::containsProj(const cvc5::Term& term) const
{
  return minProjIndex(term) != INT_MAX;
}

int TheoryDD::minProjIndex(const cvc5::Term& term) const
{
  int direct = projIndex(term);
  if (direct >= 0)
  {
    return direct;
  }
  const std::string key = term.toString();
  auto cached = d_projCache.find(key);
  if (cached != d_projCache.end())
  {
    return cached->second;
  }
  int best = INT_MAX;
  for (std::size_t i = 0, n = term.getNumChildren(); i < n; ++i)
  {
    best = std::min(best, minProjIndex(term[i]));
  }
  d_projCache.emplace(key, best);
  return best;
}

// ---------------------------------------------------------------------------
// Theory-leaf bookkeeping
// ---------------------------------------------------------------------------
void TheoryDD::collectAssoc(const cvc5::Term& term, cvc5::Kind op,
                            std::vector<std::string>& keys) const
{
  if (term.getKind() == op)
  {
    for (std::size_t i = 0, n = term.getNumChildren(); i < n; ++i)
    {
      collectAssoc(term[i], op, keys);
    }
  }
  else
  {
    keys.push_back(canonKey(term));
  }
}

// A canonical string for a theory constraint so that conjunctions/disjunctions
// that differ only in operand order (or that repeat operands) hash-cons to the
// same leaf. Atoms below an AND/OR are flattened, canonicalized, sorted and
// de-duplicated; everything else falls back to cvc5's term printer.
std::string TheoryDD::canonKey(const cvc5::Term& term) const
{
  cvc5::Kind k = term.getKind();
  if (k == cvc5::Kind::AND || k == cvc5::Kind::OR)
  {
    std::vector<std::string> keys;
    collectAssoc(term, k, keys);
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    std::string out = (k == cvc5::Kind::AND) ? "(and" : "(or";
    for (const std::string& s : keys)
    {
      out += ' ';
      out += s;
    }
    out += ")";
    return out;
  }
  return term.toString();
}

bool TheoryDD::leafSat(const cvc5::Term& constraint)
{
  const std::string key = canonKey(constraint);
  auto it = d_satCache.find(key);
  if (it != d_satCache.end())
  {
    return it->second;
  }
  d_solver.push();
  d_solver.assertFormula(constraint);
  cvc5::Result res = d_solver.checkSat();
  d_solver.pop();
  ++d_smtCalls;
  if (!res.isSat() && !res.isUnsat())
  {
    throw std::runtime_error("tdd: solver returned unknown on a theory leaf");
  }
  bool sat = res.isSat();
  d_satCache.emplace(key, sat);
  return sat;
}

// ---------------------------------------------------------------------------
// Node construction (hash-consing + reduction)
// ---------------------------------------------------------------------------
int TheoryDD::makeLeaf(const cvc5::Term& constraint)
{
  if (constraint.getKind() == cvc5::Kind::CONST_BOOLEAN)
  {
    return constraint.getBooleanValue() ? kTrue : kFalse;
  }
  const std::string key = canonKey(constraint);
  auto it = d_leafTable.find(key);
  if (it != d_leafTable.end())
  {
    return it->second;
  }
  // Theory pruning: an infeasible residual constraint collapses to FALSE.
  if (!leafSat(constraint))
  {
    d_leafTable.emplace(key, kFalse);
    return kFalse;
  }
  int id = static_cast<int>(d_nodes.size());
  d_nodes.push_back(Node{true, -1, kFalse, kFalse, constraint});
  d_leafTable.emplace(key, id);
  return id;
}

int TheoryDD::makeNode(int var, int lo, int hi)
{
  if (lo == hi)
  {
    return lo;  // redundant-node elimination
  }
  auto key = std::make_tuple(var, lo, hi);
  auto it = d_nodeTable.find(key);
  if (it != d_nodeTable.end())
  {
    return it->second;
  }
  int id = static_cast<int>(d_nodes.size());
  d_nodes.push_back(Node{false, var, lo, hi, cvc5::Term()});
  d_nodeTable.emplace(key, id);
  return id;
}

// ---------------------------------------------------------------------------
// Boolean operators over diagrams
// ---------------------------------------------------------------------------
int TheoryDD::negate(int node)
{
  if (node == kFalse) return kTrue;
  if (node == kTrue) return kFalse;
  auto it = d_negMemo.find(node);
  if (it != d_negMemo.end())
  {
    return it->second;
  }
  const Node n = d_nodes[node];
  int result;
  if (n.leaf)
  {
    auto& tm = ttc::getTermBuilder(d_solver);
    result = makeLeaf(tm.mkTerm(cvc5::Kind::NOT, {n.constraint}));
  }
  else
  {
    result = makeNode(n.var, negate(n.lo), negate(n.hi));
  }
  d_negMemo.emplace(node, result);
  return result;
}

// Classic apply: recursively combine two smaller diagrams into a larger one,
// branching on the smaller of the two top decision variables and combining the
// theory leaves once both operands have bottomed out.
int TheoryDD::apply(cvc5::Kind op, int f, int g)
{
  // Terminal short-circuits.
  if (op == cvc5::Kind::AND)
  {
    if (f == kFalse || g == kFalse) return kFalse;
    if (f == kTrue) return g;
    if (g == kTrue) return f;
  }
  else  // OR
  {
    if (f == kTrue || g == kTrue) return kTrue;
    if (f == kFalse) return g;
    if (g == kFalse) return f;
  }
  if (f == g) return f;

  auto memoKey = std::make_tuple(static_cast<int>(op), std::min(f, g),
                                 std::max(f, g));
  auto memoIt = d_applyMemo.find(memoKey);
  if (memoIt != d_applyMemo.end())
  {
    return memoIt->second;
  }

  const Node nf = d_nodes[f];
  const Node ng = d_nodes[g];

  int result;
  if (nf.leaf && ng.leaf)
  {
    // Combine the residual theory constraints of two leaves.
    auto& tm = ttc::getTermBuilder(d_solver);
    cvc5::Term combined = tm.mkTerm(op, {nf.constraint, ng.constraint});
    result = makeLeaf(combined);
  }
  else
  {
    int vf = nf.leaf ? INT_MAX : nf.var;
    int vg = ng.leaf ? INT_MAX : ng.var;
    int v = std::min(vf, vg);
    int flo = (vf == v) ? nf.lo : f;
    int fhi = (vf == v) ? nf.hi : f;
    int glo = (vg == v) ? ng.lo : g;
    int ghi = (vg == v) ? ng.hi : g;
    int lo = apply(op, flo, glo);
    int hi = apply(op, fhi, ghi);
    result = makeNode(v, lo, hi);
  }
  d_applyMemo.emplace(memoKey, result);
  return result;
}

// ---------------------------------------------------------------------------
// Compile one assertion term into a diagram (bottom-up via apply)
// ---------------------------------------------------------------------------
int TheoryDD::compileTerm(const cvc5::Term& term)
{
  int idx = projIndex(term);
  if (idx >= 0)
  {
    // A bare Boolean decision variable: node(var, FALSE, TRUE).
    return makeNode(idx, kFalse, kTrue);
  }

  auto& tm = ttc::getTermBuilder(d_solver);
  cvc5::Kind k = term.getKind();

  if (k == cvc5::Kind::CONST_BOOLEAN)
  {
    return term.getBooleanValue() ? kTrue : kFalse;
  }
  if (k == cvc5::Kind::NOT)
  {
    return negate(compileTerm(term[0]));
  }
  if (k == cvc5::Kind::AND)
  {
    int acc = kTrue;
    for (std::size_t i = 0, n = term.getNumChildren(); i < n; ++i)
    {
      acc = apply(cvc5::Kind::AND, acc, compileTerm(term[i]));
    }
    return acc;
  }
  if (k == cvc5::Kind::OR)
  {
    int acc = kFalse;
    for (std::size_t i = 0, n = term.getNumChildren(); i < n; ++i)
    {
      acc = apply(cvc5::Kind::OR, acc, compileTerm(term[i]));
    }
    return acc;
  }
  if (k == cvc5::Kind::IMPLIES)
  {
    return apply(cvc5::Kind::OR, negate(compileTerm(term[0])),
                 compileTerm(term[1]));
  }
  // Boolean equivalence / xor over Boolean operands.
  if ((k == cvc5::Kind::EQUAL || k == cvc5::Kind::XOR)
      && term.getNumChildren() == 2 && term[0].getSort().isBoolean())
  {
    int a = compileTerm(term[0]);
    int b = compileTerm(term[1]);
    int iff = apply(cvc5::Kind::OR, apply(cvc5::Kind::AND, a, b),
                    apply(cvc5::Kind::AND, negate(a), negate(b)));
    return k == cvc5::Kind::EQUAL ? iff : negate(iff);
  }
  // Boolean if-then-else.
  if (k == cvc5::Kind::ITE && term.getSort().isBoolean())
  {
    int c = compileTerm(term[0]);
    int t = compileTerm(term[1]);
    int e = compileTerm(term[2]);
    return apply(cvc5::Kind::OR, apply(cvc5::Kind::AND, c, t),
                 apply(cvc5::Kind::AND, negate(c), e));
  }

  // Default: a theory atom (or any sub-term with no Boolean structure left).
  // If it still mentions a decision variable, Shannon-expand on the smallest
  // such variable so the diagram keeps its decision order; otherwise it is a
  // pure theory leaf.
  if (!containsProj(term))
  {
    return makeLeaf(term);
  }
  int sv = minProjIndex(term);
  cvc5::Term var = d_boolVars[sv];
  int hi = compileTerm(term.substitute({var}, {tm.mkBoolean(true)}));
  int lo = compileTerm(term.substitute({var}, {tm.mkBoolean(false)}));
  return makeNode(sv, lo, hi);
}

void TheoryDD::compile(const std::vector<cvc5::Term>& assertions)
{
  int acc = kTrue;
  for (const cvc5::Term& a : assertions)
  {
    acc = apply(cvc5::Kind::AND, acc, compileTerm(a));
    if (acc == kFalse)
    {
      break;  // whole formula is theory-unsatisfiable
    }
  }
  d_root = acc;
}

// ---------------------------------------------------------------------------
// Weighted model counting over the reduced diagram
// ---------------------------------------------------------------------------
int TheoryDD::levelOf(int node) const
{
  if (node <= kTrue || d_nodes[node].leaf)
  {
    return static_cast<int>(d_boolVars.size());
  }
  return d_nodes[node].var;
}

long double TheoryDD::wpos(int varIdx) const
{
  if (!d_hasWeights) return 1.0L;
  auto it = d_weights.find(d_boolVars[varIdx]);
  const TTCParser::LiteralWeight def;
  return static_cast<long double>(it == d_weights.end() ? def.positive
                                                        : it->second.positive);
}

long double TheoryDD::wneg(int varIdx) const
{
  if (!d_hasWeights) return 1.0L;
  auto it = d_weights.find(d_boolVars[varIdx]);
  const TTCParser::LiteralWeight def;
  return static_cast<long double>(it == d_weights.end() ? def.negative
                                                        : it->second.negative);
}

long double TheoryDD::factor(int varIdx) const
{
  return wpos(varIdx) + wneg(varIdx);
}

TddResult TheoryDD::result() const
{
  TddResult r;
  r.numVars = d_boolVars.size();
  r.hasWeights = d_hasWeights;
  r.smtCalls = d_smtCalls;
  std::size_t nodes = 0;
  std::size_t leaves = 0;
  for (std::size_t id = 2; id < d_nodes.size(); ++id)
  {
    if (d_nodes[id].leaf)
      ++leaves;
    else
      ++nodes;
  }
  r.numNodes = nodes;
  r.numLeaves = leaves;

  const int n = static_cast<int>(d_boolVars.size());

  // Exact unweighted projected model count (cpp_int). A skipped variable is
  // free, contributing a factor of two.
  std::vector<cpp_int> mc(d_nodes.size(), cpp_int(0));
  mc[kTrue] = 1;
  auto skipPow2 = [&](int from, int to) {
    cpp_int p = 1;
    for (int i = from; i < to; ++i) p <<= 1;
    return p;
  };
  for (std::size_t id = 2; id < d_nodes.size(); ++id)
  {
    const Node& nd = d_nodes[id];
    if (nd.leaf)
    {
      mc[id] = 1;  // every retained leaf is feasible
      continue;
    }
    int L = nd.var;
    mc[id] = skipPow2(L + 1, levelOf(nd.lo)) * mc[nd.lo]
             + skipPow2(L + 1, levelOf(nd.hi)) * mc[nd.hi];
  }
  r.modelCount = skipPow2(0, levelOf(d_root)) * mc[d_root];

  if (d_hasWeights)
  {
    std::vector<long double> w(d_nodes.size(), 0.0L);
    w[kTrue] = 1.0L;
    auto skip = [&](int from, int to) {
      long double p = 1.0L;
      for (int i = from; i < to; ++i) p *= factor(i);
      return p;
    };
    for (std::size_t id = 2; id < d_nodes.size(); ++id)
    {
      const Node& nd = d_nodes[id];
      if (nd.leaf)
      {
        w[id] = 1.0L;
        continue;
      }
      int L = nd.var;
      w[id] = wneg(L) * skip(L + 1, levelOf(nd.lo)) * w[nd.lo]
              + wpos(L) * skip(L + 1, levelOf(nd.hi)) * w[nd.hi];
    }
    r.weightedCount = skip(0, levelOf(d_root)) * w[d_root];
  }
  else
  {
    r.weightedCount = static_cast<long double>(r.modelCount);
  }
  (void)n;
  return r;
}

// ---------------------------------------------------------------------------
// Graphviz output
// ---------------------------------------------------------------------------
namespace
{
std::string dotEscape(const std::string& s)
{
  std::string out;
  out.reserve(s.size() + 8);
  for (char c : s)
  {
    if (c == '"' || c == '\\') out.push_back('\\');
    out.push_back(c);
  }
  return out;
}
}  // namespace

void TheoryDD::writeDot(std::ostream& os) const
{
  os << "digraph tdd {\n";
  os << "  rankdir=TB;\n";
  os << "  node [fontname=\"monospace\"];\n";
  os << "  n" << kFalse << " [shape=box,label=\"FALSE\"];\n";
  os << "  n" << kTrue << " [shape=box,label=\"TRUE\"];\n";

  bool falseReachable = (d_root == kFalse);
  for (std::size_t id = 2; id < d_nodes.size(); ++id)
  {
    const Node& nd = d_nodes[id];
    if (nd.leaf)
    {
      // The leaf "contains the range of each real variable" as the conjunction
      // of LRA atoms defining its feasible region.
      os << "  n" << id << " [shape=box,label=\""
         << dotEscape(canonKey(nd.constraint)) << "\"];\n";
    }
    else
    {
      os << "  n" << id << " [shape=ellipse,label=\""
         << dotEscape(d_boolVars[nd.var].toString()) << "\"];\n";
    }
  }
  for (std::size_t id = 2; id < d_nodes.size(); ++id)
  {
    const Node& nd = d_nodes[id];
    if (nd.leaf) continue;
    if (nd.lo == kFalse || nd.hi == kFalse) falseReachable = true;
    // low (variable = false): dashed; high (variable = true): solid.
    os << "  n" << id << " -> n" << nd.lo << " [style=dashed];\n";
    os << "  n" << id << " -> n" << nd.hi << " [style=solid];\n";
  }
  os << "  root [shape=point];\n";
  os << "  root -> n" << d_root << ";\n";
  (void)falseReachable;
  os << "}\n";
}

}  // namespace ttc::tdd
