#include "tdd/tdd.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <unordered_set>

#include "features.hpp"
#if defined(TTC_ENABLE_DDNNF)
#include "var_order.hpp"
#endif

namespace ttc::tdd
{
using boost::multiprecision::cpp_int;

TheoryDD::TheoryDD(cvc5::Solver& solver,
                   std::vector<cvc5::Term> boolVars,
                   std::vector<cvc5::Term> realVars,
                   std::unordered_map<cvc5::Term, TTCParser::LiteralWeight> weights,
                   bool hasWeights,
                   Mode mode)
    : d_solver(solver),
      d_boolVars(std::move(boolVars)),
      d_realVars(std::move(realVars)),
      d_weights(std::move(weights)),
      d_hasWeights(hasWeights),
      d_mode(mode)
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
  ++d_leafCount;
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

// Theory-pruned Shannon expansion. The full formula is already asserted on the
// solver, so we only push the current decision literal and ask whether the
// prefix is still theory-satisfiable. An UNSAT child is the FALSE terminal,
// which prunes the entire subtree below it; makeNode then reduces the result
// into a feasibility BDD with TRUE/FALSE terminals.
//
// Precondition: the current solver context (formula + asserted prefix) is SAT.
int TheoryDD::buildShannon(int level)
{
  if (level == static_cast<int>(d_boolVars.size()))
  {
    ++d_paths;  // a complete, theory-feasible assignment
    return kTrue;
  }
  auto& tm = ttc::getTermBuilder(d_solver);
  const cvc5::Term& v = d_boolVars[level];

  if (d_progress)
  {
    d_progress(Progress{level, d_nodes.size(), d_smtCalls, d_paths});
  }

  d_solver.push();
  d_solver.assertFormula(v);
  ++d_smtCalls;
  int hi = d_solver.checkSat().isSat() ? buildShannon(level + 1) : kFalse;
  d_solver.pop();

  d_solver.push();
  d_solver.assertFormula(tm.mkTerm(cvc5::Kind::NOT, {v}));
  ++d_smtCalls;
  int lo = d_solver.checkSat().isSat() ? buildShannon(level + 1) : kFalse;
  d_solver.pop();

  return makeNode(level, lo, hi);
}

// ---------------------------------------------------------------------------
// Planning phase
// ---------------------------------------------------------------------------
void TheoryDD::collectVars(const cvc5::Term& term, std::vector<int>& out) const
{
  std::unordered_set<cvc5::Term> visited;
  std::vector<cvc5::Term> stack{term};
  std::unordered_set<int> seen;
  while (!stack.empty())
  {
    cvc5::Term t = stack.back();
    stack.pop_back();
    if (!visited.insert(t).second)
    {
      continue;
    }
    auto it = d_vertexOf.find(t);
    if (it != d_vertexOf.end() && seen.insert(it->second).second)
    {
      out.push_back(it->second);
    }
    for (std::size_t i = 0, n = t.getNumChildren(); i < n; ++i)
    {
      stack.push_back(t[i]);
    }
  }
}

int TheoryDD::inducedWidth(const std::vector<std::vector<int>>& adj,
                           const std::vector<int>& order) const
{
  const int n = static_cast<int>(adj.size());
  std::vector<std::unordered_set<int>> g(n);
  for (int u = 0; u < n; ++u)
  {
    g[u].insert(adj[u].begin(), adj[u].end());
  }
  std::vector<int> rank(n);
  for (int i = 0; i < static_cast<int>(order.size()); ++i)
  {
    rank[order[i]] = i;
  }
  int width = 0;
  for (int v : order)
  {
    // Neighbours not yet eliminated form this vertex's bag.
    std::vector<int> nb;
    for (int u : g[v])
    {
      if (rank[u] > rank[v])
      {
        nb.push_back(u);
      }
    }
    width = std::max(width, static_cast<int>(nb.size()));
    for (std::size_t i = 0; i < nb.size(); ++i)
    {
      for (std::size_t j = i + 1; j < nb.size(); ++j)
      {
        g[nb[i]].insert(nb[j]);
        g[nb[j]].insert(nb[i]);
      }
    }
  }
  return width;
}

std::vector<int> TheoryDD::greedyOrder(const std::vector<std::vector<int>>& adj,
                                       bool byFill, int& widthOut) const
{
  const int n = static_cast<int>(adj.size());
  std::vector<std::unordered_set<int>> g(n);
  for (int u = 0; u < n; ++u)
  {
    g[u].insert(adj[u].begin(), adj[u].end());
  }
  std::vector<bool> done(n, false);
  std::vector<int> order;
  order.reserve(n);
  widthOut = 0;
  for (int step = 0; step < n; ++step)
  {
    int best = -1;
    long bestCost = LONG_MAX;
    for (int v = 0; v < n; ++v)
    {
      if (done[v])
      {
        continue;
      }
      long cost;
      if (byFill)
      {
        // Number of non-adjacent neighbour pairs that eliminating v would link.
        std::vector<int> nb(g[v].begin(), g[v].end());
        long fill = 0;
        for (std::size_t i = 0; i < nb.size(); ++i)
        {
          for (std::size_t j = i + 1; j < nb.size(); ++j)
          {
            if (g[nb[i]].find(nb[j]) == g[nb[i]].end())
            {
              ++fill;
            }
          }
        }
        cost = fill;
      }
      else
      {
        cost = static_cast<long>(g[v].size());  // min-degree
      }
      if (cost < bestCost)
      {
        bestCost = cost;
        best = v;
      }
    }
    order.push_back(best);
    done[best] = true;
    std::vector<int> nb(g[best].begin(), g[best].end());
    widthOut = std::max(widthOut, static_cast<int>(nb.size()));
    for (std::size_t i = 0; i < nb.size(); ++i)
    {
      g[nb[i]].erase(best);
      for (std::size_t j = i + 1; j < nb.size(); ++j)
      {
        g[nb[i]].insert(nb[j]);
        g[nb[j]].insert(nb[i]);
      }
    }
    g[best].clear();
  }
  return order;
}

void TheoryDD::plan(const std::vector<cvc5::Term>& assertions)
{
  // Vertices: Boolean decision variables first, then real variables. The
  // Boolean order chosen here becomes the BDD decision order.
  const int numBool = static_cast<int>(d_boolVars.size());
  d_vertexOf.clear();
  std::vector<cvc5::Term> vertexTerm;
  for (int i = 0; i < numBool; ++i)
  {
    d_vertexOf.emplace(d_boolVars[i], static_cast<int>(vertexTerm.size()));
    vertexTerm.push_back(d_boolVars[i]);
  }
  for (const cvc5::Term& r : d_realVars)
  {
    if (d_vertexOf.emplace(r, static_cast<int>(vertexTerm.size())).second)
    {
      vertexTerm.push_back(r);
    }
  }
  const int n = static_cast<int>(vertexTerm.size());

  // Interaction graph: a clique over the variables of each assertion.
  std::vector<std::vector<int>> assertVars(assertions.size());
  std::vector<std::unordered_set<int>> adjSet(n);
  for (std::size_t a = 0; a < assertions.size(); ++a)
  {
    collectVars(assertions[a], assertVars[a]);
    const std::vector<int>& vs = assertVars[a];
    for (std::size_t i = 0; i < vs.size(); ++i)
    {
      for (std::size_t j = i + 1; j < vs.size(); ++j)
      {
        adjSet[vs[i]].insert(vs[j]);
        adjSet[vs[j]].insert(vs[i]);
      }
    }
  }
  std::vector<std::vector<int>> adj(n);
  for (int u = 0; u < n; ++u)
  {
    adj[u].assign(adjSet[u].begin(), adjSet[u].end());
  }

  // Candidate elimination orders, scored by induced width (a proxy for the
  // number of distinct residual regions, i.e. theory checks, the build issues).
  struct Cand
  {
    std::string name;
    std::vector<int> order;
    int width = 0;
    bool usable = true;
  };
  std::vector<Cand> cands;

  int wFill = 0;
  std::vector<int> oFill = greedyOrder(adj, /*byFill=*/true, wFill);
  cands.push_back({"min-fill", oFill, wFill, true});

  int wDeg = 0;
  std::vector<int> oDeg = greedyOrder(adj, /*byFill=*/false, wDeg);
  cands.push_back({"min-degree", oDeg, wDeg, true});

  std::vector<int> oDecl(n);
  std::iota(oDecl.begin(), oDecl.end(), 0);
  cands.push_back({"declaration", oDecl, inducedWidth(adj, oDecl), true});

#if defined(TTC_ENABLE_DDNNF)
  // Flowcutter tree-decomposition order over the Boolean projection variables.
  // Extend it to the full vertex set by placing each real variable right after
  // the last Boolean (in that order) it co-occurs with.
  try
  {
    auto& tm = ttc::getTermBuilder(d_solver);
    cvc5::Term formula = assertions.empty()
                             ? tm.mkBoolean(true)
                             : (assertions.size() == 1
                                    ? assertions[0]
                                    : tm.mkTerm(cvc5::Kind::AND, assertions));
    std::vector<cvc5::Term> ordered =
        computeProjVarOrder(formula, d_boolVars, d_solver, /*printTD=*/false);
    std::vector<int> boolRank(n, INT_MAX);
    int rk = 0;
    for (const cvc5::Term& t : ordered)
    {
      auto it = d_vertexOf.find(t);
      if (it != d_vertexOf.end())
      {
        boolRank[it->second] = rk++;
      }
    }
    // Key each vertex: Booleans by their flowcutter rank; reals by the max rank
    // of the Booleans they share an assertion with (so a real is eliminated just
    // after its last Boolean). Stable-sort vertices by key for the full order.
    std::vector<int> key(n, 0);
    for (int v = 0; v < numBool; ++v)
    {
      key[v] = boolRank[v] == INT_MAX ? rk : boolRank[v];
    }
    for (const auto& vs : assertVars)
    {
      int mx = 0;
      for (int v : vs)
      {
        if (v < numBool && boolRank[v] != INT_MAX)
        {
          mx = std::max(mx, boolRank[v]);
        }
      }
      for (int v : vs)
      {
        if (v >= numBool)
        {
          key[v] = std::max(key[v], mx);
        }
      }
    }
    std::vector<int> oFc(n);
    std::iota(oFc.begin(), oFc.end(), 0);
    std::stable_sort(oFc.begin(), oFc.end(),
                     [&](int a, int b) { return key[a] < key[b]; });
    cands.push_back({"flowcutter", oFc, inducedWidth(adj, oFc), true});
  }
  catch (const std::exception&)
  {
    // computeProjVarOrder may throw on degenerate inputs; just skip it.
  }
#endif

  // Choose the minimum-width candidate.
  int chosen = 0;
  for (int i = 1; i < static_cast<int>(cands.size()); ++i)
  {
    if (cands[i].width < cands[chosen].width)
    {
      chosen = i;
    }
  }
  d_planCandidates.clear();
  for (const Cand& c : cands)
  {
    d_planCandidates.push_back({c.name, c.width, c.usable});
  }
  d_planChosen = chosen;
  const std::vector<int>& elim = cands[chosen].order;

  // Decision order: Boolean variables sorted by their elimination rank. Reorder
  // d_boolVars / d_varIndex so the new index is the BDD level.
  std::vector<int> rank(n);
  for (int i = 0; i < n; ++i)
  {
    rank[elim[i]] = i;
  }
  std::vector<cvc5::Term> reordered = d_boolVars;
  std::stable_sort(reordered.begin(), reordered.end(),
                   [&](const cvc5::Term& a, const cvc5::Term& b) {
                     return rank[d_vertexOf.at(a)] < rank[d_vertexOf.at(b)];
                   });
  d_boolVars = std::move(reordered);
  d_varIndex.clear();
  for (std::size_t i = 0; i < d_boolVars.size(); ++i)
  {
    d_varIndex.emplace(d_boolVars[i], static_cast<int>(i));
  }

  // Apply schedule: conjoin each assertion once the last of its variables enters
  // play -- i.e. sorted by the maximum elimination rank among its variables
  // (bucket elimination). Empty assertions go first.
  d_applySchedule.resize(assertions.size());
  std::iota(d_applySchedule.begin(), d_applySchedule.end(), std::size_t{0});
  std::vector<int> assertKey(assertions.size(), -1);
  for (std::size_t a = 0; a < assertions.size(); ++a)
  {
    int mx = -1;
    for (int v : assertVars[a])
    {
      mx = std::max(mx, rank[v]);
    }
    assertKey[a] = mx;
  }
  std::stable_sort(d_applySchedule.begin(), d_applySchedule.end(),
                   [&](std::size_t a, std::size_t b) {
                     return assertKey[a] < assertKey[b];
                   });

  // Projection schedule (Stage 2): a real variable can be projected out of the
  // leaves once the last (in apply order) assertion mentioning it has been
  // applied. Record, per schedule step, the real variables eliminated after it.
  std::vector<int> schedPos(assertions.size(), 0);
  for (std::size_t s = 0; s < d_applySchedule.size(); ++s)
  {
    schedPos[d_applySchedule[s]] = static_cast<int>(s);
  }
  std::vector<int> elimStep(n, -1);
  for (std::size_t a = 0; a < assertions.size(); ++a)
  {
    for (int v : assertVars[a])
    {
      if (v >= numBool)  // real variable
      {
        elimStep[v] = std::max(elimStep[v], schedPos[a]);
      }
    }
  }
  d_projectAfter.assign(d_applySchedule.size(), {});
  for (int v = numBool; v < n; ++v)
  {
    if (elimStep[v] >= 0)
    {
      d_projectAfter[elimStep[v]].push_back(vertexTerm[v]);
    }
  }
  d_planned = true;
}

void TheoryDD::runPlanning(const std::vector<cvc5::Term>& assertions)
{
  if (!d_planned)
  {
    plan(assertions);
  }
}

std::vector<std::string> TheoryDD::decisionOrderNames() const
{
  std::vector<std::string> names;
  names.reserve(d_boolVars.size());
  for (const cvc5::Term& v : d_boolVars)
  {
    names.push_back(v.toString());
  }
  return names;
}

// A lazily-created helper solver with a quantifier-capable logic, used only for
// quantifier elimination on leaf regions. It shares the term manager with the
// main solver so leaf-constraint terms can be passed directly.
cvc5::Solver& TheoryDD::qeSolver()
{
  if (!d_qeSolver)
  {
    // Construct in place from the shared term manager (cvc5::Solver is not
    // movable, so it cannot be returned-by-value into the unique_ptr).
    d_qeSolver = std::make_unique<cvc5::Solver>(d_solver.getTermManager());
    try { d_qeSolver->setLogic("ALL"); }
    catch (const cvc5::CVC5ApiException&) {}
  }
  return *d_qeSolver;
}

// Project the given real variables out of a leaf region: return a quantifier-
// free constraint equivalent to (exists vars. constraint). On any failure the
// original constraint is returned (sound -- it just forgoes the merge).
cvc5::Term TheoryDD::projectConstraint(const cvc5::Term& constraint,
                                       const std::vector<cvc5::Term>& vars)
{
  if (vars.empty() || constraint.getKind() == cvc5::Kind::CONST_BOOLEAN)
  {
    return constraint;
  }
  auto& tm = ttc::getTermBuilder(d_solver);
  ++d_qeCalls;
  try
  {
    std::vector<cvc5::Term> bound;
    bound.reserve(vars.size());
    for (std::size_t i = 0; i < vars.size(); ++i)
    {
      bound.push_back(tm.mkVar(vars[i].getSort(), "__tdd_qe_" + std::to_string(i)));
    }
    cvc5::Term body = constraint.substitute(vars, bound);
    cvc5::Term varList = tm.mkTerm(cvc5::Kind::VARIABLE_LIST, bound);
    cvc5::Term existsTerm = tm.mkTerm(cvc5::Kind::EXISTS, {varList, body});
    cvc5::Solver& qe = qeSolver();
    cvc5::Term result = qe.getQuantifierElimination(existsTerm);
    try { result = qe.simplify(result); }
    catch (const cvc5::CVC5ApiException&) {}
    return result;
  }
  catch (const cvc5::CVC5ApiException&)
  {
    ++d_qeFails;
    return constraint;
  }
}

int TheoryDD::rebuildProject(int node, const std::vector<cvc5::Term>& vars,
                             std::unordered_map<int, int>& memo)
{
  if (node <= kTrue)
  {
    return node;  // FALSE / TRUE terminals are unaffected
  }
  auto it = memo.find(node);
  if (it != memo.end())
  {
    return it->second;
  }
  // Copy: makeLeaf/makeNode below may reallocate d_nodes.
  const Node n = d_nodes[node];
  int res;
  if (n.leaf)
  {
    res = makeLeaf(projectConstraint(n.constraint, vars));
  }
  else
  {
    int lo = rebuildProject(n.lo, vars, memo);
    int hi = rebuildProject(n.hi, vars, memo);
    res = makeNode(n.var, lo, hi);
  }
  memo.emplace(node, res);
  return res;
}

int TheoryDD::projectVars(int node, const std::vector<cvc5::Term>& vars)
{
  std::unordered_map<int, int> memo;
  return rebuildProject(node, vars, memo);
}

// Bottom-up apply over the planned schedule. The decision order is fixed by
// plan(); each assertion is compiled from the pre-built atom/variable nodes and
// conjoined into the accumulator one at a time. Stage 2: after the last
// assertion mentioning a real variable is applied, that variable is projected
// out of every leaf so equivalent residual regions merge (bounding the leaf
// count by the frontier rather than the path count).
int TheoryDD::compilePlanned(const std::vector<cvc5::Term>& assertions)
{
  // Bottom nodes first: compile each assertion once into its own diagram (this
  // materializes a leaf per theory atom and an elementary node per decision
  // variable -- the shared building blocks -- and issues each atom's feasibility
  // check up front). Cache the per-assertion node so the apply phase reuses it.
  std::vector<int> assertionNode(assertions.size(), kTrue);
  for (std::size_t i = 0; i < assertions.size(); ++i)
  {
    assertionNode[i] = compileTerm(assertions[i]);
  }

  int acc = kTrue;
  for (std::size_t s = 0; s < d_applySchedule.size(); ++s)
  {
    if (d_progress)
    {
      d_progress(Progress{static_cast<int>(s), d_nodes.size(),
                          d_smtCalls, d_leafCount});
    }
    acc = apply(cvc5::Kind::AND, acc, assertionNode[d_applySchedule[s]]);
    if (acc == kFalse)
    {
      break;  // whole formula is theory-unsatisfiable
    }
    if (d_project && s < d_projectAfter.size() && !d_projectAfter[s].empty())
    {
      acc = projectVars(acc, d_projectAfter[s]);
    }
  }
  return acc;
}

void TheoryDD::compile(const std::vector<cvc5::Term>& assertions)
{
  if (d_mode == Mode::Shannon)
  {
    // The assertions are already on the solver; isolate the exploration in a
    // push so the solver is left untouched afterwards.
    d_solver.push();
    ++d_smtCalls;
    d_root = d_solver.checkSat().isSat() ? buildShannon(0) : kFalse;
    d_solver.pop();
    return;
  }

  if (d_mode == Mode::Planned)
  {
    if (!d_planned)
    {
      plan(assertions);
    }
    d_root = compilePlanned(assertions);
    return;
  }

  // Regions mode: bottom-up apply(AND, .) with residual LRA leaves.
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
  r.qeCalls = d_qeCalls;
  r.qeFails = d_qeFails;
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
  r.feasiblePaths = d_paths;

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
