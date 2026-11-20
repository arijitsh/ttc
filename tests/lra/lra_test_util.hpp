#pragma once

// Helpers shared by the unit tests of the LRA volume engine (tests/lra/).

#include <cmath>
#include <cstdint>
#include <functional>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <cvc5/cvc5.h>

#include "doctest.h"
#include "features.hpp"

namespace lra_test
{

// A cvc5 context for building formulas.
struct Smt
{
  cvc5::TermManager tm;
  cvc5::Solver solver{tm};

  cvc5::Term real(const std::string& name) { return tm.mkConst(tm.getRealSort(), name); }
  cvc5::Term integer(const std::string& name) { return tm.mkConst(tm.getIntegerSort(), name); }
  cvc5::Term boolean(const std::string& name) { return tm.mkConst(tm.getBooleanSort(), name); }
  cvc5::Term num(const std::string& value) { return tm.mkReal(value); }
  cvc5::Term num(std::int64_t value) { return tm.mkReal(value); }
  cvc5::Term mk(cvc5::Kind kind, const std::vector<cvc5::Term>& children)
  {
    return tm.mkTerm(kind, children);
  }
  cvc5::Term leq(const cvc5::Term& a, const cvc5::Term& b) { return mk(cvc5::Kind::LEQ, {a, b}); }
  cvc5::Term geq(const cvc5::Term& a, const cvc5::Term& b) { return mk(cvc5::Kind::GEQ, {a, b}); }
  cvc5::Term eq(const cvc5::Term& a, const cvc5::Term& b) { return mk(cvc5::Kind::EQUAL, {a, b}); }
  cvc5::Term lnot(const cvc5::Term& a) { return mk(cvc5::Kind::NOT, {a}); }
  cvc5::Term land(const std::vector<cvc5::Term>& c) { return mk(cvc5::Kind::AND, c); }
  cvc5::Term lor(const std::vector<cvc5::Term>& c) { return mk(cvc5::Kind::OR, c); }

  // lo <= x <= hi for every (x, lo, hi).
  cvc5::Term box(const std::vector<cvc5::Term>& vars, const std::vector<std::int64_t>& lo,
                 const std::vector<std::int64_t>& hi)
  {
    std::vector<cvc5::Term> c;
    for (std::size_t j = 0; j < vars.size(); ++j)
    {
      c.push_back(geq(vars[j], num(lo[j])));
      c.push_back(leq(vars[j], num(hi[j])));
    }
    return land(c);
  }
};

// The truth value of a Boolean combination of atoms (And, Or, Not, Implies,
// Ite, Xor, Boolean =) under `atomValue`.
inline bool evalFormula(const cvc5::Term& t, const std::unordered_map<cvc5::Term, bool>& atomValue)
{
  using cvc5::Kind;
  auto ev = [&](const cvc5::Term& c) { return evalFormula(c, atomValue); };
  switch (t.getKind())
  {
    case Kind::CONST_BOOLEAN: return t.getBooleanValue();
    case Kind::NOT: return !ev(t[0]);
    case Kind::AND:
      for (const auto& c : t)
        if (!ev(c)) return false;
      return true;
    case Kind::OR:
      for (const auto& c : t)
        if (ev(c)) return true;
      return false;
    case Kind::IMPLIES: return !ev(t[0]) || ev(t[1]);
    case Kind::ITE: return ev(t[0]) ? ev(t[1]) : ev(t[2]);
    case Kind::XOR: return ev(t[0]) != ev(t[1]);
    case Kind::EQUAL:
      if (t[0].getSort().isBoolean()) return ev(t[0]) == ev(t[1]);
      [[fallthrough]];
    default: return atomValue.at(t);
  }
}

// The output of an AIG for the input values `inputValue` (variable -> value).
inline bool evalAig(const ttc::BooleanAbstractionAigData& aig,
                    const std::unordered_map<std::uint32_t, bool>& inputValue)
{
  std::unordered_map<std::uint32_t, const ttc::BooleanAbstractionAigData::AndGate*> gateOf;
  for (const auto& g : aig.andGates) gateOf[g.lhs >> 1U] = &g;
  std::unordered_map<std::uint32_t, bool> memo;
  std::function<bool(std::uint32_t)> lit = [&](std::uint32_t l) -> bool {
    if (l <= 1U) return l == 1U;
    const std::uint32_t v = l >> 1U;
    bool value;
    if (auto m = memo.find(v); m != memo.end())
      value = m->second;
    else if (auto g = gateOf.find(v); g != gateOf.end())
      value = memo[v] = lit(g->second->rhs0) && lit(g->second->rhs1);
    else
      value = memo[v] = inputValue.at(v);
    return (l & 1U) ? !value : value;
  };
  REQUIRE(aig.outputs.size() == 1);
  return lit(aig.outputs[0]);
}

// The axis-aligned box prod [lo_j, hi_j] as {x : A x <= b}.
struct Box
{
  std::vector<double> lo, hi;
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  double volume() const
  {
    double v = 1.0;
    for (std::size_t j = 0; j < lo.size(); ++j) v *= hi[j] - lo[j];
    return v;
  }
  bool contains(const Eigen::VectorXd& x) const
  {
    for (std::size_t j = 0; j < lo.size(); ++j)
      if (x(j) < lo[j] || x(j) > hi[j]) return false;
    return true;
  }
};

inline Box box(const std::vector<double>& lo, const std::vector<double>& hi)
{
  const long n = static_cast<long>(lo.size());
  Box B{lo, hi, Eigen::MatrixXd::Zero(2 * n, n), Eigen::VectorXd::Zero(2 * n)};
  for (long j = 0; j < n; ++j)
  {
    B.A(2 * j, j) = 1.0;
    B.b(2 * j) = hi[j];
    B.A(2 * j + 1, j) = -1.0;
    B.b(2 * j + 1) = -lo[j];
  }
  return B;
}

// The standard simplex {x >= 0, sum x <= 1} in R^n (volume 1/n!).
inline void simplex(int n, Eigen::MatrixXd& A, Eigen::VectorXd& b)
{
  A = Eigen::MatrixXd::Zero(n + 1, n);
  b = Eigen::VectorXd::Zero(n + 1);
  for (int j = 0; j < n; ++j) A(j, j) = -1.0;
  A.row(n).setOnes();
  b(n) = 1.0;
}

// Exactly uniform points of the axis-aligned box {x : A x <= b} (A, b of the
// form box() builds, possibly with b moved outwards): a stand-in for the
// sampling walk with no mixing error.
struct UniformBoxSampler
{
  std::mt19937 rng{2024};
  std::vector<Eigen::VectorXd> operator()(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                          long long count, double&)
  {
    const long n = A.cols();
    std::vector<Eigen::VectorXd> out;
    for (long long i = 0; i < count; ++i)
    {
      Eigen::VectorXd x(n);
      for (long j = 0; j < n; ++j)
      {
        // Rows 2j and 2j+1 bound coordinate j: x_j <= b(2j)/A(2j,j), -x_j <= ...
        const double hi = b(2 * j) / A(2 * j, j);
        const double lo = b(2 * j + 1) / A(2 * j + 1, j);
        x(j) = std::uniform_real_distribution<double>(lo, hi)(rng);
      }
      out.push_back(x);
    }
    return out;
  }
};

}  // namespace lra_test
