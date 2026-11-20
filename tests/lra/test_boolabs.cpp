// Line 1 of Algorithm 1, F_bc, M <- BoolAbs(F)  (src/volume/boolean_abstraction.hpp).
//
// The AIG must compute F when every input is read as its atom under M.

#include <cstdint>
#include <random>
#include <unordered_map>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/boolean_abstraction.hpp"
#include "volume/linear_constraints.hpp"

using lra_test::Smt;

namespace
{

// Checks that `aig` computes `formula` on every assignment of `atoms`.
void requireSameFunction(const ttc::BooleanAbstractionAigData& aig, const cvc5::Term& formula,
                         const std::vector<cvc5::Term>& atoms)
{
  std::unordered_map<cvc5::Term, std::uint32_t> inputOf;
  for (const auto& [var, atom] : aig.mapping) inputOf[atom] = var;
  for (unsigned mask = 0; mask < (1U << atoms.size()); ++mask)
  {
    std::unordered_map<cvc5::Term, bool> atomValue;
    std::unordered_map<std::uint32_t, bool> inputValue;
    for (std::size_t j = 0; j < atoms.size(); ++j)
    {
      const bool v = (mask >> j) & 1U;
      atomValue[atoms[j]] = v;
      if (inputOf.count(atoms[j])) inputValue[inputOf.at(atoms[j])] = v;
    }
    CHECK(lra_test::evalAig(aig, inputValue) == lra_test::evalFormula(formula, atomValue));
  }
}

}  // namespace

TEST_CASE("BoolAbs: every theory atom becomes one AIG input, and M maps it back")
{
  Smt s;
  const auto x = s.real("x"), y = s.real("y");
  const auto a = s.leq(x, s.num(1)), b = s.geq(y, s.num(0));
  // a appears twice: it must still be a single input.
  const auto F = s.land({a, s.lor({b, s.lnot(a)}), a});
  const auto aig = ttc::buildBooleanAigFromTerm(F);

  REQUIRE(aig.inputs.size() == 2);
  REQUIRE(aig.mapping.size() == 2);
  for (std::size_t i = 0; i < aig.inputs.size(); ++i)
    CHECK(aig.mapping[i].first == aig.inputs[i] >> 1U);  // both sorted by variable
  const bool mapsA = aig.mapping[0].second == a || aig.mapping[1].second == a;
  const bool mapsB = aig.mapping[0].second == b || aig.mapping[1].second == b;
  CHECK(mapsA);
  CHECK(mapsB);
  requireSameFunction(aig, F, {a, b});
}

TEST_CASE("BoolAbs: Implies, Ite, Xor and Boolean = compute the formula")
{
  Smt s;
  const auto x = s.real("x"), y = s.real("y"), z = s.real("z");
  const auto a = s.leq(x, s.num(1)), b = s.leq(y, s.num(2)), c = s.eq(z, s.num(3));
  const std::vector<cvc5::Term> formulas = {
      s.mk(cvc5::Kind::IMPLIES, {a, b}),
      s.mk(cvc5::Kind::ITE, {a, b, c}),
      s.mk(cvc5::Kind::XOR, {a, s.lnot(c)}),
      s.eq(a, s.lor({b, c})),
      s.lor({s.land({a, b}), s.land({s.lnot(a), c}), s.mk(cvc5::Kind::XOR, {b, c})}),
  };
  for (const auto& F : formulas)
  {
    CAPTURE(F);
    requireSameFunction(ttc::buildBooleanAigFromTerm(F), F, {a, b, c});
  }
}

TEST_CASE("BoolAbs: constants fold and equal sub-circuits share one gate")
{
  Smt s;
  const auto x = s.real("x");
  const auto a = s.leq(x, s.num(1)), b = s.geq(x, s.num(0));

  // a | !a is the constant true (AIGER literal 1), a & !a the constant false (0).
  CHECK(ttc::buildBooleanAigFromTerm(s.lor({a, s.lnot(a)})).outputs.at(0) == 1U);
  CHECK(ttc::buildBooleanAigFromTerm(s.land({a, s.lnot(a)})).outputs.at(0) == 0U);
  CHECK(ttc::buildBooleanAigFromTerm(s.land({s.tm.mkTrue(), a})).andGates.empty());

  // (a & b) | (b & a): one AND gate, and the Or of a gate with itself is the gate.
  const auto aig = ttc::buildBooleanAigFromTerm(s.lor({s.land({a, b}), s.land({b, a})}));
  CHECK(aig.andGates.size() == 1);
}

TEST_CASE("BoolAbs: cvc5's abstraction (when built in) defines the same set of points")
{
  // cvc5 may rewrite atoms (y <= 2 as not (y > 2), ...), so compare the two
  // AIGs pointwise: read every atom at a point through its linear constraints.
  Smt s;
  const auto x = s.real("x"), y = s.real("y");
  const ttc::RealVariables vars({x, y});
  const auto a = s.leq(x, s.num(1)), b = s.leq(y, s.num(2));
  const auto c = s.geq(s.mk(cvc5::Kind::ADD, {x, y}), s.num(0));
  const auto F = s.lor({s.land({a, b}), s.land({s.lnot(a), c})});
  ttc::BooleanAbstractionAigData fromCvc5;
  try
  {
    fromCvc5 = ttc::getBooleanAbstractionAig(s.solver, F);
  }
  catch (const std::exception& e)
  {
    MESSAGE("cvc5 AIG abstraction unavailable: " << std::string(e.what()));
    return;
  }
  const auto fallback = ttc::buildBooleanAigFromTerm(F);

  auto holds = [&](const cvc5::Term& atom, const Eigen::Vector2d& p) {
    const auto P = ttc::cubeConstraints(ttc::Cube{{{atom, true}}}, vars);
    return ((P.A * p - P.b).array() <= 0.0).all();
  };
  auto inputs = [&](const ttc::BooleanAbstractionAigData& aig, const Eigen::Vector2d& p) {
    std::unordered_map<std::uint32_t, bool> value;
    for (const auto& [var, atom] : aig.mapping) value[var] = holds(atom, p);
    return value;
  };
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> coord(-5.0, 5.0);
  for (int i = 0; i < 500; ++i)
  {
    const Eigen::Vector2d p(coord(rng), coord(rng));
    CAPTURE(p.transpose());
    CHECK(lra_test::evalAig(fromCvc5, inputs(fromCvc5, p))
          == lra_test::evalAig(fallback, inputs(fallback, p)));
  }
}
