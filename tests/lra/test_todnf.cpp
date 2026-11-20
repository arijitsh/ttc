// Line 2 of Algorithm 1, F_bd <- ToDNF(F_bc)  (src/volume/dnf_enumeration.hpp).
//
// The cubes must cover exactly the models of F (each cube implies F, and
// every model of F satisfies some cube), with no repeats and no cube whose
// literals strictly contain another's.

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/boolean_abstraction.hpp"
#include "volume/dnf_enumeration.hpp"

using lra_test::Smt;

namespace
{

bool satisfies(const ttc::Cube& cube, const std::unordered_map<cvc5::Term, bool>& atomValue)
{
  return std::all_of(cube.literals.begin(), cube.literals.end(), [&](const ttc::TermLiteral& l) {
    return atomValue.at(l.term) == l.positive;
  });
}

// Checks that `cubes` cover exactly the models of F over `atoms`.
void requireExactCover(const std::vector<ttc::Cube>& cubes, const cvc5::Term& F,
                       const std::vector<cvc5::Term>& atoms)
{
  for (unsigned mask = 0; mask < (1U << atoms.size()); ++mask)
  {
    std::unordered_map<cvc5::Term, bool> value;
    for (std::size_t j = 0; j < atoms.size(); ++j) value[atoms[j]] = (mask >> j) & 1U;
    const bool covered = std::any_of(cubes.begin(), cubes.end(),
                                     [&](const ttc::Cube& c) { return satisfies(c, value); });
    CHECK(covered == lra_test::evalFormula(F, value));
  }
}

bool literalsContain(const ttc::Cube& larger, const ttc::Cube& smaller)
{
  return std::all_of(smaller.literals.begin(), smaller.literals.end(), [&](const auto& l) {
    return std::any_of(larger.literals.begin(), larger.literals.end(), [&](const auto& k) {
      return k.term == l.term && k.positive == l.positive;
    });
  });
}

void requireIrredundant(const std::vector<ttc::Cube>& cubes)
{
  for (std::size_t i = 0; i < cubes.size(); ++i)
    for (std::size_t j = 0; j < cubes.size(); ++j)
      if (i != j) CHECK_FALSE(literalsContain(cubes[i], cubes[j]));
}

}  // namespace

TEST_CASE("ToDNF: a single atom is one cube with one positive literal")
{
  Smt s;
  const auto a = s.leq(s.real("x"), s.num(1));
  const auto cubes = ttc::enumerateCubes(ttc::buildBooleanAigFromTerm(a));
  REQUIRE(cubes.size() == 1);
  REQUIRE(cubes[0].literals.size() == 1);
  CHECK(cubes[0].literals[0].term == a);
  CHECK(cubes[0].literals[0].positive);
}

TEST_CASE("ToDNF: cubes cover exactly the models of F, without repeats or subsumed cubes")
{
  Smt s;
  const auto x = s.real("x"), y = s.real("y"), z = s.real("z");
  const auto a = s.leq(x, s.num(1)), b = s.leq(y, s.num(2)), c = s.leq(z, s.num(3));
  const std::vector<cvc5::Term> formulas = {
      s.lor({a, b}),
      s.lor({s.land({a, b}), a}),  // a & b is subsumed by a
      s.mk(cvc5::Kind::XOR, {a, b}),
      s.land({s.lor({a, b}), s.lor({s.lnot(a), c})}),
      s.mk(cvc5::Kind::ITE, {a, s.land({b, c}), s.lnot(c)}),
  };
  for (const auto& F : formulas)
  {
    CAPTURE(F);
    const auto cubes = ttc::enumerateCubes(ttc::buildBooleanAigFromTerm(F));
    requireExactCover(cubes, F, {a, b, c});
    requireIrredundant(cubes);
  }
}

TEST_CASE("ToDNF: an unsatisfiable F has no cubes")
{
  Smt s;
  const auto x = s.real("x"), y = s.real("y");
  const auto a = s.leq(x, s.num(1)), b = s.leq(y, s.num(1));
  // Not constant-folded by the builder: (a | b) & !a & !b.
  const auto F = s.land({s.lor({a, b}), s.lnot(a), s.lnot(b)});
  CHECK(ttc::enumerateCubes(ttc::buildBooleanAigFromTerm(F)).empty());
}

TEST_CASE("ToDNF: two inputs that M maps to one atom must agree")
{
  Smt s;
  const auto atom = s.leq(s.real("x"), s.num(1));
  // Inputs v1, v2 (literals 2, 4) both mean `atom`; gate v3 = v1 & v2 or v1 & !v2.
  auto aigFor = [&](std::uint32_t rhs1) {
    ttc::BooleanAbstractionAigData aig;
    aig.maxVariable = 3;
    aig.inputs = {2, 4};
    aig.andGates = {{6, 2, rhs1}};
    aig.outputs = {6};
    aig.mapping = {{1, atom}, {2, atom}};
    return aig;
  };
  const auto agree = ttc::enumerateCubes(aigFor(4));     // v1 & v2
  REQUIRE(agree.size() == 1);
  REQUIRE(agree[0].literals.size() == 1);
  CHECK(agree[0].literals[0].positive);
  CHECK(ttc::enumerateCubes(aigFor(5)).empty());  // v1 & !v2: contradictory
}
