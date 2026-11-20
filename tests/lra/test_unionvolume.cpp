// Algorithm 1 end to end (src/volume/unionvolume.hpp): BoolAbs and ToDNF on a
// cvc5 formula, then computeLraVolume, as main.cpp runs them.

#include <stdexcept>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/boolean_abstraction.hpp"
#include "volume/dnf_enumeration.hpp"
#include "volume/unionvolume.hpp"

using lra_test::Smt;

namespace
{

std::vector<ttc::Cube> cubesOf(const cvc5::Term& F)
{
  return ttc::enumerateCubes(ttc::buildBooleanAigFromTerm(F));
}

}  // namespace

TEST_CASE("Algorithm 1: the area of two overlapping squares is 7")
{
  Smt s;
  const auto x = s.real("x"), y = s.real("y");
  const auto F = s.lor({s.box({x, y}, {0, 0}, {2, 2}), s.box({x, y}, {1, 1}, {3, 3})});
  const auto cubes = cubesOf(F);

  int setups = 0;
  ttc::VolumeSetup setup;
  std::vector<ttc::VolumeComputationRow> rows;
  const auto result = ttc::computeLraVolume(
      cubes, {x, y}, {}, [&](const ttc::VolumeComputationRow& row) { rows.push_back(row); },
      [&](const ttc::VolumeSetup& s) {
        ++setups;
        setup = s;
      });

  CHECK(result.volumeEstimate == doctest::Approx(7.0).epsilon(0.8 / 4));
  CHECK(setups == 1);
  CHECK(setup.dimension == 2);
  CHECK(setup.keptPolytopes == rows.size());
  CHECK(rows.size() == result.rows.size());
  CHECK(setup.latticePrecision == result.latticePrecision);
  CHECK(setup.gmpMode == ttc::GmpMode::None);  // Auto: few lattice digits
  CHECK(result.samplingPrecision == 0);
}

TEST_CASE("Algorithm 1: empty and flat cubes contribute nothing")
{
  Smt s;
  const auto x = s.real("x"), y = s.real("y");
  const auto square = s.box({x, y}, {0, 0}, {1, 1});
  const auto segment = s.land({s.eq(x, s.num(5)), s.geq(y, s.num(0)), s.leq(y, s.num(1))});
  const auto F = s.lor({square, segment});
  ttc::VolumeSetup setup;
  const auto result =
      ttc::computeLraVolume(cubesOf(F), {x, y}, {}, {}, [&](const ttc::VolumeSetup& s) { setup = s; });
  CHECK(setup.keptPolytopes == 1);
  CHECK(result.volumeEstimate == doctest::Approx(1.0).epsilon(0.2));

  // Only flat cubes: volume 0, and the loop never runs.
  const auto flat = ttc::computeLraVolume(cubesOf(segment), {x, y});
  CHECK(flat.volumeEstimate == 0.0);
  CHECK(flat.rows.empty());
}

TEST_CASE("Algorithm 1: deterministic per seed; double and GMP agree on a small instance")
{
  Smt s;
  const auto x = s.real("x"), y = s.real("y"), z = s.real("z");
  const auto F = s.lor({s.box({x, y, z}, {0, 0, 0}, {2, 1, 1}), s.box({x, y, z}, {1, 0, 0}, {3, 2, 1})});
  const auto cubes = cubesOf(F);
  ttc::VolumeOptions options;
  options.gmpMode = ttc::GmpMode::None;
  const double a = ttc::computeLraVolume(cubes, {x, y, z}, options).volumeEstimate;
  CHECK(ttc::computeLraVolume(cubes, {x, y, z}, options).volumeEstimate == a);
  CHECK(a == doctest::Approx(5.0).epsilon(0.2));

  options.gmpMode = ttc::GmpMode::PointRepr;
  CHECK(ttc::computeLraVolume(cubes, {x, y, z}, options).volumeEstimate == a);

  options.gmpMode = ttc::GmpMode::None;
  options.seed = 99;
  CHECK(ttc::computeLraVolume(cubes, {x, y, z}, options).volumeEstimate != a);
}

TEST_CASE("Algorithm 1: invalid options and inputs are rejected")
{
  Smt s;
  const auto x = s.real("x");
  const auto cubes = cubesOf(s.box({x}, {0}, {1}));
  auto with = [](auto set) {
    ttc::VolumeOptions o;
    set(o);
    return o;
  };
  CHECK_THROWS_AS(ttc::computeLraVolume(cubes, {x}, with([](auto& o) { o.epsilon = 0.0; })),
                  std::invalid_argument);
  CHECK_THROWS_AS(ttc::computeLraVolume(cubes, {x}, with([](auto& o) { o.delta = 1.0; })),
                  std::invalid_argument);
  CHECK_THROWS_AS(ttc::computeLraVolume(cubes, {x}, with([](auto& o) { o.sampleWalkLength = 0; })),
                  std::invalid_argument);
  // No real variable to measure in.
  CHECK_THROWS_AS(ttc::computeLraVolume(cubes, {s.integer("k")}), std::runtime_error);
}
