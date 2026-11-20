// Lines 4-5 of Algorithm 1 (eps', delta', Thresh) and the scale
// c = Thresh 10^{bn} of lines 8 and 22  (src/volume/union_parameters.hpp),
// plus the run settings fixed before the loop: random streams
// (seed_streams.hpp) and sampling arithmetic (volume_options.hpp).

#include <cmath>
#include <set>
#include <stdexcept>

#include "lra_test_util.hpp"
#include "volume/seed_streams.hpp"
#include "volume/union_parameters.hpp"
#include "volume/volume_options.hpp"

TEST_CASE("parameters: eps' = eps/12 and delta' = delta/(4m)")
{
  const auto p = ttc::unionParameters(0.6, 0.1, 5, 3, 2);
  CHECK(p.epsPrime == doctest::Approx(0.05));
  CHECK(p.deltaPrime == doctest::Approx(0.005));
  // m = 0 (no cubes) counts as one polytope.
  CHECK(ttc::unionParameters(0.6, 0.1, 0, 3, 2).deltaPrime == doctest::Approx(0.025));
}

TEST_CASE("parameters: Thresh is the larger of line 5's two bounds")
{
  for (double eps : {0.05, 0.3, 0.8, 1.0})
    for (double delta : {0.01, 0.2, 0.9})
      for (std::size_t m : {1u, 10u, 100000u})
      {
        const double e = eps / 12.0;
        const double first = 24.0 * std::log(24.0 / delta) / ((1.0 - e) * e * e);
        const double second = 15.0 * (std::log(6.0 / delta) + std::log(double(m)));
        CHECK(ttc::lraSampleThreshold(eps, delta, m) == doctest::Approx(std::max(first, second)));
      }
  // The value ttc prints for its defaults (eps = 0.8, delta = 0.2) on two polytopes.
  CHECK(ttc::lraSampleThreshold(0.8, 0.2, 2) == doctest::Approx(27699.1).epsilon(1e-5));
}

TEST_CASE("parameters: log2 c = log2 Thresh + b n log2 10, and absurd Thresh is rejected")
{
  const auto p = ttc::unionParameters(0.8, 0.2, 2, 7, 30);  // c = Thresh 10^210 overflows nothing
  CHECK(p.log2c == doctest::Approx(std::log2(p.thresh) + 210 * std::log2(10.0)));
  CHECK(ttc::unionParameters(0.8, 0.2, 2, 7, 0).log2c == doctest::Approx(std::log2(p.thresh)));
  CHECK_THROWS_AS(ttc::unionParameters(1e-9, 0.2, 2, 2, 0), std::invalid_argument);
}

TEST_CASE("settings: seed streams are deterministic and distinct per call and per polytope")
{
  ttc::SeedStream s1(42), s2(42), s3(43);
  std::set<unsigned> seen;
  for (int i = 0; i < 1000; ++i)
  {
    const unsigned a = s1.next();
    CHECK(a == s2.next());
    CHECK(a != s3.next());
    seen.insert(a);
  }
  CHECK(seen.size() == 1000);

  std::set<unsigned> oracle;
  for (std::size_t polytope = 1; polytope <= 50; ++polytope)
    for (unsigned attempt = 0; attempt < 3; ++attempt) oracle.insert(ttc::oracleSeed(42, polytope, attempt));
  CHECK(oracle.size() == 150);
  CHECK(ttc::oracleSeed(42, 1, 0) != ttc::oracleSeed(43, 1, 0));
}

TEST_CASE("settings: Auto arithmetic is double up to 15 lattice digits, then GMP")
{
  CHECK(ttc::latticeDigits(6, 3.0) == doctest::Approx(6 + std::log10(3.0)));
  CHECK(ttc::latticeDigits(6, 0.2) == doctest::Approx(6.0));  // |x| < 1 adds no digits
  CHECK(ttc::resolveGmpMode(ttc::GmpMode::Auto, 15.0) == ttc::GmpMode::None);
  CHECK(ttc::resolveGmpMode(ttc::GmpMode::Auto, 15.01) == ttc::GmpMode::PointRepr);
  for (auto mode : {ttc::GmpMode::None, ttc::GmpMode::PointRepr, ttc::GmpMode::Full})
    CHECK(ttc::resolveGmpMode(mode, 40.0) == mode);
  CHECK(ttc::gmpDigits(0, 6) == 38);
  CHECK(ttc::gmpDigits(100, 6) == 100);
}
