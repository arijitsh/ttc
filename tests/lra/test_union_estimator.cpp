// Lines 6-22 of Algorithm 1, the union loop (src/volume/union_estimator.hpp),
// driven by exact oracles: VolEst returns |P_i ∩ L| / 10^{bn} for the lattice
// L = 10^-b Z^n, and GenerateSamples draws exactly uniform points of P_i ∩ L.
// The loop must then estimate |∪ P_i ∩ L| / 10^{bn}.

#include <cmath>
#include <random>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/mp_float.hpp"
#include "volume/union_estimator.hpp"
#include "volume/union_parameters.hpp"

namespace
{

// An axis-aligned box with its lattice points at precision b.
struct LatticeBox
{
  std::vector<double> lo, hi;

  long first(std::size_t j, int b) const { return std::lround(std::ceil(lo[j] * std::pow(10.0, b))); }
  long last(std::size_t j, int b) const { return std::lround(std::floor(hi[j] * std::pow(10.0, b))); }
  double count(int b) const
  {
    double c = 1.0;
    for (std::size_t j = 0; j < lo.size(); ++j) c *= double(last(j, b) - first(j, b) + 1);
    return c;
  }
};

LatticeBox intersect(const LatticeBox& x, const LatticeBox& y)
{
  LatticeBox z = x;
  for (std::size_t j = 0; j < x.lo.size(); ++j)
  {
    z.lo[j] = std::max(x.lo[j], y.lo[j]);
    z.hi[j] = std::min(x.hi[j], y.hi[j]);
  }
  return z;
}

// |∪ boxes ∩ L| by inclusion-exclusion.
double unionCount(const std::vector<LatticeBox>& boxes, int b)
{
  double total = 0.0;
  for (unsigned mask = 1; mask < (1U << boxes.size()); ++mask)
  {
    LatticeBox meet;
    int size = 0;
    bool empty = false;
    for (std::size_t i = 0; i < boxes.size(); ++i)
      if ((mask >> i) & 1U) meet = size++ == 0 ? boxes[i] : intersect(meet, boxes[i]);
    for (std::size_t j = 0; j < meet.lo.size(); ++j) empty |= meet.first(j, b) > meet.last(j, b);
    if (!empty) total += (size % 2 ? 1.0 : -1.0) * meet.count(b);
  }
  return total;
}

// The loop's inputs for `boxes` at lattice precision b.
struct Instance
{
  std::vector<LatticeBox> boxes;
  int b;
  std::vector<ttc::PolytopeVolume> P;

  Instance(std::vector<LatticeBox> bs, int precision) : boxes(std::move(bs)), b(precision)
  {
    for (const auto& B : boxes)
    {
      const auto H = lra_test::box(B.lo, B.hi);
      P.push_back({H.A, H.b, B.count(b) / std::pow(10.0, b * double(B.lo.size()))});
    }
  }
  std::size_t n() const { return boxes[0].lo.size(); }
  double exact() const { return unionCount(boxes, b) / std::pow(10.0, b * double(n())); }

  // Exactly uniform lattice points of the box whose A is `A`, in Coord.
  template <typename Coord>
  auto sampler(unsigned seed)
  {
    return [this, rng = std::mt19937(seed)](const Eigen::MatrixXd& A, const Eigen::VectorXd&, long long N,
                                           double&) mutable {
      std::size_t i = 0;
      while (&P[i].A != &A) ++i;  // the loop passes P[i].A itself
      const LatticeBox& B = boxes[i];
      std::vector<Eigen::Matrix<Coord, Eigen::Dynamic, 1>> S;
      for (long long s = 0; s < N; ++s)
      {
        Eigen::Matrix<Coord, Eigen::Dynamic, 1> x(n());
        for (std::size_t j = 0; j < n(); ++j)
        {
          const long k = std::uniform_int_distribution<long>(B.first(j, b), B.last(j, b))(rng);
          x(j) = Coord(double(k)) / Coord(std::pow(10.0, b));
        }
        S.push_back(x);
      }
      return S;
    };
  }

  template <typename Coord = double>
  ttc::VolumeComputationResult run(unsigned seed, double eps = 0.8, double delta = 0.2)
  {
    const auto params = ttc::unionParameters(eps, delta, P.size(), n(), b);
    ttc::VolumeComputationResult result;
    ttc::estimateUnionVolume<Coord>(P, params.thresh, params.log2c, seed, ttc::storedVolumes(P),
                                    sampler<Coord>(seed + 1), result);
    return result;
  }
};

}  // namespace

TEST_CASE("union loop: estimates the lattice volume of disjoint, nested and overlapping unions")
{
  const std::vector<Instance> instances = {
      Instance({{{0, 0}, {1, 1}}, {{2, 0}, {3, 2}}}, 1),                       // disjoint
      Instance({{{0, 0}, {4, 4}}, {{1, 1}, {2, 2}}}, 1),                       // nested, big first
      Instance({{{1, 1}, {2, 2}}, {{0, 0}, {4, 4}}}, 1),                       // nested, small first
      Instance({{{0, 0}, {2, 2}}, {{1, 1}, {3, 3}}, {{0.5, 0}, {1.5, 3}}}, 2),  // overlapping
      Instance({{{0, 0, 0}, {1, 1, 1}}, {{0.5, 0.5, 0.5}, {2, 1, 1}}}, 1),     // 3D
  };
  for (auto instance : instances)
  {
    const double exact = instance.exact();
    for (unsigned seed : {1u, 2u, 3u})
    {
      const auto result = instance.run(seed);
      CAPTURE(exact);
      CAPTURE(seed);
      CHECK(result.volumeEstimate == doctest::Approx(exact).epsilon(0.05));
    }
  }
}

TEST_CASE("union loop: |X| stays within Thresh, and the rows account for every sample")
{
  Instance instance({{{0, 0}, {3, 3}}, {{1, 1}, {4, 4}}, {{2, 0}, {5, 2}}}, 3);  // c = Thresh 10^6
  const double thresh = ttc::lraSampleThreshold(0.8, 0.2, 3);
  const auto result = instance.run(9);
  REQUIRE(result.rows.size() == 3);
  std::size_t generated = 0, deleted = 0;
  for (std::size_t i = 0; i < result.rows.size(); ++i)
  {
    const auto& row = result.rows[i];
    CHECK(row.index == i + 1);
    CHECK(double(row.totalSamples) <= thresh);
    CHECK(row.volume == instance.P[i].volume);
    generated += row.samplesGenerated;
    deleted += row.samplesDeleted;
  }
  CHECK(result.totalSamplesGenerated == generated);
  CHECK(result.totalSamplesDeleted == deleted);
  CHECK(result.finalSampleCount == result.rows.back().totalSamples);
  CHECK(result.volumeEstimate == doctest::Approx(instance.exact()).epsilon(0.05));
}

TEST_CASE("union loop: lines 9-11 remove what a later polytope covers, and nothing else")
{
  // Identical boxes: the second deletes every sample the first left.
  Instance same({{{0, 0}, {1, 1}}, {{0, 0}, {1, 1}}}, 1);
  const auto r = same.run(4);
  CHECK(r.rows[1].samplesDeleted == r.rows[0].totalSamples);
  // Disjoint boxes: nothing is ever deleted.
  Instance apart({{{0, 0}, {1, 1}}, {{5, 5}, {6, 6}}}, 1);
  CHECK(apart.run(4).totalSamplesDeleted == 0);
}

TEST_CASE("union loop: deterministic per seed, and independent of the arithmetic of X")
{
  // Integer corners at b = 0: every point and bound is exact in double and
  // GMP, so both make the same membership decisions on the same random stream.
  Instance instance({{{0, 0}, {200, 150}}, {{100, 50}, {300, 300}}, {{250, 0}, {400, 100}}}, 0);
  const auto a = instance.run<double>(17);
  const auto b = instance.run<double>(17);
  CHECK(a.volumeEstimate == b.volumeEstimate);
  CHECK(instance.run<double>(18).volumeEstimate != a.volumeEstimate);

  ttc::MpFloat::default_precision(40);
  const auto g = instance.run<ttc::MpFloat>(17);
  CHECK(g.volumeEstimate == a.volumeEstimate);
  CHECK(g.finalSampleCount == a.finalSampleCount);
}
