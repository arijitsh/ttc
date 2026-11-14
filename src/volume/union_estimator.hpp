#pragma once

// The union-volume estimator of the ttc algorithm (KR25 manuscript, Alg. 1
// lines 6-22) and its parameters (Thresh, GetPrecision), used by
// lra_volume.cpp and by the tests in tests/volesti/.  The loop is a template
// on the stored coordinate type (double or GMP) and on the sample generator,
// so tests can drive it with exact volumes and exact uniform samplers.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "logger.hpp"
#include "volume/volume_result.hpp"

namespace ttc
{

// Thresh of Alg. 1 line 5 for epsilon/delta and the number of polytopes m:
// max(24 ln(24/delta) / ((1-e') e'^2), 15 (ln(6/delta) + ln m)), e' = epsilon/12.
double lraSampleThreshold(double epsilon, double delta, std::size_t numPolytopes);

// GetPrecision (Alg. 2): lattice precision b such that every polytope's inner
// radius is >= 10^-b * 64 m n / eta * sqrt(ln(4 r / eta)), eta = epsilon/96;
// facetCount is r (total facets), minRadius the smallest inner radius.
int precisionFromCubes(std::size_t dimension, std::size_t facetCount,
                       std::size_t numCubes, double minRadius, double epsilon);

// Stores the running set of sample points "X" used by the union-volume
// (Karp-Luby style) estimator.  Mirrors the behaviour of the reference Python
// implementation in src/cube_processor_nondis.py.  Templated on the coordinate
// type: `double` for the --nogmp path, `MpFloat` for the GMP paths (the
// membership test A x <= b is then evaluated in GMP, avoiding the double
// cancellation that can misclassify points near a facet).
template <typename Coord>
class SampleStore
{
 public:
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;
  using Mat = Eigen::Matrix<Coord, Eigen::Dynamic, Eigen::Dynamic>;

  explicit SampleStore(double tolerance) : d_tolerance(tolerance) {}

  void addSample(const Vec& sample)
  {
    d_samples.push_back(sample);
  }

  // Removes every stored point that lies inside the polytope {x : A x <= b}.
  // Equivalent to the Python:  X = [s for s in X if not poly.is_in(s)].  The
  // half-space matrices are double (built once); they are cast to the
  // coordinate type so the dot products run at the stored precision.
  std::size_t removeInside(const Eigen::MatrixXd& Ad, const Eigen::VectorXd& bd)
  {
    if (Ad.rows() == 0)
    {
      return 0;
    }
    Mat A = Ad.template cast<Coord>();
    Vec b = bd.template cast<Coord>();
    Coord tol(d_tolerance);
    auto oldSize = d_samples.size();
    auto it = std::remove_if(d_samples.begin(), d_samples.end(),
                             [&](const Vec& sample) {
                               Vec diff = A * sample - b;
                               return (diff.array() <= tol).all();
                             });
    d_samples.erase(it, d_samples.end());
    return oldSize - d_samples.size();
  }

  // Keeps each stored point independently with probability 1/2.  Mirrors the
  // Python:  X = [s for s in X if random.random() >= 0.5].
  template <typename URNG>
  std::size_t downsampleHalf(URNG& rng)
  {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto oldSize = d_samples.size();
    auto it = std::remove_if(d_samples.begin(), d_samples.end(),
                             [&](const Vec&) {
                               return unit(rng) < 0.5;
                             });
    d_samples.erase(it, d_samples.end());
    return oldSize - d_samples.size();
  }

  std::size_t size() const
  {
    return d_samples.size();
  }

 private:
  double d_tolerance;
  std::vector<Vec> d_samples;
};

// A polytope kept after phase 1, with its (double) half-space matrices and the
// volesti volume estimate.
struct PolytopeVolume
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  double volume = 0.0;
};

// Phase 2: streaming union-volume (Karp-Luby style) estimation, templated on the
// coordinate type used to store the surviving sample set `X`.  `p` is the
// running sampling probability; the union volume is |X| / (p*C), C=T*10^(b*n).
// `generate(A, b, n, samplingTime)` draws `n` fresh points inside {x : A x <= b}
// (the only step that differs between the double and GMP walk backends).
template <typename Coord, typename GenFn>
inline void runUnionAlgorithm(const std::vector<PolytopeVolume>& kept, double thresh,
                       unsigned seed, double log2Scale, GenFn&& generate,
                       VolumeComputationResult& result,
                       const std::function<void(const VolumeComputationRow&)>& onRow)
{
  std::mt19937 auxRng(seed);  // drives Poisson draws and down-sampling
  SampleStore<Coord> store(0.0);

  // Track log2(p*C), avoiding overflow of C and underflow of p.
  double log2Density = log2Scale;
  std::size_t rowIndex = 0;
  std::size_t totalGenerated = 0;
  std::size_t totalDeleted = 0;
  double totalSamplingTime = 0.0;
  double totalFilterTime = 0.0;

  for (const auto& entry : kept)
  {
    const double volume = entry.volume;

    // Drop stored points that lie inside the current polytope.
    const double filterStart = Log.elapsed();
    std::size_t deleted = store.removeInside(entry.A, entry.b);

    // Shrink the sampling probability until p*C*volume < thresh.
    while (log2Density + std::log2(volume) >= std::log2(thresh))
    {
      store.downsampleHalf(auxRng);
      log2Density -= 1.0;
    }

    // Draw the number of fresh samples; keep |X| + N within the threshold.
    std::poisson_distribution<long long> poisson(std::exp2(log2Density + std::log2(volume)));
    long long n = poisson(auxRng);
    while (static_cast<double>(n) + static_cast<double>(store.size()) > thresh)
    {
      store.downsampleHalf(auxRng);
      log2Density -= 1.0;
      n = std::binomial_distribution<long long>(n, 0.5)(auxRng);
    }

    totalFilterTime += Log.elapsed() - filterStart;

    std::size_t sampleCount = 0;
    double samplingTime = 0.0;
    if (n > 0)
    {
      auto samples = generate(entry.A, entry.b, n, samplingTime);
      const double storeStart = Log.elapsed();
      for (auto& sample : samples)
      {
        store.addSample(sample);
      }
      totalFilterTime += Log.elapsed() - storeStart;
      sampleCount = samples.size();
    }

    totalGenerated += sampleCount;
    totalDeleted += deleted;
    totalSamplingTime += samplingTime;
    VolumeComputationRow row;
    row.index = ++rowIndex;
    row.volume = volume;
    row.samplesDeleted = deleted;
    row.totalSamples = store.size();
    row.samplesGenerated = sampleCount;
    result.rows.push_back(row);
    if (onRow)
    {
      onRow(row);
    }
  }

  result.finalSampleCount = store.size();
  result.volumeEstimate =
      store.size() ? std::exp2(std::log2(double(store.size())) - log2Density) : 0.0;
  result.totalSamplesGenerated = totalGenerated;
  result.totalSamplesDeleted = totalDeleted;
  result.samplingTime = totalSamplingTime;
  result.filterTime = totalFilterTime;
}

}  // namespace ttc
