#pragma once

// Lines 6-22 of Algorithm 1, the streaming union estimator:
//
//    6  p <- 1;  X <- {}
//    7  for i = 1 to m
//    8    t <- c VolEst(P_i, eps', delta')                    c = Thresh 10^{bn}
//    9      for s in X:
//   10        if s in P_i: remove s from X
//   12    while p >= Thresh / t
//   13      remove every element of X with prob. 1/2
//   14      p <- p/2
//   15    N_i <- Poisson(t p)
//   16    while N_i + |X| > Thresh
//   17      remove every element of X with prob. 1/2
//   18      N_i <- Binomial(N_i, 1/2);  p <- p/2
//   19    S <- GenerateSamples(P_i, N_i, b, eps', delta')
//   20    X.Append(S)
//   22  output |X| / (p c)
//
// The loop is generic in its oracles so tests can drive it with exact
// volumes and exact uniform samplers:
//   volEst(i)                           line 8's VolEst factor for P[i]
//   generateSamples(A, b, N, seconds)   line 19 for P[i] = {x : A x <= b}
// and in Coord, the arithmetic of X (double or MpFloat).

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "logger.hpp"
#include "volume/sample_multiset.hpp"
#include "volume/volume_result.hpp"

namespace ttc
{

// p of Algorithm 1, kept as log2(p c): c = Thresh 10^{bn} overflows a double
// and p underflows one long before their product does.
class ScaledProbability
{
 public:
  explicit ScaledProbability(double log2c) : d_log2pc(log2c) {}  // line 6: p = 1

  // log2(t p) for t = c vol.
  double log2tp(double vol) const { return d_log2pc + std::log2(vol); }
  void halve() { d_log2pc -= 1.0; }
  // |X| / (p c), line 22.
  double estimate(std::size_t X) const
  {
    return X ? std::exp2(std::log2(double(X)) - d_log2pc) : 0.0;
  }

 private:
  double d_log2pc;
};

template <typename Coord, typename Polytopes, typename VolEstFn, typename GenerateFn>
void estimateUnionVolume(const Polytopes& P, double thresh, double log2c, unsigned seed,
                         VolEstFn&& volEst, GenerateFn&& generateSamples,
                         VolumeComputationResult& result, const RowCallback& onRow = {})
{
  std::mt19937 rng(seed);       // Poisson, Binomial and the halvings
  ScaledProbability p(log2c);   // line 6
  SampleMultiset<Coord> X;      // line 6
  std::size_t totalGenerated = 0;
  std::size_t totalRemoved = 0;
  double samplingSeconds = 0.0;
  double filterSeconds = 0.0;   // lines 9-18 and 20

  for (std::size_t i = 0; i < P.size(); ++i)  // line 7
  {
    const double vol = volEst(i);  // line 8 (t = c vol)

    const double filterStart = Log.elapsed();
    const std::size_t removed = X.removeInside(P[i].A, P[i].b);  // lines 9-11
    while (p.log2tp(vol) >= std::log2(thresh))                     // line 12: t p >= Thresh
    {
      X.halve(rng);  // line 13
      p.halve();     // line 14
    }
    long long N = std::poisson_distribution<long long>(std::exp2(p.log2tp(vol)))(rng);  // line 15
    while (double(N) + double(X.size()) > thresh)  // line 16
    {
      X.halve(rng);                                             // line 17
      N = std::binomial_distribution<long long>(N, 0.5)(rng);  // line 18
      p.halve();
    }
    filterSeconds += Log.elapsed() - filterStart;

    std::size_t generated = 0;
    if (N > 0)
    {
      double seconds = 0.0;
      const auto S = generateSamples(P[i].A, P[i].b, N, seconds);  // line 19
      samplingSeconds += seconds;
      const double appendStart = Log.elapsed();
      X.append(S);  // line 20
      filterSeconds += Log.elapsed() - appendStart;
      generated = S.size();
    }

    totalGenerated += generated;
    totalRemoved += removed;
    VolumeComputationRow row;
    row.index = i + 1;
    row.volume = vol;
    row.samplesDeleted = removed;
    row.samplesGenerated = generated;
    row.totalSamples = X.size();
    result.rows.push_back(row);
    if (onRow) onRow(row);
  }

  result.finalSampleCount = X.size();
  result.volumeEstimate = p.estimate(X.size());  // line 22
  result.totalSamplesGenerated = totalGenerated;
  result.totalSamplesDeleted = totalRemoved;
  result.samplingTime = samplingSeconds;
  result.filterTime = filterSeconds;
}

// A polytope with a known volume, for driving the loop without VolEst.
struct PolytopeVolume
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  double volume = 0.0;
};

// volEst for estimateUnionVolume that answers from P[i].volume.
template <typename Polytopes>
auto storedVolumes(const Polytopes& P)
{
  return [&P](std::size_t i) { return P[i].volume; };
}

}  // namespace ttc
