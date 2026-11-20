#pragma once

// Lines 4-5 of Algorithm 1, and the scale c of lines 8 and 22:
//
//   4  eps' <- eps/12,  delta' <- delta/(4m)
//   5  Thresh <- max( 24 ln(24/delta) / ((1 - eps') eps'^2),  15 (ln(6/delta) + ln m) )
//   8  t <- c VolEst(P(c_i), eps', delta')          c = Thresh 10^{bn}
//  22  output |X| / (p c)

#include <cstddef>

namespace ttc
{

struct UnionParameters
{
  double epsPrime = 0.0;    // eps', the per-polytope VolEst and GenerateSamples accuracy
  double deltaPrime = 0.0;  // delta', their confidence
  double thresh = 0.0;      // Thresh, the bound on |X|
  double log2c = 0.0;       // log2(Thresh 10^{bn}); c itself overflows a double for large b n
};

// m polytopes in R^n with lattice precision b.  Throws std::invalid_argument
// when Thresh is not a usable sample bound.
UnionParameters unionParameters(double epsilon, double delta, std::size_t m, std::size_t n, int b);

// Thresh of line 5 alone (m = 0 counts as 1).
double lraSampleThreshold(double epsilon, double delta, std::size_t m);

}  // namespace ttc
