// Skinny / badly scaled bodies: failure envelope of ttc's single-polytope
// uniform sampler (src/volume/polytope_sampler.cpp) and a rounding-based fix.
//
// ttc's sampleRestartWalk restarts the walk at the Chebyshev centre for every
// sample and runs 10+n^2 steps on the *unrounded* body.  Every volesti walk
// steps in units of the inner radius r (billiard trajectory length
// 2 sqrt(n) r, compute_diameter.hpp:40; ball walk delta 4r/sqrt(n); hit-and-run
// chords are fine along the long axis but CDHR/RDHR need ~(R/r)^2 steps to
// forget the start), so on a body with R/r >> sqrt(n) the endpoints stay near
// the centre.
//
// The fix (roundedProposals below): round P with volesti's
// inscribed_ellipsoid_rounding (max inscribed ellipsoid -> unit ball), sample
// the rounded body with ttc's own sampler, map back with x = shift + T y.
// Affine maps preserve uniformity, so the samples stay uniform on P.
//
// Uniformity tests (all samples are i.i.d.: every one is a fresh restart):
//  * boxes: z = Q^-1 (x - t) are box coordinates; each coordinate is uniform
//    on its side.  One-sample KS on every coordinate, Bonferroni-corrected:
//    reject if max_i D_i > c(alpha / d) / sqrt(N),  c(a) = sqrt(-ln(a/2)/2).
//  * simplices {z >= 0, sum z <= 1} in dim d: z_1 has CDF 1-(1-t)^d and
//    sum z has CDF t^d; KS on both (Bonferroni over 2).
//  * lattice samples: chi-square over 10 equal slabs of the long axis.
// Family-wise alpha = 0.001 for every CHECK.  "bug" checks FAIL on the
// current sampler; "fixed" checks must pass.
// Slow sweeps are skipped by default:  bin/test_sampling_skinny -ts=slow --no-skip
#include "doctest.h"
#include "common.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"

#include <chrono>
#include <cstdio>
#include <limits>

using namespace vt;
using ttc::SamplerWalk;

namespace
{

constexpr double kAlpha = 0.001;

// ---------------------------------------------------------------- the fix
// Rounding of {A x <= b}: x = shift + T y maps the rounded body {Ar y <= br}
// onto P.  Same call and same failure checks as roundForVolume
// (src/volume/polytope_volume.cpp:272-304): the rounding centre is the
// certified Chebyshev centre, 5 iterations, eigenvalue ratio 6.
struct Rounded
{
  MT Ar, T;
  VT br, shift;
  bool applied = false;
  double seconds = 0;
};

Rounded roundBody(const MT& A, const VT& b)
{
  const auto t0 = std::chrono::steady_clock::now();
  Rounded R;
  const ttc::InnerBall ball = ttc::certifiedInnerBall(A, b);
  if (ball.radius > 0.0)
  {
    HPoly P(static_cast<unsigned>(A.cols()), A, b);
    P.normalize();
    try
    {
      auto [T, shift, det] =
          inscribed_ellipsoid_rounding<MT, VT, double>(P, Point(ball.center), 5, 6.0);
      if (std::isfinite(det) && det > 0.0 && P.get_mat().allFinite()
          && P.get_vec().allFinite() && (P.get_vec().array() > 0.0).all())
      {
        R.Ar = P.get_mat();
        R.br = P.get_vec();
        R.T = T;
        R.shift = shift;
        R.applied = true;
      }
    }
    catch (const std::exception&)
    {
    }
  }
  R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return R;
}

// ttc's sampler on the rounded body, mapped back (falls back to the unrounded
// sampler if rounding fails, exactly as the volume path does).
std::vector<VT> roundedProposals(const MT& A, const VT& b, long long n, unsigned seed,
                                 SamplerWalk walk = SamplerWalk::Billiard, int walkLength = -1,
                                 double* roundSeconds = nullptr)
{
  const Rounded R = roundBody(A, b);
  if (roundSeconds) *roundSeconds += R.seconds;
  if (!R.applied) return ttcProposals(A, b, n, seed, walk, walkLength);
  auto ys = ttcProposals(R.Ar, R.br, n, seed, walk, walkLength);
  for (auto& y : ys) y = R.shift + R.T * y;
  return ys;
}

// ttcLatticeSamples with the rounded generator (lattice rounding in x-space).
std::vector<VT> roundedLatticeSamples(const MT& A, const VT& b, long long n, int precision,
                                      unsigned seed, SamplerWalk walk = SamplerWalk::Billiard,
                                      double epsilon = 0.8, double delta = 0.2)
{
  unsigned walkSeed = seed * 2654435761u + 1u;
  auto generate = [&](const MT& Ai, const VT& bi, long long k, double& t) {
    t = 0;
    return roundedProposals(Ai, bi, k, walkSeed++, walk);
  };
  double samplingTime = 0;
  return ttc::sampleLattice<double>(A, b, n, precision, epsilon / 12.0, delta, seed, generate,
                                    samplingTime);
}

// ---------------------------------------------------------------- bodies
// A test body: P = {A x <= b} is the image x = Q z + t of a reference body
// (a box [lo,hi] or the standard simplex) in z-coordinates.
struct Case
{
  std::string name;
  Body B;
  MT Qinv;
  VT t;
  bool simplex = false;
  std::vector<double> lo, hi;
};

Case boxCase(const std::string& name, const std::vector<double>& lo,
             const std::vector<double>& hi)
{
  Case C;
  C.name = name;
  C.B = box(lo, hi);
  const int d = int(lo.size());
  C.Qinv = MT::Identity(d, d);
  C.t = VT::Zero(d);
  C.lo = lo;
  C.hi = hi;
  return C;
}

// One long axis (the last), sides 1 otherwise.
Case oneLong(int d, double a)
{
  std::vector<double> lo(d, 0.0), hi(d, 1.0);
  hi[d - 1] = a;
  return boxCase("oneLong d=" + std::to_string(d), lo, hi);
}
// ceil(d/2) long axes.
Case severalLong(int d, double a)
{
  std::vector<double> lo(d, 0.0), hi(d, 1.0);
  for (int i = d / 2; i < d; ++i) hi[i] = a;
  return boxCase("sevLong d=" + std::to_string(d), lo, hi);
}
// One thin axis of width 1/a.
Case oneThin(int d, double a)
{
  std::vector<double> lo(d, 0.0), hi(d, 1.0);
  hi[0] = 1.0 / a;
  return boxCase("oneThin d=" + std::to_string(d), lo, hi);
}
// oneLong rotated by a random rotation and translated.
Case rotatedLong(int d, double a, unsigned seed = 11)
{
  Case C = oneLong(d, a);
  const MT Q = randomRotation(d, seed);
  VT t = VT::LinSpaced(d, 1.0, double(d));
  C.B = linearImage(C.B, Q, t, "rot");
  C.Qinv = Q.transpose();
  C.t = t;
  C.name = "rotLong d=" + std::to_string(d);
  return C;
}
// Standard simplex with axis 0 stretched by a.
Case stretchedSimplex(int d, double a)
{
  Case C;
  C.simplex = true;
  MT Q = MT::Identity(d, d);
  Q(0, 0) = a;
  C.B = linearImage(simplex(d), Q, VT::Zero(d), "stretch");
  C.Qinv = Q.inverse();
  C.t = VT::Zero(d);
  C.name = "simplex d=" + std::to_string(d);
  return C;
}
// Cube [c - s, c + s]^d.
Case scaledCube(int d, double s, double c = 0.0)
{
  return boxCase((c == 0.0 ? "cube d=" : "farCube d=") + std::to_string(d),
                 std::vector<double>(d, c - s), std::vector<double>(d, c + s));
}

// ---------------------------------------------------------------- statistics
double ksCrit(double alpha, long long n) { return std::sqrt(-0.5 * std::log(alpha / 2.0)) / std::sqrt(double(n)); }

// One-sample KS statistic of u against the continuous CDF F.
template <typename F>
double ksStat(std::vector<double> u, F cdf)
{
  std::sort(u.begin(), u.end());
  const double n = double(u.size());
  double D = 0;
  for (std::size_t i = 0; i < u.size(); ++i)
  {
    const double f = cdf(u[i]);
    D = std::max({D, (i + 1) / n - f, f - i / n});
  }
  return D;
}

struct Verdict
{
  double ratio = 0;     // max_i D_i / critical value (reject iff > 1)
  double endFrac = 0;   // fraction in the first 10% of the reference axis (uniform: 0.1)
  int outside = 0;      // samples violating A x <= b (+16 ulp of max|b|)
  bool ok() const { return ratio <= 1.0 && outside == 0; }
};

// KS uniformity test of samples x of case C at family-wise level alpha.
Verdict judge(const Case& C, const std::vector<VT>& xs, double alpha = kAlpha)
{
  Verdict V;
  const int d = int(C.B.A.cols());
  const long long N = (long long)xs.size();
  for (const auto& x : xs)
  {
    const VT s = C.B.A * x - C.B.b;
    // A few ulps of the largest |b| (rows of A have unit norm here).
    const double tol = 16 * std::numeric_limits<double>::epsilon()
                       * (1.0 + C.B.b.cwiseAbs().maxCoeff());
    V.outside += (s.maxCoeff() > tol);
  }
  std::vector<VT> zs;
  for (const auto& x : xs) zs.push_back(C.Qinv * (x - C.t));
  if (C.simplex)
  {
    const double crit = ksCrit(alpha / 2.0, N);
    std::vector<double> z0, sum;
    for (const auto& z : zs) z0.push_back(z(0)), sum.push_back(z.sum());
    const double D1 = ksStat(z0, [&](double t) {
      t = std::clamp(t, 0.0, 1.0);
      return 1.0 - std::pow(1.0 - t, d);
    });
    const double D2 = ksStat(sum, [&](double t) { return std::pow(std::clamp(t, 0.0, 1.0), d); });
    V.ratio = std::max(D1, D2) / crit;
    // 10% quantile of z_0: 1-(1-q)^d = 0.1 -> q = 1 - 0.9^(1/d)
    const double q = 1.0 - std::pow(0.9, 1.0 / d);
    V.endFrac = double(std::count_if(z0.begin(), z0.end(), [&](double v) { return v <= q; })) / N;
    return V;
  }
  const double crit = ksCrit(alpha / d, N);
  for (int i = 0; i < d; ++i)
  {
    std::vector<double> u;
    for (const auto& z : zs) u.push_back((z(i) - C.lo[i]) / (C.hi[i] - C.lo[i]));
    V.ratio = std::max(V.ratio, ksStat(u, [](double t) { return std::clamp(t, 0.0, 1.0); }) / crit);
    if (i == d - 1)
      V.endFrac =
          double(std::count_if(u.begin(), u.end(), [](double v) { return v <= 0.1; })) / N;
  }
  return V;
}

const char* walkName(SamplerWalk w)
{
  switch (w)
  {
    case SamplerWalk::Billiard: return "billiard";
    case SamplerWalk::AcceleratedBilliard: return "accbill";
    case SamplerWalk::Ball: return "ball";
    case SamplerWalk::RDHR: return "rdhr";
    case SamplerWalk::CDHR: return "cdhr";
  }
  return "?";
}

void print(const Case& C, double a, SamplerWalk w, const char* tag, const Verdict& V,
           double secs)
{
  std::printf("%-16s a=%-8.0e %-8s %-7s KS/crit=%7.2f end10%%=%.3f out=%d %s  %.3fs\n",
              C.name.c_str(), a, walkName(w), tag, V.ratio, V.endFrac, V.outside,
              V.ok() ? "ok  " : "FAIL", secs);
}

double secondsSince(std::chrono::steady_clock::time_point t0)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// The representative set used by the default (fast) tests.
std::vector<std::pair<Case, double>> representative()
{
  return {{oneLong(2, 500), 500},          {oneLong(5, 100), 100},
          {severalLong(5, 100), 100},      {oneThin(5, 100), 100},
          {rotatedLong(5, 100), 100},      {stretchedSimplex(3, 100), 100},
          {scaledCube(3, 1e-6), 1e-6},     {scaledCube(3, 1e6), 1e6},
          {scaledCube(3, 1e-6, 1e6), 1e-6}};
}

constexpr long long kN = 2000;

}  // namespace

// ------------------------------------------------------------------ default
TEST_CASE("skinny bug: ttc billiard sampler is not uniform on skinny bodies")
{
  // Expected to FAIL on the current sampler for every skinny body (the three
  // cubes are controls: uniform scaling is harmless).  Alpha = 0.001 per body.
  for (auto& [C, a] : representative())
  {
    const auto t0 = std::chrono::steady_clock::now();
    auto xs = ttcProposals(C.B.A, C.B.b, kN, 7);
    const Verdict V = judge(C, xs);
    print(C, a, SamplerWalk::Billiard, "ttc", V, secondsSince(t0));
    CHECK_MESSAGE(V.ok(), C.name << " a=" << a << " KS/crit=" << V.ratio
                                 << " end10%=" << V.endFrac);
  }
}

TEST_CASE("skinny fixed: rounded billiard sampler is uniform on the same bodies")
{
  for (auto& [C, a] : representative())
  {
    const auto t0 = std::chrono::steady_clock::now();
    double rs = 0;
    auto xs = roundedProposals(C.B.A, C.B.b, kN, 7, SamplerWalk::Billiard, -1, &rs);
    const Verdict V = judge(C, xs);
    print(C, a, SamplerWalk::Billiard, "rounded", V, secondsSince(t0));
    CHECK_MESSAGE(V.ok(), C.name << " a=" << a << " KS/crit=" << V.ratio
                                 << " end10%=" << V.endFrac);
  }
}

TEST_CASE("skinny lattice: GenerateSamples on [0,2]x[0,1000] covers all ten slabs")
{
  // Chi-square over 10 equal slabs of y (9 df): reject at alpha = 0.001 when
  // X^2 > 27.877.  The lattice 10^-2 Z is fine enough that each slab holds
  // the same number of lattice points (10000 * 201), so expected counts are equal.
  const Body B = box({0, 0}, {2, 1000});
  const long long n = 2000;
  auto chi2 = [&](const std::vector<VT>& xs) {
    std::vector<int> cnt(10, 0);
    for (const auto& x : xs) cnt[std::min(9, int(x(1) / 100.0))]++;
    double X = 0, e = double(xs.size()) / 10.0;
    for (int c : cnt) X += (c - e) * (c - e) / e;
    return X;
  };
  auto t0 = std::chrono::steady_clock::now();
  const double before = chi2(ttcLatticeSamples(B.A, B.b, n, 2, 5));
  const double tb = secondsSince(t0);
  t0 = std::chrono::steady_clock::now();
  const double after = chi2(roundedLatticeSamples(B.A, B.b, n, 2, 5));
  const double ta = secondsSince(t0);
  std::printf("lattice [0,2]x[0,1000] b=2: chi2 ttc=%.1f (%.3fs) rounded=%.1f (%.3fs), crit 27.88\n",
              before, tb, after, ta);
  CHECK_MESSAGE(before <= 27.877, "ttc lattice sampler chi2=" << before);  // bug: FAILS
  CHECK_MESSAGE(after <= 27.877, "rounded lattice sampler chi2=" << after);
}

TEST_CASE("skinny far: tiny cube [1e6-1e-6, 1e6+1e-6]^d, every walk's samples lie in P")
{
  // Not a mixing problem: in x-space the slacks b - A x ~ 1e-6 are computed
  // next to |b| ~ 1e6, so the walks' boundary intersections lose ~4 of 16
  // digits, and their incrementally updated A x drifts; points land outside
  // P (CDHR even returned an infinite coordinate).  Membership is exact up to
  // 16 ulp of |b| = 3.6e-9 (0.4% of the half-width), and any sample outside
  // is a bug, so the assertion is "zero outside".  Rounding moves the
  // Chebyshev centre to the origin and the body to unit scale, so the walks
  // run in well-conditioned coordinates and only x = shift + T y is rounded.
  for (int d : {2, 5, 10})
    for (SamplerWalk w : {SamplerWalk::Billiard, SamplerWalk::AcceleratedBilliard,
                          SamplerWalk::Ball, SamplerWalk::RDHR, SamplerWalk::CDHR})
    {
      const Case C = scaledCube(d, 1e-6, 1e6);
      const Verdict V0 = judge(C, ttcProposals(C.B.A, C.B.b, 1000, 3, w));
      const Verdict V1 = judge(C, roundedProposals(C.B.A, C.B.b, 1000, 3, w));
      print(C, 1e-6, w, "ttc", V0, 0);
      print(C, 1e-6, w, "rounded", V1, 0);
      CHECK_MESSAGE(V0.outside == 0, C.name << " " << std::string(walkName(w)) << " ttc: " << V0.outside
                                            << "/1000 outside P");  // bug: FAILS
      CHECK_MESSAGE(V1.outside == 0, C.name << " " << std::string(walkName(w)) << " rounded: " << V1.outside
                                            << "/1000 outside P");
    }
}

// ------------------------------------------------------------------ slow
TEST_SUITE("slow")
{
  TEST_CASE("skinny envelope: every walk x body kind x dim x aspect, before/after rounding"
            * doctest::skip())
  {
    // Prints the failure envelope.  Asserted: ttc's default walk (billiard)
    // after rounding passes on every body, Bonferroni over the sweep.
    const std::vector<SamplerWalk> walks = {SamplerWalk::Billiard,
                                            SamplerWalk::AcceleratedBilliard, SamplerWalk::Ball,
                                            SamplerWalk::RDHR, SamplerWalk::CDHR};
    std::vector<std::pair<Case, double>> bodies;
    for (int d : {2, 3, 5, 10})
      for (double a : {1.0, 10.0, 100.0, 1e3, 1e4})
      {
        bodies.push_back({oneLong(d, a), a});
        if (d >= 3) bodies.push_back({severalLong(d, a), a});
        bodies.push_back({oneThin(d, a), a});
        bodies.push_back({rotatedLong(d, a), a});
        bodies.push_back({stretchedSimplex(d, a), a});
      }
    for (int d : {2, 5, 10})
      for (double s : {1e-6, 1e-3, 1e3, 1e6})
      {
        bodies.push_back({scaledCube(d, s), s});
        bodies.push_back({scaledCube(d, s, 1e6), s});
      }
    const double alpha = kAlpha / double(bodies.size());
    const long long N = 1000;
    for (SamplerWalk w : walks)
    {
      int failBefore = 0, failAfter = 0;
      double tBefore = 0, tAfter = 0, tRound = 0;
      for (auto& [C, a] : bodies)
      {
        auto t0 = std::chrono::steady_clock::now();
        const Verdict V0 = judge(C, ttcProposals(C.B.A, C.B.b, N, 3, w), alpha);
        const double s0 = secondsSince(t0);
        t0 = std::chrono::steady_clock::now();
        const Verdict V1 = judge(C, roundedProposals(C.B.A, C.B.b, N, 3, w, -1, &tRound), alpha);
        const double s1 = secondsSince(t0);
        tBefore += s0, tAfter += s1;
        failBefore += !V0.ok(), failAfter += !V1.ok();
        print(C, a, w, "ttc", V0, s0);
        print(C, a, w, "rounded", V1, s1);
        if (w == SamplerWalk::Billiard)
          CHECK_MESSAGE(V1.ok(), C.name << " a=" << a << " rounded KS/crit=" << V1.ratio);
      }
      std::printf("== %s: failing bodies ttc %d/%zu, rounded %d/%zu; time ttc %.2fs, rounded "
                  "%.2fs (of which rounding %.2fs)\n",
                  walkName(w), failBefore, bodies.size(), failAfter, bodies.size(), tBefore,
                  tAfter, tRound);
    }
  }

  TEST_CASE("skinny walk length: is 10+n^2 enough after rounding?" * doctest::skip())
  {
    // Rounded billiard / CDHR on rotated long boxes and stretched simplices
    // in d = 10 for walk lengths 1..10+n^2; prints KS/crit (no assertion
    // below the default length; the default length is asserted).
    for (SamplerWalk w : {SamplerWalk::Billiard, SamplerWalk::Ball, SamplerWalk::CDHR})
      for (auto& [C, a] : std::vector<std::pair<Case, double>>{
               {rotatedLong(10, 1e3), 1e3}, {stretchedSimplex(10, 1e3), 1e3},
               {severalLong(10, 1e3), 1e3}})
        for (int wl : {1, 2, 5, 10, 30, 110, 500})
        {
          auto t0 = std::chrono::steady_clock::now();
          const Verdict V = judge(C, roundedProposals(C.B.A, C.B.b, kN, 9, w, wl), kAlpha / 63);
          char tag[32];
          std::snprintf(tag, sizeof tag, "wl=%d", wl);
          print(C, a, w, tag, V, secondsSince(t0));
          if (wl == 110 && w == SamplerWalk::Billiard)
            CHECK_MESSAGE(V.ok(), C.name << " wl=110 KS/crit=" << V.ratio);
        }
  }
}
