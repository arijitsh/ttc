// Mixing of ttc's single-polytope sampler (src/volume/polytope_sampler.cpp).
//
// sampleRestartWalk restarts the walk at the certified Chebyshev centre c for
// EVERY sample and keeps the endpoint after L steps, so its output law is the
// L-step law of the walk from the point mass at c -- never the stationary
// (uniform) law, only close to it once L exceeds the mixing time from c.
//
// Statistic.  For a convex body P with c in its interior, the gauge
//   g(x) = max_i a_i.(x - c) / (b_i - a_i.c)
// of a uniform point X satisfies P(g(X) <= t) = vol(tP)/vol(P) = t^n exactly.
// A walk that has not left the neighbourhood of c has too much mass at small
// g, so a one-sample Kolmogorov-Smirnov test of {g(X_k)} against F(t) = t^n is
// a direct test of centre bias.  All tests use N = 5000 samples and reject at
// alpha = 0.001 (critical D = 1.949/sqrt(N) = 0.0276).  For boxes the
// coordinate marginals are exactly uniform, so we also KS-test the long-axis
// marginal.  Seeds are fixed, so results are deterministic.
//
// Default suite: ttc's default walk length L = 10 + n^2 for every walk on
// round and elongated bodies; independence across consecutive samples and
// across per-polytope seed streams; and the candidate fix (round P with its
// max inscribed ellipsoid, sample there, map back).  test_suite("slow"): sweep L to find the
// smallest passing walk length per (walk, body), and the one-long-chain
// alternative (burn-in + thinning).
#include "doctest.h"
#include "common.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"

#include <cstdint>
#include <cstdio>
#include <list>

using namespace vt;
using ttc::SamplerWalk;

namespace
{

constexpr int kN = 5000;
constexpr double kAlpha = 0.001;

const char* walkName(SamplerWalk w)
{
  switch (w)
  {
    case SamplerWalk::Billiard: return "billiard";
    case SamplerWalk::AcceleratedBilliard: return "accel-billiard";
    case SamplerWalk::Ball: return "ball";
    case SamplerWalk::RDHR: return "rdhr";
    case SamplerWalk::CDHR: return "cdhr";
  }
  return "?";
}
const SamplerWalk kWalks[] = {SamplerWalk::Billiard, SamplerWalk::AcceleratedBilliard,
                              SamplerWalk::Ball, SamplerWalk::RDHR, SamplerWalk::CDHR};

// Kolmogorov survival function Q(lambda) = 2 sum (-1)^{k-1} exp(-2 k^2 lambda^2).
double kolmogorovQ(double lambda)
{
  if (lambda < 0.2) return 1.0;
  double s = 0;
  for (int k = 1; k <= 100; ++k)
    s += (k % 2 ? 2.0 : -2.0) * std::exp(-2.0 * k * k * lambda * lambda);
  return std::min(1.0, std::max(0.0, s));
}

struct KS
{
  double D = 1, p = 0;
};

// One-sample KS of xs against the continuous CDF F.
template <typename F>
KS ksTest(std::vector<double> xs, F cdf)
{
  std::sort(xs.begin(), xs.end());
  const double N = xs.size();
  double D = 0;
  for (std::size_t i = 0; i < xs.size(); ++i)
  {
    const double f = cdf(xs[i]);
    D = std::max({D, f - i / N, (i + 1) / N - f});
  }
  const double sn = std::sqrt(N);
  return {D, kolmogorovQ((sn + 0.12 + 0.11 / sn) * D)};
}

// Gauge (Minkowski functional) of x w.r.t. c; slack_i = b_i - a_i.c > 0.
double gauge(const MT& A, const VT& b, const VT& c, const VT& x)
{
  const VT slack = b - A * c;
  const VT ax = A * (x - c);
  double g = 0;
  for (Eigen::Index i = 0; i < A.rows(); ++i) g = std::max(g, ax(i) / slack(i));
  return g;
}

std::vector<double> gauges(const Body& B, const VT& c, const std::vector<VT>& pts)
{
  std::vector<double> g;
  g.reserve(pts.size());
  for (const auto& p : pts) g.push_back(gauge(B.A, B.b, c, p));
  return g;
}

VT chebyshevCentre(const Body& B) { return ttc::certifiedInnerBall(B.A, B.b).center; }

KS gaugeKS(const Body& B, const VT& c, const std::vector<VT>& pts)
{
  const double n = B.A.cols();
  return ksTest(gauges(B, c, pts), [n](double t) { return std::pow(std::clamp(t, 0.0, 1.0), n); });
}

// KS of coordinate j against uniform on [lo, hi] (exact for axis-aligned boxes).
KS marginalKS(const std::vector<VT>& pts, int j, double lo, double hi)
{
  std::vector<double> v;
  for (const auto& p : pts) v.push_back(p(j));
  return ksTest(v, [lo, hi](double t) { return std::clamp((t - lo) / (hi - lo), 0.0, 1.0); });
}

// Pearson correlation of (x_k, y_{k+lag}).
double lagCorr(const std::vector<double>& x, const std::vector<double>& y, int lag)
{
  const std::size_t m = std::min(x.size(), y.size()) - lag;
  double mx = 0, my = 0;
  for (std::size_t k = 0; k < m; ++k) mx += x[k], my += y[k + lag];
  mx /= m, my /= m;
  double sxy = 0, sxx = 0, syy = 0;
  for (std::size_t k = 0; k < m; ++k)
  {
    const double dx = x[k] - mx, dy = y[k + lag] - my;
    sxy += dx * dy, sxx += dx * dx, syy += dy * dy;
  }
  return sxy / std::sqrt(sxx * syy);
}

// ------------------------------------------------------------------ bodies
// Ball-like: m = 8n random tangent halfspaces of the unit ball, intersected
// with [-2,2]^n so the body is bounded.
Body randomFacets(int n, unsigned seed)
{
  Body C = cube(n, -2, 2);
  const int m = 8 * n;
  std::mt19937 g(seed);
  std::normal_distribution<double> N01(0, 1);
  Body B;
  B.name = "randfacet" + std::to_string(n);
  B.A = MT(C.A.rows() + m, n);
  B.b = VT(C.A.rows() + m);
  B.A.topRows(C.A.rows()) = C.A;
  B.b.head(C.A.rows()) = C.b;
  for (int i = 0; i < m; ++i)
  {
    VT a(n);
    for (int j = 0; j < n; ++j) a(j) = N01(g);
    B.A.row(C.A.rows() + i) = a.normalized().transpose();
    B.b(C.A.rows() + i) = 1.0;
  }
  return B;
}

// [-1,1]^n with the last axis stretched by `aspect` (optionally rotated).
Body stretchedCube(int n, double aspect, bool rotate)
{
  MT D = MT::Identity(n, n);
  D(n - 1, n - 1) = aspect;
  MT Q = rotate ? MT(randomRotation(n, 11) * D) : D;
  char tag[64];
  std::snprintf(tag, sizeof tag, "stretch%g%s", aspect, rotate ? "rot" : "");
  Body B = linearImage(cube(n), Q, VT::Zero(n), tag);
  B.name = "cube" + std::to_string(n) + "-" + tag;
  return B;
}

Body namedBox(const std::vector<double>& lo, const std::vector<double>& hi, const char* name)
{
  Body B = box(lo, hi);
  B.name = name;
  return B;
}

std::vector<Body> roundBodies()
{
  std::vector<Body> v = {cube(2), cube(5), cube(10), crossPolytope(6), randomFacets(8, 5)};
  v[0].name = "cube2";
  v[1].name = "cube5";
  v[2].name = "cube10";
  return v;
}

// --------------------------------------------- long chain (alternative design)
// One chain from the centre: `burnIn` steps, then a sample every `thin` steps.
std::vector<VT> chainSamples(const Body& B, SamplerWalk walk, unsigned burnIn, unsigned thin,
                             long long n, unsigned seed)
{
  const ttc::InnerBall ball = ttc::certifiedInnerBall(B.A, B.b);
  HPoly P(static_cast<unsigned>(B.A.cols()), B.A, B.b);
  P.normalize();
  P.set_InnerBall({Point(ball.center), ball.radius});
  TtcSampleRNG rng(static_cast<int>(B.A.cols()));
  rng.set_seed(seed);
  Point p(ball.center);
  std::list<Point> pts;
  PushBackWalkPolicy policy;
  auto run = [&](auto tag) {
    using W = typename decltype(tag)::template Walk<HPoly, TtcSampleRNG>;
    W w(P, p, rng);
    w.apply(P, p, burnIn, rng);
    for (long long k = 0; k < n; ++k)
    {
      w.apply(P, p, thin, rng);
      policy.apply(pts, p);
    }
  };
  switch (walk)
  {
    case SamplerWalk::Billiard: run(BilliardWalk()); break;
    case SamplerWalk::AcceleratedBilliard: run(AcceleratedBilliardWalk()); break;
    case SamplerWalk::Ball: run(BallWalk()); break;
    case SamplerWalk::RDHR: run(RDHRWalk()); break;
    case SamplerWalk::CDHR: run(CDHRWalk()); break;
  }
  std::vector<VT> out;
  for (const auto& q : pts) out.push_back(q.getCoefficients());
  return out;
}

// Candidate fix (no src change): round P with volesti's max-inscribed-
// ellipsoid rounding exactly as polytope_volume.cpp roundForVolume() does
// (5 iterations, eigenvalue ratio 6, centred at the certified Chebyshev
// centre), run ttc's own sampleRestartWalk on the rounded body
// {y : A' y <= b'}, and map back x = shift + T y.  Affine maps preserve
// uniformity, so the output is uniform on P iff the rounded sampler is.
std::vector<VT> roundedProposals(const Body& B, long long n, unsigned seed, SamplerWalk w,
                                 int walkLength = -1)
{
  const ttc::InnerBall ball = ttc::certifiedInnerBall(B.A, B.b);
  HPoly P(static_cast<unsigned>(B.A.cols()), B.A, B.b);
  P.normalize();
  auto [T, shift, det] =
      inscribed_ellipsoid_rounding<MT, VT, double>(P, Point(ball.center), 5, 6.0);
  (void)det;
  auto ys = ttcProposals(P.get_mat(), P.get_vec(), n, seed, w, walkLength);
  std::vector<VT> xs;
  xs.reserve(ys.size());
  for (const auto& y : ys) xs.push_back(shift + T * y);
  return xs;
}

// splitmix64 per-polytope stream seeds exactly as lra_volume.cpp streamSeed().
unsigned streamSeed(unsigned base, unsigned index)
{
  std::uint64_t z = base + 0x9E3779B97F4A7C15ULL * index;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return static_cast<unsigned>(z ^ (z >> 31));
}

}  // namespace

// ===================================================================== default

TEST_CASE("mixing: gauge KS at ttc's default walk length, round bodies")
{
  for (const Body& B : roundBodies())
  {
    const int n = int(B.A.cols());
    const VT c = chebyshevCentre(B);
    for (SamplerWalk w : kWalks)
    {
      auto pts = ttcProposals(B.A, B.b, kN, 1234, w);
      const KS ks = gaugeKS(B, c, pts);
      std::printf("[default L=%u] %-12s %-14s gauge KS D=%.4f p=%.2e\n", ttcWalkLength(n),
                  B.name.c_str(), walkName(w), ks.D, ks.p);
      CHECK_MESSAGE(ks.p > kAlpha, B.name << " " << std::string(walkName(w)) << " D=" << ks.D);
    }
  }
}

TEST_CASE("mixing: gauge + long-axis KS at ttc's default walk length, elongated bodies")
{
  struct Case
  {
    Body B;
    double lo, hi;  // long-axis range if axis-aligned box, else lo > hi
  };
  std::vector<Case> cases = {{namedBox({0, 0}, {2, 1000}, "box2x1000"), 0, 1000},
                             {namedBox({0, 0}, {1, 10}, "box1x10"), 0, 10},
                             {stretchedCube(5, 10, false), -10, 10},
                             {stretchedCube(5, 10, true), 1, 0}};
  for (auto& [B, lo, hi] : cases)
  {
    const int n = int(B.A.cols());
    const VT c = chebyshevCentre(B);
    for (SamplerWalk w : kWalks)
    {
      auto pts = ttcProposals(B.A, B.b, kN, 1234, w);
      const KS ks = gaugeKS(B, c, pts);
      KS mk;
      if (lo < hi) mk = marginalKS(pts, n - 1, lo, hi);
      std::printf("[default L=%u] %-22s %-14s gauge D=%.4f p=%.2e", ttcWalkLength(n),
                  B.name.c_str(), walkName(w), ks.D, ks.p);
      if (lo < hi) std::printf("  long-axis D=%.4f p=%.2e", mk.D, mk.p);
      std::printf("\n");
      CHECK_MESSAGE(ks.p > kAlpha, B.name << " " << std::string(walkName(w)) << " gauge D=" << ks.D);
      if (lo < hi)
        CHECK_MESSAGE(mk.p > kAlpha, B.name << " " << std::string(walkName(w)) << " long-axis D=" << mk.D);
    }
  }
}

// Candidate fix: sample in rounded coordinates (roundedProposals above) at
// ttc's default walk length.  Billiard (ttc's default walk) must pass the
// gauge KS and the long-axis KS on every elongated body; other walks printed.
TEST_CASE("mixing: rounded-coordinate sampling at default walk length, elongated bodies")
{
  struct Case
  {
    Body B;
    double lo, hi;
  };
  std::vector<Case> cases = {{namedBox({0, 0}, {2, 1000}, "box2x1000"), 0, 1000},
                             {namedBox({0, 0}, {1, 10}, "box1x10"), 0, 10},
                             {stretchedCube(5, 10, false), -10, 10},
                             {stretchedCube(5, 10, true), 1, 0},
                             {stretchedCube(10, 100, true), 1, 0}};
  for (auto& [B, lo, hi] : cases)
  {
    const int n = int(B.A.cols());
    const VT c = chebyshevCentre(B);
    for (SamplerWalk w : kWalks)
    {
      auto pts = roundedProposals(B, kN, 1234, w);
      const KS ks = gaugeKS(B, c, pts);
      KS mk;
      if (lo < hi) mk = marginalKS(pts, n - 1, lo, hi);
      std::printf("[rounded L=%u] %-20s %-14s gauge D=%.4f p=%.2e", ttcWalkLength(n),
                  B.name.c_str(), walkName(w), ks.D, ks.p);
      if (lo < hi) std::printf("  long-axis D=%.4f p=%.2e", mk.D, mk.p);
      std::printf("\n");
      if (w == SamplerWalk::Billiard)
      {
        CHECK_MESSAGE(ks.p > kAlpha, B.name << " gauge D=" << ks.D);
        if (lo < hi) CHECK_MESSAGE(mk.p > kAlpha, B.name << " long-axis D=" << mk.D);
      }
    }
  }
}

// Independence of consecutive samples of one call (one RNG stream).  Under
// independence r_lag * sqrt(N) ~ N(0,1); reject |z| > 3.29 (alpha = 0.001
// two-sided per lag).  Uses the gauge and the first coordinate, lags 1..3.
TEST_CASE("mixing: consecutive samples of one call are uncorrelated")
{
  Body C5 = cube(5);
  C5.name = "cube5";
  for (const Body& B : {C5, stretchedCube(5, 10, true)})
  {
    const VT c = chebyshevCentre(B);
    for (SamplerWalk w : kWalks)
    {
      auto pts = ttcProposals(B.A, B.b, kN, 99, w);
      const auto g = gauges(B, c, pts);
      std::vector<double> x0;
      for (const auto& p : pts) x0.push_back(p(0));
      for (int lag = 1; lag <= 3; ++lag)
      {
        const double zg = lagCorr(g, g, lag) * std::sqrt(double(kN - lag));
        const double zx = lagCorr(x0, x0, lag) * std::sqrt(double(kN - lag));
        if (lag == 1)
          std::printf("[lag] %-20s %-14s z(gauge,1)=%+.2f z(x0,1)=%+.2f\n", B.name.c_str(),
                      walkName(w), zg, zx);
        CHECK_MESSAGE(std::abs(zg) < 3.29, B.name << " " << std::string(walkName(w)) << " lag " << lag);
        CHECK_MESSAGE(std::abs(zx) < 3.29, B.name << " " << std::string(walkName(w)) << " lag " << lag);
      }
    }
  }
}

// Per-polytope seed streams: the same seed reproduces the sample set; two
// different streams (lra_volume's splitmix seeds, and adjacent raw seeds as
// in tests' ttcLatticeSamples) give sample sets whose index-aligned pairs are
// uncorrelated (|z| < 3.29) and share no point.
TEST_CASE("mixing: per-polytope seed streams give independent sample sets")
{
  const Body B = cube(4);
  const VT c = chebyshevCentre(B);
  auto same1 = ttcProposals(B.A, B.b, 200, 17);
  auto same2 = ttcProposals(B.A, B.b, 200, 17);
  bool identical = true;
  for (int k = 0; k < 200; ++k) identical &= same1[k] == same2[k];
  CHECK(identical);

  const std::pair<unsigned, unsigned> pairs[] = {
      {streamSeed(0, 1), streamSeed(0, 2)}, {17u, 18u}, {17u * 2654435761u + 1u, 17u * 2654435761u + 2u}};
  for (auto [s1, s2] : pairs)
  {
    for (SamplerWalk w : {SamplerWalk::Billiard, SamplerWalk::Ball, SamplerWalk::CDHR})
    {
      auto P1 = ttcProposals(B.A, B.b, kN, s1, w);
      auto P2 = ttcProposals(B.A, B.b, kN, s2, w);
      const auto g1 = gauges(B, c, P1), g2 = gauges(B, c, P2);
      std::vector<double> x1, x2;
      for (int k = 0; k < kN; ++k) x1.push_back(P1[k](0)), x2.push_back(P2[k](0));
      int shared = 0;
      for (int k = 0; k < kN; ++k) shared += (P1[k] == P2[k]);
      for (int lag = 0; lag <= 1; ++lag)
      {
        const double zg = lagCorr(g1, g2, lag) * std::sqrt(double(kN - lag));
        const double zx = lagCorr(x1, x2, lag) * std::sqrt(double(kN - lag));
        if (lag == 0)
          std::printf("[seeds %u,%u] %-14s z(gauge)=%+.2f z(x0)=%+.2f shared=%d\n", s1, s2,
                      walkName(w), zg, zx, shared);
        CHECK(std::abs(zg) < 3.29);
        CHECK(std::abs(zx) < 3.29);
      }
      CHECK(shared == 0);
    }
  }
}

// ======================================================================== slow

// Smallest L on a 1.5x geometric grid such that the gauge KS (and the
// long-axis KS for axis-aligned boxes) passes at alpha = 0.001 at L and at the
// next grid point.  Expected for the paper's theory (from a point mass, no
// warm start, so theory is an upper-bound guide only): ball walk and
// hit-and-run ~ n^2 (R/r)^2 steps; billiard empirically ~ O(n) on round bodies.
TEST_CASE("mixing: smallest passing walk length per walk and body" * doctest::test_suite("slow"))
{
  struct Case
  {
    Body B;
    double lo, hi;
    bool rounded = false;
  };
  std::vector<Case> cases;
  for (const Body& B : roundBodies()) cases.push_back({B, 1, 0});
  {
    Body C20 = cube(20);
    C20.name = "cube20";
    cases.push_back({C20, 1, 0});
  }
  cases.push_back({namedBox({0, 0}, {1, 10}, "box1x10"), 0, 10});
  cases.push_back({stretchedCube(5, 10, false), -10, 10});
  cases.push_back({stretchedCube(5, 10, true), 1, 0});
  cases.push_back({stretchedCube(5, 3, true), 1, 0});
  cases.push_back({stretchedCube(5, 30, true), 1, 0});
  // Candidate fix: the same elongated bodies sampled in rounded coordinates.
  cases.push_back({namedBox({0, 0}, {2, 1000}, "box2x1000"), 0, 1000, true});
  cases.push_back({namedBox({0, 0}, {1, 10}, "box1x10"), 0, 10, true});
  cases.push_back({stretchedCube(5, 10, true), 1, 0, true});
  cases.push_back({stretchedCube(5, 30, true), 1, 0, true});
  cases.push_back({stretchedCube(10, 100, true), 1, 0, true});
  const unsigned cap = 3000;
  for (auto& [B, lo, hi, rounded] : cases)
  {
    const int n = int(B.A.cols());
    const VT c = chebyshevCentre(B);
    const std::string label = B.name + (rounded ? "(R)" : "");
    for (SamplerWalk w : kWalks)
    {
      unsigned found = 0, passStreak = 0, firstPass = 0;
      double lastD = 1;
      for (double Ld = 1; Ld <= cap; Ld = std::max(Ld + 1, std::floor(Ld * 1.5)))
      {
        const unsigned L = unsigned(Ld);
        auto pts = rounded ? roundedProposals(B, kN, 4321 + L, w, int(L))
                           : ttcProposals(B.A, B.b, kN, 4321 + L, w, int(L));
        KS ks = gaugeKS(B, c, pts);
        bool pass = ks.p > kAlpha;
        if (lo < hi) pass = pass && marginalKS(pts, n - 1, lo, hi).p > kAlpha;
        lastD = ks.D;
        if (pass)
        {
          if (passStreak++ == 0) firstPass = L;
          if (passStreak == 2) { found = firstPass; break; }
        }
        else
          passStreak = 0;
      }
      if (found)
        std::printf("[Lmin] %-22s n=%-2d %-14s Lmin=%-5u default=%-4u ratio=%.2f\n",
                    label.c_str(), n, walkName(w), found, ttcWalkLength(n),
                    double(found) / ttcWalkLength(n));
      else
        std::printf("[Lmin] %-22s n=%-2d %-14s Lmin>%u (D at cap %.4f) default=%u\n",
                    label.c_str(), n, walkName(w), cap, lastD, ttcWalkLength(n));
      std::fflush(stdout);
    }
  }
}

// One chain per polytope instead of N restarts: burn-in then thinning.
// Measurement only (no assertions; informs the design recommendation):
// gauge KS and lag-1 z of the thinned chain, and its total step count
// against the restart design's N * (10 + n^2).
TEST_CASE("mixing: one long chain with burn-in and thinning" * doctest::test_suite("slow"))
{
  struct Case
  {
    Body B;
    unsigned burnIn, thin;
  };
  Body C10 = cube(10);
  C10.name = "cube10";
  for (const Body& B : {C10, stretchedCube(5, 10, true)})
  {
    const int n = int(B.A.cols());
    const VT c = chebyshevCentre(B);
    for (SamplerWalk w : kWalks)
    {
      for (unsigned thin : {1u, 2u, 5u, 10u})
      {
        const unsigned burnIn = 10 * ttcWalkLength(n);
        auto pts = chainSamples(B, w, burnIn, thin, kN, 777);
        const KS ks = gaugeKS(B, c, pts);
        const auto g = gauges(B, c, pts);
        const double z = lagCorr(g, g, 1) * std::sqrt(double(kN - 1));
        std::printf("[chain] %-20s %-14s burn=%u thin=%-2u steps=%u (restart %u) "
                    "gauge D=%.4f p=%.2e lag1 z=%+.2f\n",
                    B.name.c_str(), walkName(w), burnIn, thin, burnIn + thin * kN,
                    kN * ttcWalkLength(n), ks.D, ks.p, z);
        std::fflush(stdout);
      }
    }
  }
}
