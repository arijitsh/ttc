// GaussianBallWalk as an MCMC sampler (the walk ttc's LRA engine hands to
// volume_cooling_gaussians), checked against analytic targets.
//
// The oracle shifts P so that its Chebyshev centre is the origin and then runs
// a ball walk whose stationary density is  pi_a(x) ∝ exp(-a |x|^2) 1_P(x).
// For axis-aligned boxes pi_a factorises into 1-D truncated normals
// (sd 1/sqrt(2a)), so every marginal moment / CDF is known in closed form.
//
// What is checked here:
//   * GetPointInDsphere (sampling/sphere.hpp) is uniform in the d-ball.
//   * HPolytope::is_in convention (-1 = inside incl. boundary, 0 = outside),
//     which the walk's acceptance test `P.is_in(y) == -1` relies on.
//   * the step size delta = 4 r / sqrt(max(1,a) n) of the Walk ctor, and the
//     (by-value, hence no-op) update_delta helper in volume_cooling_gaussians.
//   * stationary marginals (mean / variance / KS) for a in {large, medium, 0}.
//   * reversibility (detailed balance) via an antisymmetric pair statistic.
//   * mixing: walk_length 10+n^2 decorrelates on well-rounded bodies, but not
//     on the user's long box [0,2]x[0,1000]; diffusion-time measurement.
//   * the cooling schedule the oracle builds from the chain's samples: the
//     exact relative variance of every ratio estimator (computed analytically
//     for boxes) must stay bounded.  On the long box the oracle (seed 3)
//     jumps a = 0.517 -> 0, whose relative variance is ~10^224000; even its
//     k=1 halvings have ~10^2.7 because for n=2 ratio=1-1/n gives 2a_{i+1}-a_i=0.
#include "doctest.h"
#include "common.hpp"

#include <cstdio>
#include <iostream>

using namespace vt;

namespace
{
using Rng = TtcSampleRNG;
using GBWalk = GaussianBallWalk::Walk<HPoly, Rng>;

// ------------------------------------------------------------ 1-D analytics
double Phi(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
double phi(double x) { return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI); }

// Density ∝ exp(-a x^2) on [l,u] (a = 0: uniform).
struct TruncNormal
{
  double a, l, u;
  double mean() const
  {
    if (a <= 0) return 0.5 * (l + u);
    const double s = 1.0 / std::sqrt(2 * a), al = l / s, be = u / s;
    const double Z = Phi(be) - Phi(al);
    return s * (phi(al) - phi(be)) / Z;
  }
  double var() const
  {
    if (a <= 0) return (u - l) * (u - l) / 12.0;
    const double s = 1.0 / std::sqrt(2 * a), al = l / s, be = u / s;
    const double Z = Phi(be) - Phi(al);
    const double t = (phi(al) - phi(be)) / Z;
    return s * s * (1 + (al * phi(al) - be * phi(be)) / Z - t * t);
  }
  double cdf(double x) const
  {
    if (x <= l) return 0;
    if (x >= u) return 1;
    if (a <= 0) return (x - l) / (u - l);
    const double s = 1.0 / std::sqrt(2 * a);
    return (Phi(x / s) - Phi(l / s)) / (Phi(u / s) - Phi(l / s));
  }
};

// log of  int_l^u exp(-c y^2) dy  for any real c (log-sum-exp trapezoid).
double logZ1D(double c, double l, double u)
{
  const int N = 400000;
  const double h = (u - l) / N;
  double mx = -INFINITY;
  for (int k = 0; k <= N; ++k)
  {
    const double y = l + k * h;
    mx = std::max(mx, -c * y * y);
  }
  double s = 0;
  for (int k = 0; k <= N; ++k)
  {
    const double y = l + k * h;
    s += (k == 0 || k == N ? 0.5 : 1.0) * std::exp(-c * y * y - mx);
  }
  return mx + std::log(s * h);
}

// --------------------------------------------------------------- statistics
double mean(const std::vector<double>& x)
{
  return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

// Batch-means standard error of the mean of a (correlated) chain.
double batchSE(const std::vector<double>& x, int nb = 100)
{
  const size_t bs = x.size() / nb;
  std::vector<double> m(nb);
  for (int b = 0; b < nb; ++b)
    m[b] = std::accumulate(x.begin() + b * bs, x.begin() + (b + 1) * bs, 0.0) / bs;
  const double mm = mean(m);
  double v = 0;
  for (double t : m) v += (t - mm) * (t - mm);
  return std::sqrt(v / (nb - 1) / nb);
}

double lagAutocorr(const std::vector<double>& x, size_t lag)
{
  const double m = mean(x);
  double num = 0, den = 0;
  for (size_t i = 0; i < x.size(); ++i)
  {
    den += (x[i] - m) * (x[i] - m);
    if (i + lag < x.size()) num += (x[i] - m) * (x[i + lag] - m);
  }
  return num / den;
}

template <class CDF>
double ksStat(std::vector<double> x, CDF F)
{
  std::sort(x.begin(), x.end());
  const double n = x.size();
  double D = 0;
  for (size_t i = 0; i < x.size(); ++i)
  {
    const double f = F(x[i]);
    D = std::max({D, std::abs(f - i / n), std::abs((i + 1) / n - f)});
  }
  return D;
}

// Box [lo,hi] (must contain the origin) with the inner ball set to
// (origin, min distance to a facet) -- the state the oracle is in after its
// Chebyshev shift.
HPoly centredBox(const std::vector<double>& lo, const std::vector<double>& hi, double* rOut = nullptr)
{
  HPoly P = toHPoly(box(lo, hi));
  double r = INFINITY;
  for (size_t i = 0; i < lo.size(); ++i) r = std::min({r, -lo[i], hi[i]});
  setBall(P, VT::Zero(lo.size()), r);
  if (rOut) *rOut = r;
  return P;
}

// The polytope exactly as the oracle sees it: volesti inner ball, then shifted
// so the Chebyshev centre is the origin.  Returns the shifted box bounds.
HPoly oracleShiftedBox(const Body& B, std::vector<double>& lo, std::vector<double>& hi, double& r)
{
  HPoly P = toHPoly(B);
  auto ib = P.ComputeInnerBall();
  P.shift(ib.first.getCoefficients());  // shift() keeps _inner_ball.second
  r = ib.second;
  const int n = B.A.cols();
  lo.assign(n, 0);
  hi.assign(n, 0);
  for (int i = 0; i < n; ++i)  // box(): rows 2i = +e_i <= hi, 2i+1 = -e_i <= -lo
  {
    hi[i] = B.b(2 * i) - ib.first[i];
    lo[i] = -B.b(2 * i + 1) - ib.first[i];
  }
  return P;
}

double maxStep(GBWalk& w, HPoly& P, double a, int n, int steps, Rng& rng, int* moves = nullptr)
{
  Point p(n);
  double mx = 0;
  int mv = 0;
  for (int t = 0; t < steps; ++t)
  {
    Point q = p;
    w.apply(P, p, a, 1u, rng);
    const double d = std::sqrt((p - q).squared_length());
    if (d > 0) ++mv;
    mx = std::max(mx, d);
  }
  if (moves) *moves = mv;
  return mx;
}

}  // namespace

// ===========================================================================
TEST_CASE("sphere.hpp: GetPointInDsphere is uniform in the d-ball")
{
  Rng rng(1);
  rng.set_seed(12345);
  const double R = 2.5;
  const int N = 100000;
  for (int n : {1, 2, 3, 5, 10})
  {
    std::vector<double> rad(N), x0(N), x0sq(N), x01(N);
    bool inside = true;
    for (int k = 0; k < N; ++k)
    {
      Point p = GetPointInDsphere<Point>::apply(n, R, rng);
      const double r = std::sqrt(p.squared_length());
      inside &= r <= R * (1 + 1e-12);
      rad[k] = r;
      x0[k] = p[0];
      x0sq[k] = p[0] * p[0];
      x01[k] = n > 1 ? p[0] * p[1] : 0.0;
    }
    // |x| has CDF (t/R)^n.
    const double D = ksStat(rad, [&](double t) { return std::pow(std::min(t / R, 1.0), n); });
    // Uniform ball: E x_i^2 = R^2/(n+2), E x_i^4 = 3R^4/((n+2)(n+4)), E x_i^2 x_j^2 = R^4/((n+2)(n+4)).
    const double m2 = R * R / (n + 2), m4 = 3 * std::pow(R, 4) / ((n + 2) * (n + 4));
    const double seMean = std::sqrt(m2 / N), seVar = std::sqrt((m4 - m2 * m2) / N);
    const double seCov = std::sqrt(std::pow(R, 4) / ((n + 2) * (n + 4)) / N);
    std::printf("[dsphere] n=%2d KS=%.4f (crit %.4f) mean=%+.4f var=%.4f/%.4f cov01=%+.4f\n", n, D,
                1.95 / std::sqrt(N), mean(x0), mean(x0sq), m2, mean(x01));
    CHECK(inside);
    CHECK(D < 1.95 / std::sqrt(double(N)));  // alpha ~ 0.1 %
    CHECK(std::abs(mean(x0)) < 5 * seMean);
    CHECK(std::abs(mean(x0sq) - m2) < 5 * seVar);
    if (n > 1) CHECK(std::abs(mean(x01)) < 5 * seCov);
  }
}

TEST_CASE("hpolytope.h: is_in returns -1 inside/on boundary, 0 outside")
{
  HPoly P = toHPoly(box({-1, -2}, {1, 2}));
  P.normalize();
  CHECK(P.is_in(Point(VT::Zero(2))) == -1);
  VT v(2);
  v << 1.0, 2.0;  // vertex: on the boundary, b - Ax = 0 is not < -0
  CHECK(P.is_in(Point(v)) == -1);
  v << 1.0 + 1e-9, 0.0;
  CHECK(P.is_in(Point(v)) == 0);
  v << 0.0, -2.5;
  CHECK(P.is_in(Point(v)) == 0);
}

TEST_CASE("GaussianBallWalk: step size delta = 4r/sqrt(max(1,a) n); update_delta is a harmless no-op")
{
  // Huge box so nothing is rejected by the boundary; fake inner radius r.
  for (int n : {2, 3, 6})
  {
    HPoly P = toHPoly(cube(n, -1e6, 1e6));
    const double r = 0.7;
    setBall(P, VT::Zero(n), r);
    Rng rng(n);
    rng.set_seed(7);
    for (double a : {0.0, 0.5, 1.0, 4.0, 100.0})
    {
      const double want = 4 * r / std::sqrt(std::max(1.0, a) * n);
      CHECK(GBWalk::compute_delta(P, a) == doctest::Approx(want).epsilon(1e-14));
      GBWalk w(P, Point(n), a, rng);
      // volume_cooling_gaussians.hpp: update_delta<WalkType>::apply takes the
      // walk BY VALUE, so this call cannot change w.  It is harmless only
      // because the requested delta equals what the ctor already set.
      update_delta<GBWalk>::apply(w, 1e-3 * want);
      int moves = 0;
      const double mx = maxStep(w, P, a, n, 40000, rng, &moves);
      std::printf("[delta] n=%d a=%-5g delta=%.4f maxStep=%.4f moved=%d/40000\n", n, a, want, mx,
                  moves);
      CHECK(mx <= want * (1 + 1e-12));
      CHECK(mx >= want * 0.97);  // proposal radius really is delta (not 1e-3*delta)
    }
  }
}

TEST_CASE("GaussianBallWalk: stationary marginals are the truncated normals of exp(-a|x|^2) on a box")
{
  struct Case
  {
    std::vector<double> lo, hi;
    double a;
  };
  const std::vector<Case> cases = {
      {{-1, -1}, {1, 1}, 8.0},          {{-1, -1}, {1, 1}, 0.5},
      {{-1, -1}, {1, 1}, 0.0},          {{-1, -1}, {3, 1}, 2.0},
      {{-1, -1}, {3, 1}, 0.25},         {{-1, -1}, {3, 1}, 0.0},
      {{-1, -1, -1}, {1, 2, 1.5}, 1.0}, {{-1, -1, -1}, {1, 2, 1.5}, 0.0},
  };
  const int burn = 5000, M = 400000, thin = 100;
  Rng rng(3);
  rng.set_seed(2024);
  for (const Case& c : cases)
  {
    const int n = c.lo.size();
    double r;
    HPoly P = centredBox(c.lo, c.hi, &r);
    GBWalk w(P, Point(n), c.a, rng);
    Point p(n);
    for (int t = 0; t < burn; ++t) w.apply(P, p, c.a, 1u, rng);
    std::vector<std::vector<double>> xs(n, std::vector<double>(M));
    bool alwaysInside = true;
    int moved = 0;
    for (int t = 0; t < M; ++t)
    {
      Point q = p;
      w.apply(P, p, c.a, 1u, rng);
      moved += (p - q).squared_length() > 0;
      alwaysInside &= P.is_in(p) == -1;
      for (int i = 0; i < n; ++i) xs[i][t] = p[i];
    }
    CHECK(alwaysInside);
    for (int i = 0; i < n; ++i)
    {
      TruncNormal T{c.a, c.lo[i], c.hi[i]};
      const auto& x = xs[i];
      std::vector<double> dev2(M), thinned;
      const double mu = T.mean(), var = T.var();
      for (int t = 0; t < M; ++t) dev2[t] = (x[t] - mu) * (x[t] - mu);
      for (int t = 0; t < M; t += thin) thinned.push_back(x[t]);
      const double m = mean(x), v = mean(dev2);
      const double seM = batchSE(x), seV = batchSE(dev2);
      const double D = ksStat(thinned, [&](double y) { return T.cdf(y); });
      const double tau = (seM * seM * M) / var;  // integrated autocorrelation time (steps)
      std::printf("[stationary] box n=%d a=%-4g coord %d: mean %+.4f/%+.4f (se %.4f) "
                  "var %.4f/%.4f (se %.4f) KS %.4f tau~%.1f acc=%.2f\n",
                  n, c.a, i, m, mu, seM, v, var, seV, D, tau, double(moved) / M);
      CHECK(std::abs(m - mu) < 5 * seM + 1e-3);
      CHECK(std::abs(v - var) < 5 * seV + 1e-3);
      CHECK(D < 2.2 / std::sqrt(double(thinned.size())));  // alpha ~ 0.01 %, thinned >> tau
    }
  }
}

TEST_CASE("GaussianBallWalk: chain is reversible (antisymmetric pair statistic vanishes)")
{
  // Under a reversible kernel (X_t, X_{t+1}) is exchangeable, so
  // E[g(X_t) h(X_{t+1}) - h(X_t) g(X_{t+1})] = 0 for all g, h.
  Rng rng(2);
  rng.set_seed(99);
  for (double a : {0.25, 2.0})
  {
    HPoly P = centredBox({-1, -1}, {3, 1});
    GBWalk w(P, Point(2), a, rng);
    Point p(2);
    for (int t = 0; t < 5000; ++t) w.apply(P, p, a, 1u, rng);
    const int M = 400000;
    std::vector<double> d(M);
    auto g = [](const Point& x) { return x[0]; };
    auto h = [](const Point& x) { return x[0] * x[0] + x[1]; };
    for (int t = 0; t < M; ++t)
    {
      Point q = p;
      w.apply(P, p, a, 1u, rng);
      d[t] = g(q) * h(p) - h(q) * g(p);
    }
    const double m = mean(d), se = batchSE(d);
    std::printf("[reversible] a=%g  E[antisym]=%+.5f se=%.5f\n", a, m, se);
    CHECK(std::abs(m) < 5 * se + 1e-4);
  }
}

TEST_CASE("mixing: walk_length 10+n^2 decorrelates successive samples on well-rounded bodies")
{
  // The oracle's convergence test (sliding-window of running means) is only
  // meaningful if successive samples (walk_length steps apart) are not
  // nearly identical.  Start in stationarity (uniform start on the cube).
  Rng rng(2);
  rng.set_seed(5);
  std::mt19937 g(5);
  for (int n = 2; n <= 6; ++n)
    for (double a : {0.0, 1.0})
    {
      HPoly P = centredBox(std::vector<double>(n, -1), std::vector<double>(n, 1));
      GBWalk w(P, Point(n), a, rng);
      Point p(n);
      for (int t = 0; t < 2000; ++t) w.apply(P, p, a, 1u, rng);
      const unsigned wl = ttcWalkLength(n);
      std::vector<double> y(4000);
      for (auto& v : y)
      {
        w.apply(P, p, a, wl, rng);
        v = p[0];
      }
      const double rho = lagAutocorr(y, 1);
      std::printf("[mix-cube] n=%d a=%g walk_length=%u lag-1 autocorr of x0 = %.3f\n", n, a, wl, rho);
      CHECK(rho < 0.5);
    }
}

TEST_CASE("mixing: long box [0,2]x[0,1000] as the oracle sees it (diffusion time vs walk_length 14)"
          * doctest::may_fail())
{
  // The oracle picks the volesti Chebyshev centre, which for this box is
  // (1, 999) -- the very END of the strip -- so after the shift the chain
  // starts at the origin of [-1,1] x [-999,1] and must diffuse ~1000 units
  // with proposals of radius delta = 4r/sqrt(2) = 2.83 (r = 1).
  Body B = box({0, 0}, {2, 1000});
  std::vector<double> lo, hi;
  double r;
  HPoly P = oracleShiftedBox(B, lo, hi, r);
  std::printf("[strip] volesti Chebyshev ball r=%.3f, shifted box y in [%g, %g]\n", r, lo[1], hi[1]);
  const int n = 2;
  const double L = hi[1] - lo[1];
  Rng rng(2);
  rng.set_seed(11);
  GBWalk w(P, Point(n), 0.0, rng);

  // (1) Per-step diffusion constant of y in the uniform (a=0) target.
  Point p(n);
  std::vector<double> dy2;
  for (int t = 0; t < 200000; ++t)
  {
    Point q = p;
    w.apply(P, p, 0.0, 1u, rng);
    dy2.push_back((p[1] - q[1]) * (p[1] - q[1]));
  }
  const double s2 = mean(dy2);
  // Reflected random walk on [0,L]: slowest mode relaxes in L^2/(pi^2 s2) steps.
  const double tauRelax = L * L / (M_PI * M_PI * s2);
  const double crossSteps = L * L / s2;  // ~ steps for one end-to-end traversal
  std::printf("[strip] E[dy^2] per step = %.3f -> relaxation ~%.3g steps, traversal ~%.3g steps "
              "(= %.3g samples at walk_length 14)\n",
              s2, tauRelax, crossSteps, crossSteps / 14);

  // (2) Lag-1 autocorrelation of y at the oracle's walk length, in stationarity.
  std::mt19937 g(3);
  std::uniform_real_distribution<double> Uy(lo[1], hi[1]);
  p = Point(n);
  p.set_coord(1, Uy(g));
  std::vector<double> y(4000);
  for (auto& v : y)
  {
    w.apply(P, p, 0.0, ttcWalkLength(n), rng);
    v = p[1];
  }
  const double rho = lagAutocorr(y, 1);
  std::printf("[strip] lag-1 autocorr of y at walk_length 14 = %.5f (cube: ~0)\n", rho);

  // (3) What the oracle's chain actually does in a schedule-test phase:
  //     totalSteps = 150/((1-frac)*err)+1 samples of walk_length 14 from the origin.
  const unsigned totalSteps = unsigned(150 / (0.9 * 0.8 / 12)) + 1;
  p = Point(n);
  std::vector<double> ys(totalSteps);
  for (auto& v : ys)
  {
    w.apply(P, p, 0.0, ttcWalkLength(n), rng);
    v = p[1];
  }
  const TruncNormal Ty{0, lo[1], hi[1]};
  const double my = mean(ys);
  double vy = 0;
  for (double v : ys) vy += (v - my) * (v - my);
  vy /= ys.size();
  std::printf("[strip] %u samples x 14 steps from the start (y=0): mean y = %.1f (target %.1f), "
              "sample var / target var = %.3f, min y = %.1f\n",
              totalSteps, my, Ty.mean(), vy / Ty.var(), *std::min_element(ys.begin(), ys.end()));

  // Correct behaviour for the sampler as the oracle configures it: samples
  // walk_length apart should be weakly correlated, and one phase's worth of
  // samples should cover the body.  Both fail by orders of magnitude here; the
  // walk itself is correct (see the stationary/reversible tests) -- the body
  // is just not rounded, so delta ~ r is tiny relative to the diameter.
  CHECK(rho < 0.5);                                     // ~0.999
  CHECK(std::abs(my - Ty.mean()) < 0.1 * L);            // ~ -60 vs -499
  CHECK(vy / Ty.var() == doctest::Approx(1.0).epsilon(0.3));  // ~0.02
}

TEST_CASE("cooling schedule built from chain samples: exact relative variance of every ratio is bounded")
{
  // For a box, Z(c) = int_P exp(-c|x|^2) is a product of 1-D integrals, so the
  // true relative variance of the i-th ratio estimator
  //     f = exp(-(a_{i+1}-a_i)|x|^2),  x ~ pi_{a_i}
  // is  Z(2a_{i+1}-a_i) Z(a_i) / Z(a_{i+1})^2 - 1.  Cousins-Vempala choose
  // a_{i+1} so that this is <= C-1 = 1 (C = 2); the schedule is built from the
  // walk's samples (get_next_gaussian), so a badly mixed chain shows up here.
  // The schedule is regenerated with the oracle's exact RNG stream (fixed
  // seed 3; compute_annealing_schedule is the first consumer of the RNG).
  for (const Body& B : {cube(2), cube(4), box({0, 0}, {2, 1000})})
  {
    const int n = B.A.cols();
    std::vector<double> lo, hi;
    double r;
    HPoly P = oracleShiftedBox(B, lo, hi, r);
    TtcVolumeRNG rng(n);
    gaussian_annealing_parameters<double> par(n);
    std::vector<double> a;
    using VW = GaussianBallWalk::Walk<HPoly, TtcVolumeRNG>;
    compute_annealing_schedule<VW, GaussianRandomPointGenerator<VW>>(
        P, par.ratio, par.C, par.frac, par.N, ttcWalkLength(n), r, 0.8 / 12, a, rng);
    auto logZ = [&](double c) {
      double s = 0;
      for (int i = 0; i < n; ++i) s += logZ1D(c, lo[i], hi[i]);
      return s;
    };
    double worst = -INFINITY;  // log10(1 + relvar)
    std::string sched;
    for (size_t i = 0; i + 1 < a.size(); ++i)
    {
      const double l10 = (logZ(2 * a[i + 1] - a[i]) + logZ(a[i]) - 2 * logZ(a[i + 1])) / std::log(10.0);
      worst = std::max(worst, l10);
      char buf[96];
      std::snprintf(buf, sizeof buf, " %.3g[%.3g]", a[i], l10);
      sched += buf;
    }
    sched += " " + std::to_string(a.back());
    std::printf("[schedule] %s: a[log10(E f^2/(E f)^2) of step a_i->a_i+1] =%s\n"
                "           worst = %.4g\n",
                B.name.c_str(), sched.c_str(), worst);
    CHECK(worst < std::log10(1 + 3.0));  // relvar <= 3 (C-1 = 1, with slack)
  }
}
