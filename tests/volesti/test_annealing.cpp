// Annealing-schedule internals of volume_cooling_gaussians (Cousins-Vempala
// Gaussian cooling), in exactly ttc's configuration.
//
// The volume estimate is the telescoping product
//
//     vol(P) = Z(0) = (pi/a_0)^{n/2} * prod_{i=0}^{k-1} Z(a_{i+1}) / Z(a_i),
//     Z(a)   = int_P exp(-a |x|^2) dx,    a_0 > a_1 > ... > a_k = 0,
//
// where phase i estimates R_i = Z(a_{i+1})/Z(a_i) = E_{mu_i}[f_{i+1}/f_i] from a
// Gaussian ball walk targeting mu_i ~ exp(-a_i|x|^2) on P.  For an axis-aligned
// box (after shifting the Chebyshev center to the origin) Z(a) factorizes into
// 1-D integrals that we evaluate in closed form (erf, a > 0) or by stable
// quadrature (a <= 0), so every per-phase ratio the oracle produces can be
// compared with its exact value.
//
// traceOracle() below is a line-by-line replica of volume_cooling_gaussians
// (src/volesti/include/volume/volume_cooling_gaussians.hpp:292) that calls the
// library's compute_annealing_schedule and then records the per-phase ratios.
// It consumes the RNG identically, so it reproduces the oracle bit for bit
// (checked by the first test case).
#include "doctest.h"
#include "common.hpp"

#include <array>
#include <cstdio>
#include <iostream>
#include <list>
#include <memory>

using namespace vt;

namespace
{

using WalkT = GaussianBallWalk::Walk<HPoly, TtcVolumeRNG>;
using PointGen = GaussianRandomPointGenerator<WalkT>;

constexpr double kTtcError = 0.8 / 12.0;  // ttc's oracle error target (eps/12)

struct Trace
{
  VT center;                  // inner-ball center used by the oracle
  double radius = 0;          // inner-ball radius used by the oracle
  double scale = 1;           // oracle rescales the body by 1/scale (unit inner radius)
  MT A;                       // shifted, rescaled, row-normalized polytope (center at 0)
  VT b;
  std::vector<double> a;      // a_0 > a_1 > ... > a_k
  std::vector<double> ratio;  // per-phase estimate of Z(a_{i+1})/Z(a_i)
  std::vector<double> steps;  // walk steps used per phase
  double logVol0 = 0;         // log (pi/a_0)^{n/2}
  double vol = 0;             // final estimate (== oracle output)
};

// Replica of volume_cooling_gaussians (volume_cooling_gaussians.hpp:292-455).
template <typename RNG>
Trace traceOracle(HPoly Pin, RNG& rng, double error, unsigned walk_length)
{
  using NT = double;
  Trace T;
  HPoly P(Pin);
  const unsigned n = P.dimension();
  const unsigned m = P.num_of_hyperplanes();
  gaussian_annealing_parameters<NT> parameters(n);
  auto InnerBall = P.ComputeInnerBall();
  T.center = InnerBall.first.getCoefficients();
  T.radius = InnerBall.second;
  if (InnerBall.second < 0.0) { T.vol = -1; return T; }
  NT radius = InnerBall.second;
  P.shift(InnerBall.first.getCoefficients());
  // Rescale to unit inner radius, as volume_cooling_gaussians does.
  const NT scale = radius;
  P.normalize();
  P.set_vec(P.get_vec() / scale);
  radius = NT(1);
  P.set_InnerBall(std::pair<Point, NT>(Point(n), radius));
  T.scale = scale;
  T.A = P.get_mat();
  T.b = P.get_vec();

  std::vector<NT> a_vals;
  compute_annealing_schedule<WalkT, PointGen>(P, parameters.ratio, parameters.C, parameters.frac,
                                              parameters.N, walk_length, radius, error,
                                              a_vals, rng);
  T.a = a_vals;

  const unsigned W = parameters.W;
  const unsigned mm = a_vals.size() - 1;
  std::vector<NT> last_W2(W, 0), fn(mm, 0), its(mm, 0);
  NT vol = std::pow(M_PI / a_vals[0], NT(n) / 2.0);
  T.logVol0 = 0.5 * n * std::log(M_PI / a_vals[0]);
  Point p(n);
  for (unsigned i = 0; i < mm; ++i)
  {
    bool done = false;
    NT curr_eps = error / std::sqrt(NT(mm));
    NT min_val = std::numeric_limits<NT>::min();
    NT max_val = std::numeric_limits<NT>::max();
    unsigned min_index = W - 1, max_index = W - 1, index = 0, min_steps = 0;
    std::vector<NT> last_W = last_W2;
    WalkT walk(P, p, a_vals[i], rng);
    update_delta<WalkT>::apply(walk, 4.0 * radius / std::sqrt(std::max(NT(1.0), a_vals[i]) * NT(n)));
    while (!done || its[i] < min_steps)
    {
      walk.apply(P, p, a_vals[i], walk_length, rng);
      its[i] += 1.0;
      fn[i] += eval_exp(p, a_vals[i + 1]) / eval_exp(p, a_vals[i]);
      NT val = fn[i] / its[i];
      last_W[index] = val;
      if (val <= min_val) { min_val = val; min_index = index; }
      else if (min_index == index)
      {
        auto it = std::min_element(last_W.begin(), last_W.end());
        min_val = *it;
        min_index = std::distance(last_W.begin(), it);
      }
      if (val >= max_val) { max_val = val; max_index = index; }
      else if (max_index == index)
      {
        auto it = std::max_element(last_W.begin(), last_W.end());
        max_val = *it;
        max_index = std::distance(last_W.begin(), it);
      }
      if ((max_val - min_val) / max_val <= curr_eps / 2.0) done = true;
      index = index % W + 1;
      if (index == W) index = 0;
    }
    T.ratio.push_back(fn[i] / its[i]);
    T.steps.push_back(its[i]);
    vol *= fn[i] / its[i];
  }
  T.vol = vol * std::pow(scale, NT(n));
  return T;
}

Trace traceTtc(const HPoly& P)
{
  TtcVolumeRNG rng(P.dimension());
  return traceOracle(P, rng, kTtcError, ttcWalkLength(P.dimension()));
}

// ------------------------------------------------- exact Gaussian integrals
double logAddExp(double x, double y)
{
  if (x < y) std::swap(x, y);
  return x + std::log1p(std::exp(y - x));
}

// log int_0^t exp(-a x^2) dx, t > 0, any sign of a.
double logHalfLine(double a, double t)
{
  if (a > 0)
  {
    const double s = std::sqrt(a);
    return std::log(0.5 * std::sqrt(M_PI / a) * std::erf(s * t));
  }
  if (a == 0) return std::log(t);
  // c = -a > 0: int_0^t e^{c x^2} = e^{c t^2} int_0^t e^{-c w (2t - w)} dw.
  const double c = -a;
  const double Wmax = std::min(t, 60.0 / (c * t));
  const int K = 20000;  // Simpson intervals (even)
  const double h = Wmax / K;
  double s = 0;
  for (int k = 0; k <= K; ++k)
  {
    const double w = k * h;
    const double f = std::exp(-c * w * (2 * t - w));
    s += f * (k == 0 || k == K ? 1 : (k % 2 ? 4 : 2));
  }
  return c * t * t + std::log(s * h / 3.0);
}

// For a box with center inside: A rows are +-e_i (normalized).  Returns per
// coordinate (lo_i, hi_i) with lo_i < 0 < hi_i; empty if P is not such a box.
bool boxBounds(const Trace& T, std::vector<double>& lo, std::vector<double>& hi)
{
  const int n = T.A.cols();
  lo.assign(n, -INFINITY);
  hi.assign(n, INFINITY);
  for (int r = 0; r < T.A.rows(); ++r)
  {
    int idx = -1;
    for (int j = 0; j < n; ++j)
      if (std::abs(T.A(r, j)) > 1e-12)
      {
        if (idx >= 0) return false;
        idx = j;
      }
    if (idx < 0) return false;
    const double v = T.b(r) / T.A(r, idx);
    if (T.A(r, idx) > 0) hi[idx] = std::min(hi[idx], v);
    else lo[idx] = std::max(lo[idx], v);
  }
  for (int j = 0; j < n; ++j)
    if (!(lo[j] < 0 && hi[j] > 0 && std::isfinite(lo[j]) && std::isfinite(hi[j]))) return false;
  return true;
}

// log Z(a) for the shifted box.
double logZbox(const std::vector<double>& lo, const std::vector<double>& hi, double a)
{
  double s = 0;
  for (size_t j = 0; j < lo.size(); ++j)
    s += logAddExp(logHalfLine(a, hi[j]), logHalfLine(a, -lo[j]));
  return s;
}

// Gaussian mass of N(0, I/(2a)) outside the shifted box.
double massOutsideBox(const std::vector<double>& lo, const std::vector<double>& hi, double a)
{
  double in = 1;
  const double s = std::sqrt(a);
  for (size_t j = 0; j < lo.size(); ++j) in *= 0.5 * (std::erf(s * hi[j]) + std::erf(-s * lo[j]));
  return 1 - in;
}

// Replica of compute_annealing_schedule (volume_cooling_gaussians.hpp:189-265)
// that calls the library get_first_gaussian / get_next_gaussian and records the
// stopping test of every iteration.  Consumes the RNG exactly like the library.
struct SchedStep
{
  double a = 0, next_a = 0, stopRatio = 0, maxR2walk = 0;
  int branch = 0;  // 0 = append next_a, 1 = next_a<=0: append 0, 2 = REPLACE a_it by 0
};
template <typename RNG>
std::vector<SchedStep> traceSchedule(HPoly Pin, RNG& rng, double error, unsigned walk_length,
                                     std::vector<double>& a_vals)
{
  using NT = double;
  HPoly P(Pin);
  const unsigned n = P.dimension();
  gaussian_annealing_parameters<NT> par(n);
  auto IB = P.ComputeInnerBall();
  P.shift(IB.first.getCoefficients());
  // Rescale to unit inner radius, as volume_cooling_gaussians does.
  P.normalize();
  P.set_vec(P.get_vec() / IB.second);
  const NT radius = NT(1);
  P.set_InnerBall(std::pair<Point, NT>(Point(n), radius));
  std::vector<SchedStep> out;
  a_vals.clear();
  get_first_gaussian(P, par.frac, radius, error, a_vals);
  const NT tol = 0.001;
  unsigned it = 0;
  const unsigned totalSteps = ((int)150 / ((1.0 - par.frac) * error)) + 1;
  Point p(n);
  while (true)
  {
    SchedStep S;
    S.a = a_vals[it];
    NT next_a = get_next_gaussian<PointGen>(P, p, a_vals[it], par.N, par.ratio, par.C,
                                            walk_length, rng);
    S.next_a = next_a;
    NT curr_fn = 0, curr_its = 0;
    WalkT walk(P, p, a_vals[it], rng);
    for (unsigned j = 0; j < totalSteps; j++)
    {
      walk.apply(P, p, a_vals[it], walk_length, rng);
      curr_its += 1.0;
      curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
      S.maxR2walk = std::max(S.maxR2walk, p.squared_length());
    }
    S.stopRatio = curr_fn / curr_its;
    if (next_a > 0 && curr_fn / curr_its > (1.0 + tol))
    {
      out.push_back(S);
      a_vals.push_back(next_a);
      it++;
    }
    else if (next_a <= 0)
    {
      S.branch = 1;
      out.push_back(S);
      a_vals.push_back(0.0);
      it++;
      break;
    }
    else
    {
      S.branch = 2;
      out.push_back(S);
      a_vals[it] = 0.0;
      break;
    }
  }
  return out;
}

// Diagnostic replica of get_next_gaussian's k-doubling loop (:150-176) for the
// ORIGINAL code: returns (k, sample relvar, sample mean) per tried k.
template <typename RNG>
std::vector<std::array<double, 3>> kTrace(HPoly& P, Point p, double a, RNG& rng,
                                          unsigned walk_length)
{
  gaussian_annealing_parameters<double> par(P.dimension());
  std::list<Point> pts;
  PushBackWalkPolicy pol;
  PointGen::apply(P, p, a, par.N, walk_length, pts, pol, rng);
  std::vector<double> fn(par.N);
  std::vector<std::array<double, 3>> out;
  double last_ratio = 0.1, k = 1;
  for (;;)
  {
    const double na = a * std::pow(par.ratio, k);
    auto fit = fn.begin();
    for (auto& q : pts) *fit++ = eval_exp(q, na) / eval_exp(q, a);
    auto mv = get_mean_variance(fn);
    out.push_back({k, mv.second / (mv.first * mv.first), mv.first});
    if (mv.second / (mv.first * mv.first) >= par.C || mv.first / last_ratio < 1.0 + 0.00001) break;
    k *= 2;
    last_ratio = mv.first;
    if (k > 1e6) break;
  }
  return out;
}

// ------------------------------------------------------- exact references
// logZ(a) = log int_P exp(-a|x|^2) for the SHIFTED, RESCALED polytope the oracle
// anneals over (center at 0, unit inner radius); vol is in original units.
struct Exact
{
  std::function<double(double)> logZ;
  double vol = 0;
  bool analytic = true;  // false: Monte Carlo (adds ~0.5% noise)
  // Gaussian mass of N(0, I/(2a)) outside P.
  std::function<double(double)> massOutside;
};

Exact exactBox(const Trace& T)
{
  std::vector<double> lo, hi;
  REQUIRE(boxBounds(T, lo, hi));
  Exact E;
  E.logZ = [lo, hi](double a) { return logZbox(lo, hi, a); };
  E.vol = std::exp(logZbox(lo, hi, 0.0)) * std::pow(T.scale, double(T.A.cols()));
  E.massOutside = [lo, hi](double a) { return massOutsideBox(lo, hi, a); };
  return E;
}

// Standard simplex {x >= 0, sum x <= s} shifted by -center: Z(a) by Monte Carlo
// with two estimators (uniform-on-P importance weights, and Gaussian hit rate),
// whichever has the smaller standard error.  Common random numbers across a.
struct SimplexMC
{
  int n;
  double vol;       // original units
  double logScale;  // oracle works on P / scale
  MT U;  // uniform samples on the shifted, rescaled P, n x M
  MT G;  // standard normals, n x M
  MT A;
  VT b;
  SimplexMC(int n_, double s, const Trace& T, int M, unsigned seed)
      : n(n_), logScale(std::log(T.scale)), A(T.A), b(T.b)
  {
    vol = std::pow(s, n) / std::tgamma(n + 1.0);
    std::mt19937_64 g(seed);
    std::exponential_distribution<double> Ex(1.0);
    std::normal_distribution<double> Nd(0, 1);
    U.resize(n, M);
    G.resize(n, M);
    for (int j = 0; j < M; ++j)
    {
      double tot = 0;
      std::vector<double> e(n + 1);
      for (auto& x : e) tot += (x = Ex(g));
      for (int i = 0; i < n; ++i) U(i, j) = (s * e[i] / tot - T.center(i)) / T.scale;
      for (int i = 0; i < n; ++i) G(i, j) = Nd(g);
    }
  }
  // returns (logZ, relative standard error)
  std::pair<double, double> uniformEst(double a) const
  {
    const int M = U.cols();
    VT l = -a * U.colwise().squaredNorm().transpose();
    const double mx = l.maxCoeff();
    VT w = (l.array() - mx).exp();
    const double mean = w.mean();
    const double var = (w.array() - mean).square().mean();
    return {std::log(vol) - n * logScale + mx + std::log(mean), std::sqrt(var / M) / mean};
  }
  std::pair<double, double> gaussEst(double a) const
  {
    const int M = G.cols();
    MT Y = G / std::sqrt(2 * a);
    MT S = A * Y;
    int hit = 0;
    for (int j = 0; j < M; ++j)
      if (((S.col(j) - b).array() <= 0).all()) ++hit;
    if (hit == 0) return {0, INFINITY};
    const double pr = double(hit) / M;
    return {0.5 * n * std::log(M_PI / a) + std::log(pr), std::sqrt((1 - pr) / (pr * M))};
  }
  double logZ(double a) const
  {
    auto u = uniformEst(a);
    if (a <= 0) return u.first;
    auto g = gaussEst(a);
    return g.second < u.second ? g.first : u.first;
  }
  double massOutside(double a) const
  {
    const auto g = gaussEst(a);
    return 1 - std::exp(g.first - 0.5 * n * std::log(M_PI / a));
  }
};

Exact exactSimplex(int n, double s, const Trace& T)
{
  auto mc = std::make_shared<SimplexMC>(n, s, T, 400000, 12345);
  Exact E;
  E.logZ = [mc](double a) { return mc->logZ(a); };
  E.vol = mc->vol;
  E.analytic = false;
  E.massOutside = [mc](double a) { return mc->massOutside(a); };
  return E;
}

double relVar(const Exact& E, double a, double an)
{
  const double l = E.logZ(2 * an - a) + E.logZ(a) - 2 * E.logZ(an);
  return l > 700 ? INFINITY : std::exp(l) - 1;
}

double exactRatio(const Exact& E, double a, double an) { return std::exp(E.logZ(an) - E.logZ(a)); }

void printPhases(const std::string& name, const Trace& T, const Exact& E)
{
  std::printf("== %s: est=%.6g exact=%.6g rel=%+.4f  center=[", name.c_str(), T.vol, E.vol,
              T.vol / E.vol - 1);
  for (int j = 0; j < T.center.size(); ++j) std::printf("%s%.4g", j ? "," : "", T.center(j));
  std::printf("] r=%.4g phases=%zu  (pi/a0)^(n/2)/Z(a0)=%.6f\n", T.radius, T.ratio.size(),
              std::exp(T.logVol0 - E.logZ(T.a[0])));
  for (size_t i = 0; i < T.ratio.size(); ++i)
    std::printf("  phase %2zu: a %.4e -> %.4e  est R=%-9.6g exact R=%-9.6g est/exact=%.4f "
                "relvar=%.3g steps=%.0f\n",
                i, T.a[i], T.a[i + 1], T.ratio[i], exactRatio(E, T.a[i], T.a[i + 1]),
                T.ratio[i] / exactRatio(E, T.a[i], T.a[i + 1]), relVar(E, T.a[i], T.a[i + 1]),
                T.steps[i]);
}

// ---------------------------------------------------------- test bodies
struct Case
{
  std::string name;
  Body body;
  HPoly P;   // with inner ball set (as ttc does before calling the oracle)
  Trace T;   // traced ttc-config oracle run
  Exact E;
};

Case makeCase(const std::string& name, const Body& B)
{
  Case c{name, B, toHPoly(B), {}, {}};
  volestiBall(c.P);
  c.T = traceTtc(c.P);
  if (name.rfind("simplex", 0) == 0) c.E = exactSimplex(B.A.cols(), 1.0, c.T);
  else c.E = exactBox(c.T);
  return c;
}

const std::vector<Case>& cases()
{
  static std::vector<Case> C = [] {
    std::vector<Case> v;
    v.push_back(makeCase("longbox[0,2]x[0,1000]", box({0, 0}, {2, 1000})));
    v.push_back(makeCase("unitcube2", cube(2, 0, 1)));
    v.push_back(makeCase("unitcube5", cube(5, 0, 1)));
    v.push_back(makeCase("unitcube10", cube(10, 0, 1)));
    v.push_back(makeCase("simplex5", simplex(5)));
    for (auto& c : v) printPhases(c.name, c.T, c.E);
    return v;
  }();
  return C;
}

constexpr double kC = 2.0;     // gaussian_annealing_parameters::C
constexpr double kFrac = 0.1;  // gaussian_annealing_parameters::frac

}  // namespace

TEST_SUITE("annealing")
{
TEST_CASE("replica reproduces volume_cooling_gaussians bit for bit")
{
  for (const auto& c : cases())
  {
    CAPTURE(c.name);
    CHECK(c.T.vol == ttcOracle(c.P));
    TtcVolumeRNG rng(c.P.dimension());
    std::vector<double> av;
    traceSchedule(c.P, rng, kTtcError, ttcWalkLength(c.P.dimension()), av);
    CHECK(av == c.T.a);
  }
}

TEST_CASE("schedule shape: a_0 > a_1 > ... > a_k = 0")
{
  for (const auto& c : cases())
  {
    CAPTURE(c.name);
    const auto& a = c.T.a;
    REQUIRE(a.size() >= 2);
    CHECK(a.back() == 0.0);
    for (size_t i = 0; i + 1 < a.size(); ++i)
    {
      CAPTURE(i);
      CHECK(a[i] > a[i + 1]);
      if (i + 2 < a.size()) CHECK(a[i + 1] > 0.0);
    }
  }
}

// get_first_gaussian (:59) uses the per-facet Mills-ratio tail bound, so the
// Gaussian N(0, I/(2 a_0)) loses at most frac*error of its mass outside P and
// (pi/a_0)^{n/2} overestimates Z(a_0) by at most a factor 1/(1 - frac*error).
TEST_CASE("a_0: Gaussian mass outside P <= frac * error")
{
  for (const auto& c : cases())
  {
    CAPTURE(c.name);
    const double a0 = c.T.a[0];
    // also via a direct call of the helper on the shifted polytope
    HPoly S(c.P);
    auto ib = S.ComputeInnerBall();
    S.shift(ib.first.getCoefficients());
    S.normalize();
    S.set_vec(S.get_vec() / ib.second);  // unit inner radius, as the oracle does
    std::vector<double> av;
    get_first_gaussian(S, kFrac, 1.0, kTtcError, av);
    REQUIRE(av.size() == 1);
    CHECK(av[0] == doctest::Approx(a0).epsilon(1e-12));
    const double out = c.E.massOutside(a0);
    const double slack = c.E.analytic ? 0.0 : 4 * std::sqrt(kFrac * kTtcError / 400000.0);
    std::printf("%s: a0=%.5g  mass outside=%.3e  bound frac*err=%.3e  (pi/a0)^(n/2)/Z(a0)=%.6f\n",
                c.name.c_str(), a0, out, kFrac * kTtcError,
                std::exp(c.T.logVol0 - c.E.logZ(a0)));
    CHECK(out <= kFrac * kTtcError + slack);
    CHECK(std::exp(c.T.logVol0 - c.E.logZ(a0)) - 1 <= kFrac * kTtcError / (1 - kFrac * kTtcError) + slack);
  }
}

// Cousins-Vempala: phase i is a Monte Carlo mean of Y = f_{i+1}/f_i under mu_i;
// its relative variance Var(Y)/E[Y]^2 = Z(2a_{i+1}-a_i) Z(a_i) / Z(a_{i+1})^2 - 1
// must stay O(1) (the code targets C = 2 on a sample estimate; we allow 2C).
TEST_CASE("every annealing step satisfies the variance condition (exact relvar <= 2C)")
{
  for (const auto& c : cases())
  {
    for (size_t i = 0; i + 1 < c.T.a.size(); ++i)
    {
      CAPTURE(c.name);
      CAPTURE(i);
      CAPTURE(c.T.a[i]);
      CAPTURE(c.T.a[i + 1]);
      CHECK(relVar(c.E, c.T.a[i], c.T.a[i + 1]) <= 2 * kC);
    }
  }
}

// Each phase runs until its own convergence test with curr_eps = error/sqrt(mm)
// (:383); allow 3x that (plus MC noise of the reference for the simplex).
TEST_CASE("per-phase ratio estimates match exact Z(a_{i+1})/Z(a_i)")
{
  for (const auto& c : cases())
  {
    const double tol = 3 * kTtcError / std::sqrt(double(c.T.ratio.size()))
                       + (c.E.analytic ? 0.0 : 0.01);
    for (size_t i = 0; i < c.T.ratio.size(); ++i)
    {
      CAPTURE(c.name);
      CAPTURE(i);
      CAPTURE(c.T.a[i]);
      CAPTURE(c.T.a[i + 1]);
      const double ex = exactRatio(c.E, c.T.a[i], c.T.a[i + 1]);
      CAPTURE(c.T.ratio[i]);
      CAPTURE(ex);
      CHECK(std::abs(c.T.ratio[i] / ex - 1) <= tol);
    }
  }
}

// The last phase estimates Z(0)/Z(a_last) = vol(P)/Z(a_last) from samples of
// mu_last; it is only correct if those samples "see" the whole body.
TEST_CASE("final phase a_last -> 0 recovers vol(P)/Z(a_last)")
{
  for (const auto& c : cases())
  {
    CAPTURE(c.name);
    const size_t L = c.T.ratio.size() - 1;
    const double ex = exactRatio(c.E, c.T.a[L], 0.0);
    std::printf("%s: last phase a=%.4e  est=%.6g exact=%.6g est/exact=%.4f\n", c.name.c_str(),
                c.T.a[L], c.T.ratio[L], ex, c.T.ratio[L] / ex);
    CHECK(c.T.a[L + 1] == 0.0);
    CHECK(std::abs(c.T.ratio[L] / ex - 1) <= 3 * kTtcError / std::sqrt(double(L + 1))
                                                 + (c.E.analytic ? 0.0 : 0.01));
  }
}

// Decomposition of the total log-error into the a_0 truncation and the phases.
TEST_CASE("telescoping product: total error and its decomposition")
{
  for (const auto& c : cases())
  {
    CAPTURE(c.name);
    double worst = 0;
    size_t worstI = 0;
    double sumLog = c.T.logVol0 - c.E.logZ(c.T.a[0]);
    for (size_t i = 0; i < c.T.ratio.size(); ++i)
    {
      const double e = std::log(c.T.ratio[i] / exactRatio(c.E, c.T.a[i], c.T.a[i + 1]));
      sumLog += e;
      if (std::abs(e) > std::abs(worst)) { worst = e; worstI = i; }
    }
    std::printf("%s: est/exact=%.5f  (a0 term %.5f, worst phase %zu factor %.5f, product of "
                "terms %.5f)\n",
                c.name.c_str(), c.T.vol / c.E.vol, std::exp(c.T.logVol0 - c.E.logZ(c.T.a[0])),
                worstI, std::exp(worst), std::exp(sumLog));
    // the per-term decomposition must account for the whole estimate
    CHECK(std::exp(sumLog) * c.E.vol == doctest::Approx(c.T.vol).epsilon(c.E.analytic ? 1e-6 : 0.03));
    CHECK(std::abs(c.T.vol / c.body.exactVolume - 1) <= 3 * kTtcError);
  }
}

// get_next_gaussian (:128) doubles k while the SAMPLE relvar of
// exp((a-a')|x|^2) stays < C.  Samples near the start (|x|^2 <= ~7 on the long
// box) bound that statistic no matter how far the body extends, so k can run
// away (k=32 at a=0.517: a' = 1.2e-10).  Call the helper directly, starting at
// the Chebyshev center as the oracle does, and check the step it proposes.
TEST_CASE("get_next_gaussian proposes a step that satisfies the exact variance condition")
{
  struct Probe
  {
    std::string name;
    Body B;
    std::vector<double> as;
  };
  std::vector<Probe> probes = {
      {"longbox", box({0, 0}, {2, 1000}), {4.1336, 0.5167, 0.05, 5e-3, 5e-4}},
      {"unitcube10", cube(10, 0, 1), {23.4, 2.85}},
      {"unitcube2", cube(2, 0, 1), {17.5, 8.78}},
  };
  for (auto& pr : probes)
  {
    HPoly P = toHPoly(pr.B);
    volestiBall(P);
    HPoly S(P);
    auto ib = S.ComputeInnerBall();
    S.shift(ib.first.getCoefficients());
    Trace geom;
    geom.A = S.get_mat();
    geom.b = S.get_vec();
    const Exact E = exactBox(geom);
    const unsigned n = S.dimension();
    gaussian_annealing_parameters<double> par(n);
    for (double a : pr.as)
    {
      TtcVolumeRNG rng(n);
      Point p(n);
      const double na = get_next_gaussian<PointGen>(S, p, a, par.N, par.ratio, par.C,
                                                    ttcWalkLength(n), rng);
      TtcVolumeRNG rng2(n);
      auto kt = kTrace(S, Point(n), a, rng2, ttcWalkLength(n));
      std::printf("%s: a=%.4g -> a'=%.4g (a'/a=%.3g)  exact relvar=%.3g  k-trace(orig):",
                  pr.name.c_str(), a, na, na / a, relVar(E, a, na));
      for (auto& t : kt) std::printf(" k=%g:%.3g", t[0], t[1]);
      std::printf("\n");
      CAPTURE(pr.name);
      CAPTURE(a);
      CAPTURE(na);
      CHECK(na > 0);
      CHECK(na < a);
      CHECK(relVar(E, a, na) <= 2 * kC);
    }
  }
}

// The stopping branch `else { a_vals[it] = a_stop; break; }` (:260) REPLACES
// the current a_it by 0 instead of appending 0, so the last phase becomes
// a_{it-1} -> 0.  This is harmless only if the dropped a_it is already
// effectively uniform on P: Z(0)/Z(a_it) - 1 small compared with the per-phase
// accuracy.
TEST_CASE("stopping test: a dropped (replaced) a_it is effectively uniform on P")
{
  for (const auto& c : cases())
  {
    TtcVolumeRNG rng(c.P.dimension());
    std::vector<double> av;
    auto steps = traceSchedule(c.P, rng, kTtcError, ttcWalkLength(c.P.dimension()), av);
    const auto& last = steps.back();
    std::printf("%s: stop at a=%.4e next=%.4e stopRatio=%.6f branch=%d max|x|^2(walk)=%.3g "
                "Z(0)/Z(a)=%.6g\n",
                c.name.c_str(), last.a, last.next_a, last.stopRatio, last.branch,
                last.maxR2walk, exactRatio(c.E, last.a, 0.0));
    CAPTURE(c.name);
    if (last.branch == 2)
      CHECK(exactRatio(c.E, last.a, 0.0) - 1
            <= kTtcError / std::sqrt(double(av.size() - 1)) + (c.E.analytic ? 0.0 : 0.01));
  }
}
}
