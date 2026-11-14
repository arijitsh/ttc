// Audit of the per-phase ratio estimation and its stopping rule in
// volume_cooling_gaussians (src/volesti/include/volume/volume_cooling_gaussians.hpp,
// the `for (viterator fnIt ...)` loop, lines ~376-443).
//
// Part (a): the window/min/max/stop logic, replicated verbatim in isolation and
//           driven by synthetic streams with a known mean.
// Part (b): the full oracle (ttc configuration: GaussianBallWalk, walk length
//           10+n^2) over many seeds, with an instrumented but RNG-identical
//           replica that records each phase; for the cube [-1,1]^n the Gaussian
//           integrals are known in closed form, so the total error can be split
//           exactly into "first Gaussian" error and per-ratio errors.
//
// Tests tagged "slow" (test-suite) are skipped by default; run them with
//   bin/test_ratio_convergence --no-skip -ts=slow
#include "doctest.h"
#include "common.hpp"

#include <cfloat>
#include <cstdio>
#include <iostream>
#include <limits>

using namespace vt;

namespace rc
{

// Relative standard error of the mean of a (possibly autocorrelated) stream,
// non-overlapping batch means with ~sqrt(N) batches (same as the proposed fix).
inline double batchMeansRelSE(const std::vector<double>& y, double mean)
{
  const size_t N = y.size();
  const size_t nb = std::max<size_t>(2, static_cast<size_t>(std::sqrt(double(N))));
  const size_t bs = N / nb;
  if (bs == 0 || mean <= 0) return std::numeric_limits<double>::max();
  std::vector<double> bm(nb, 0.0);
  double m = 0;
  for (size_t b = 0; b < nb; ++b)
  {
    for (size_t i = 0; i < bs; ++i) bm[b] += y[b * bs + i];
    bm[b] /= double(bs);
    m += bm[b] / nb;
  }
  double v = 0;
  for (double x : bm) v += (x - m) * (x - m);
  v /= double(nb - 1);
  return std::sqrt(v / nb) / mean;
}

// ===================================================================== (a)
// Verbatim replica of the convergence test inside the ratio loop.  Only the
// source of `val` differs: here it is the running mean of a synthetic stream.
struct StopRule
{
  double curr_eps;
  unsigned W;
  std::vector<double> last_W;
  double min_val = std::numeric_limits<double>::min();  // as in volesti: DBL_MIN (!)
  double max_val = std::numeric_limits<double>::max();
  unsigned min_index, max_index, index = 0;
  double fn = 0, its = 0;
  bool done = false;
  std::vector<unsigned> slotWrites;  // audit: how often each slot is written
  // Proposed fix (off by default): once the window test passes, also require
  // the batch-means relative standard error of the mean to be <= curr_eps/2.
  bool batchMeans = false;
  std::vector<double> ys;
  double nextCheck = 0;

  StopRule(double ce, unsigned W_, bool bm = false)
      : curr_eps(ce), W(W_), last_W(W_, 0.0), min_index(W_ - 1), max_index(W_ - 1),
        slotWrites(W_, 0), batchMeans(bm)
  {
  }

  double mean() const { return fn / its; }

  // One iteration of `while (!done || its < min_steps)` with min_steps == 0.
  void push(double x)
  {
    its += 1.0;
    fn += x;
    if (batchMeans) ys.push_back(x);
    double val = fn / its;
    last_W[index] = val;
    ++slotWrites[index];
    if (val <= min_val) { min_val = val; min_index = index; }
    else if (min_index == index)
    {
      auto it = std::min_element(last_W.begin(), last_W.end());
      min_val = *it;
      min_index = static_cast<unsigned>(std::distance(last_W.begin(), it));
    }
    if (val >= max_val) { max_val = val; max_index = index; }
    else if (max_index == index)
    {
      auto it = std::max_element(last_W.begin(), last_W.end());
      max_val = *it;
      max_index = static_cast<unsigned>(std::distance(last_W.begin(), it));
    }
    if ((max_val - min_val) / max_val <= curr_eps / 2.0)
    {
      if (!batchMeans) done = true;
      else if (its >= nextCheck)
      {
        if (batchMeansRelSE(ys, mean()) <= curr_eps / 2.0) done = true;
        else nextCheck = 1.2 * its;
      }
    }
    index = index % W + 1;
    if (index == W) index = 0;
  }
};

// Run the rule on a stream until it stops (or maxN).
template <class Gen>
StopRule runRule(Gen&& next, double ce, unsigned W, bool bm = false, size_t maxN = 50000000)
{
  StopRule R(ce, W, bm);
  while (!R.done && R.its < maxN) R.push(next());
  return R;
}

struct StreamStats
{
  double meanN = 0, bias = 0, rms = 0, sd = 0;
};

// Log-normal stream with mean 1 and relative sd sigma: y = exp(s z - s^2/2),
// sigma^2 = e^{s^2}-1, where z is a stationary AR(1) N(0,1) process with lag-1
// correlation rho (rho = 0: i.i.d.).  rho > 0 mimics an MCMC ratio stream with
// autocorrelation (tau_int ~ (1+rho)/(1-rho) for small s).
inline StreamStats streamExperiment(double sigma, double rho, double ce, unsigned W, int reps,
                                    unsigned seed, bool bm = false)
{
  std::mt19937_64 g(seed);
  std::normal_distribution<double> Z(0.0, 1.0);
  const double s = std::sqrt(std::log1p(sigma * sigma));
  const double q = std::sqrt(1.0 - rho * rho);
  StreamStats S;
  std::vector<double> rel;
  for (int r = 0; r < reps; ++r)
  {
    double z = Z(g);
    auto next = [&]() {
      z = rho * z + q * Z(g);
      return std::exp(s * z - 0.5 * s * s);
    };
    StopRule R = runRule(next, ce, W, bm);
    S.meanN += R.its / reps;
    rel.push_back(R.mean() - 1.0);
  }
  for (double e : rel) { S.bias += e / reps; S.rms += e * e / reps; }
  S.rms = std::sqrt(S.rms);
  S.sd = std::sqrt(std::max(0.0, S.rms * S.rms - S.bias * S.bias));
  return S;
}

// ===================================================================== (b)
struct PhaseTrace
{
  double a = 0, aNext = 0;
  double N = 0;          // samples used for this ratio
  double ratio = 0;      // estimate fn/its
  double relVar = 0;     // sample Var(Y)/E(Y)^2 of the ratio stream
  double tauInt = 1;     // integrated autocorrelation time (batch means)
};

struct Trace
{
  std::vector<double> a_vals;
  std::vector<PhaseTrace> phases;
  double firstGaussVol = 0;
  double vol = 0;
  double radius = 0;
};

// Batch-means estimate of the integrated autocorrelation time.
inline double tauBatch(const std::vector<double>& y)
{
  const size_t n = y.size();
  if (n < 400) return 1.0;
  const size_t nb = static_cast<size_t>(std::sqrt(double(n)));
  const size_t bs = n / nb;
  double m = 0, v = 0;
  for (size_t i = 0; i < nb * bs; ++i) m += y[i];
  m /= double(nb * bs);
  for (size_t i = 0; i < nb * bs; ++i) v += (y[i] - m) * (y[i] - m);
  v /= double(nb * bs - 1);
  double vb = 0;
  for (size_t b = 0; b < nb; ++b)
  {
    double s = 0;
    for (size_t i = 0; i < bs; ++i) s += y[b * bs + i];
    s /= double(bs);
    vb += (s - m) * (s - m);
  }
  vb /= double(nb - 1);
  return v > 0 ? std::max(1.0, bs * vb / v) : 1.0;
}

// Exact i.i.d. sample from exp(-a|x|^2) restricted to the cube [-1,1]^n
// (product of 1-D truncated Gaussians, by rejection).
inline void exactCubeGauss(std::mt19937_64& g, int n, double a, Point& p)
{
  std::uniform_real_distribution<double> U(-1.0, 1.0), U01(0.0, 1.0);
  std::normal_distribution<double> Z(0.0, 1.0);
  const double sd = a > 0 ? std::sqrt(0.5 / a) : 1e300;
  for (int i = 0; i < n; ++i)
  {
    double x;
    if (sd < 0.5)
      do x = sd * Z(g); while (std::abs(x) > 1.0);
    else
      do x = U(g); while (U01(g) > std::exp(-a * x * x));
    p.set_coord(i, x);
  }
}

// RNG-identical, instrumented copy of volume_cooling_gaussians<GaussianBallWalk>
// (same calls, same order, same RNG consumption).  `stopScale` multiplies the
// number of samples demanded after the original rule fires (1 = original).
// If iidSeed != 0 (cube [-1,1]^n only), the ratio phase draws exact i.i.d.
// samples instead of running the ball walk: isolates the stopping rule /
// estimator from MCMC error, on the same annealing schedule.
template <class RNG>
double tracedVolume(HPoly Pin, RNG& rng, double error, unsigned walk_length, Trace& T,
                    double stopScale = 1.0, unsigned iidSeed = 0)
{
  std::mt19937_64 iid(iidSeed);
  using WalkType = GaussianBallWalk::Walk<HPoly, RNG>;
  using RPG = GaussianRandomPointGenerator<WalkType>;
  auto P(Pin);
  const unsigned n = P.dimension();
  gaussian_annealing_parameters<NT> par(n);
  auto IB = P.ComputeInnerBall();
  if (IB.second < 0.0) return -1.0;
  Point c = IB.first;
  NT radius = IB.second;
  T.radius = radius;
  P.shift(c.getCoefficients());
  std::vector<NT> a_vals;
  compute_annealing_schedule<WalkType, RPG>(P, par.ratio, par.C, par.frac, par.N, walk_length,
                                            radius, error, a_vals, rng);
  T.a_vals = a_vals;
  const unsigned W = par.W;
  const unsigned mm = a_vals.size() - 1;
  NT vol = std::pow(M_PI / a_vals[0], NT(n) / 2.0);
  T.firstGaussVol = vol;
  Point p(n);
  for (unsigned i = 0; i < mm; ++i)
  {
    const double a = a_vals[i], an = a_vals[i + 1];
    StopRule R(error / std::sqrt(NT(mm)), W);
    WalkType walk(P, p, a, rng);
    std::vector<double> ys;
    double extra = 0;  // samples still required by stopScale
    bool fired = false;
    while (true)
    {
      if (iidSeed) exactCubeGauss(iid, n, a, p);
      else walk.apply(P, p, a, walk_length, rng);
      double y = eval_exp(p, an) / eval_exp(p, a);
      ys.push_back(y);
      R.push(y);
      if (R.done && !fired)
      {
        fired = true;
        extra = (stopScale - 1.0) * R.its;
      }
      if (fired)
      {
        if (extra <= 0) break;
        extra -= 1.0;
      }
    }
    PhaseTrace ph;
    ph.a = a;
    ph.aNext = an;
    ph.N = R.its;
    ph.ratio = R.mean();
    auto mv = get_mean_variance(ys);
    ph.relVar = mv.second / (mv.first * mv.first);
    ph.tauInt = tauBatch(ys);
    T.phases.push_back(ph);
    vol *= R.mean();
  }
  T.vol = vol;
  return vol;
}

// Closed-form Gaussian integral over the cube [-1,1]^n: (sqrt(pi/a) erf(sqrt a))^n.
inline double cubeGaussLog(int n, double a)
{
  if (a <= 0) return n * std::log(2.0);
  return n * std::log(std::sqrt(M_PI / a) * std::erf(std::sqrt(a)));
}

inline double runSeed(const Body& B, unsigned seed, double error)
{
  HPoly P = toHPoly(B);
  volestiBall(P);
  return gaussianOracle(P, seed, error, ttcWalkLength(static_cast<int>(B.A.cols())));
}

struct SeedSweep
{
  Stats st;
  double meanPhases = 0, meanSamples = 0;
  // cube-only decomposition (log-scale), averaged over seeds
  double firstErr = 0, ratioErrRms = 0, ratioErrBudget = 0, ratioErrMean = 0, sumRatioLogErr = 0;
};

inline SeedSweep sweep(const Body& B, double error, int seeds, unsigned seed0, bool cubeDecomp,
                       double stopScale = 1.0, bool iid = false)
{
  const int n = static_cast<int>(B.A.cols());
  std::vector<double> est;
  SeedSweep S;
  double rr = 0;
  int cnt = 0;
  for (int s = 0; s < seeds; ++s)
  {
    HPoly P = toHPoly(B);
    volestiBall(P);
    TtcSampleRNG rng(n);
    rng.set_seed(seed0 + s);
    Trace T;
    est.push_back(tracedVolume(P, rng, error, ttcWalkLength(n), T, stopScale, iid ? 7919u * (seed0 + s) + 1 : 0u));
    S.meanPhases += double(T.phases.size()) / seeds;
    double tot = 0;
    for (auto& ph : T.phases) tot += ph.N;
    S.meanSamples += tot / seeds;
    if (cubeDecomp)
    {
      S.firstErr += (std::log(T.firstGaussVol) - cubeGaussLog(n, T.a_vals[0])) / seeds;
      for (auto& ph : T.phases)
      {
        double e = std::log(ph.ratio) - (cubeGaussLog(n, ph.aNext) - cubeGaussLog(n, ph.a));
        rr += e * e;
        S.ratioErrMean += e;
        S.sumRatioLogErr += e / seeds;
        ++cnt;
      }
      S.ratioErrBudget = error / std::sqrt(double(T.phases.size()));
    }
  }
  if (cnt) { S.ratioErrRms = std::sqrt(rr / cnt); S.ratioErrMean /= cnt; }
  S.st = relStats(est, B.exactVolume, error);
  return S;
}

inline void printSweep(const char* tag, const Body& B, double error, const SeedSweep& S)
{
  std::printf("  %-8s %-14s err=%.4f phases=%.1f samples=%.0f | %s", tag, B.name.substr(0, 14).c_str(),
              error, S.meanPhases, S.meanSamples, fmt(S.st).c_str());
  if (S.ratioErrRms > 0)
    std::printf(" | firstGauss=%+.4f ratio: mean=%+.4f rms=%.4f (budget err/sqrt(mm)=%.4f) "
                "sum-of-ratio-bias=%+.4f",
                S.firstErr, S.ratioErrMean, S.ratioErrRms, S.ratioErrBudget, S.sumRatioLogErr);
  std::printf("\n");
}

}  // namespace rc

using namespace rc;

// ------------------------------------------------------------------ part (a)
TEST_CASE("stop rule: index wraparound writes every slot 0..W-1 exactly once per W pushes")
{
  for (unsigned W : {1u, 2u, 7u, 950u})
  {
    StopRule R(1e-300, W);  // never stops
    for (unsigned k = 0; k < 5 * W; ++k) R.push(1.0 + 1e-3 * (k % 3));
    for (unsigned s = 0; s < W; ++s) CHECK_MESSAGE(R.slotWrites[s] == 5, "W=" << W << " slot " << s);
  }
}

TEST_CASE("stop rule: DBL_MIN/DBL_MAX sentinels => never stops before the window is full")
{
  // numeric_limits<double>::min() is the smallest *positive* double, not the
  // most negative one, so `val <= min_val` is never true for a ratio >= 1 and
  // min/max are first computed when the sentinel index W-1 is overwritten,
  // i.e. after exactly W samples.  The effect is a hidden min_steps == W.
  const unsigned W = 950;
  StopRule C = runRule([] { return 1.25; }, 0.02, W);  // constant stream
  CHECK(C.its == W);
  CHECK(C.min_val == doctest::Approx(1.25));

  std::mt19937 g(1);
  std::uniform_real_distribution<double> U(1.0, 1.0001);  // tiny spread: would stop at once
  for (int r = 0; r < 20; ++r)
  {
    StopRule R = runRule([&] { return U(g); }, 0.02, W);
    CHECK(R.its == W);
  }
}

TEST_CASE("stop rule: tracked min/max equal brute-force extrema of the last W running means")
{
  const unsigned W = 97;
  std::mt19937 g(7);
  std::exponential_distribution<double> E(1.0);
  StopRule R(1e-300, W);
  std::vector<double> hist;
  for (int k = 0; k < 20 * int(W); ++k)
  {
    R.push(E(g));
    hist.push_back(R.mean());
    if (hist.size() >= W)
    {
      double mn = *std::min_element(hist.end() - W, hist.end());
      double mx = *std::max_element(hist.end() - W, hist.end());
      REQUIRE(R.min_val == mn);
      REQUIRE(R.max_val == mx);
    }
  }
}

TEST_CASE("stop rule on synthetic streams: window rule misses the curr_eps budget; "
          "window + batch-means meets it")
{
  // Each ratio is allotted curr_eps = error/sqrt(mm); for the product of mm
  // ratios to have sd <= error/2 (=> ~95% within `error`) each ratio needs
  // relative sd <= curr_eps/2.  The window rule stops when the running mean
  // moved by <= curr_eps/2 over the last W steps, i.e. N ~ sigma*sqrt(tau*W)/curr_eps,
  // so its delivered sd ~ sqrt(curr_eps*sigma*sqrt(tau/W)): it scales like the
  // *square root* of the target and ignores autocorrelation (tau) almost
  // entirely.  This test characterises the replica of the current rule (it
  // asserts the flaw) and validates the proposed batch-means rule.
  const unsigned W = 6 * 25 + 800;  // n = 5
  const int reps = 100;
  std::printf("[synthetic] lognormal AR(1) streams, mean 1, W=%u, reps=%d\n", W, reps);
  std::printf("  %5s %5s %7s | %9s %9s %7s | %9s %9s %7s\n", "sigma", "rho", "ce", "N(win)",
              "sd(win)", "sd/ce", "N(win+bm)", "sd(w+bm)", "sd/ce");
  int winOk = 0, bmOk = 0, total = 0, bmTotal = 0;
  double winSdOverCe[2][2][3] = {};
  int ri = 0;
  for (double rho : {0.0, 0.8})  // 0.8: tau_int ~ 9
  {
    int si = 0;
    for (double sigma : {0.5, 1.4})  // 1.4: Var/E^2 ~ 2 = volesti's C, as seen in real phases
    {
      int ci = 0;
      for (double ce : {0.06, 0.02, 0.01})
      {
        const unsigned seed = 1234u + unsigned(1000 * ce) + 17u * ri + 5u * si;
        StreamStats A = streamExperiment(sigma, rho, ce, W, reps, seed, false);
        winSdOverCe[ri][si][ci] = A.sd / ce;
        ++total;
        winOk += A.sd <= ce / 2;
        CHECK(std::abs(A.bias) <= 3 * A.sd / std::sqrt(double(reps)) + 1e-3);
        // batch-means needs ~4 sigma^2 tau / ce^2 samples; skip the costliest cells
        const bool runBm = rho == 0.0 || ce >= 0.02;
        StreamStats B;
        if (runBm)
        {
          B = streamExperiment(sigma, rho, ce, W, reps, seed, true);
          ++bmTotal;
          bmOk += B.sd <= 0.6 * ce;  // ce/2 target + slack for 100 reps & batch-means noise
          CHECK_MESSAGE(B.sd <= 0.6 * ce,
                        "batch-means rule: sigma=" << sigma << " rho=" << rho << " ce=" << ce);
        }
        std::printf("  %5.2f %5.2f %7.4f | %9.0f %9.5f %7.3f | %9.0f %9.5f %7.3f\n", sigma, rho, ce,
                    A.meanN, A.sd, A.sd / ce, B.meanN, B.sd, B.sd / ce);
        ++ci;
      }
      ++si;
    }
    ++ri;
  }
  std::printf("  meeting budget: window rule %d/%d, window+batch-means %d/%d\n", winOk, total, bmOk,
              bmTotal);
  // Flaw of the current rule: relative to its budget the delivered error grows
  // as the budget shrinks (sd/ce at ce=0.01 well above sd/ce at ce=0.06), and
  // autocorrelation inflates it further.
  for (int r = 0; r < 2; ++r)
    for (int s2 = 0; s2 < 2; ++s2) CHECK(winSdOverCe[r][s2][2] > 2.0 * winSdOverCe[r][s2][0]);
  CHECK(winSdOverCe[1][1][1] > 1.0);  // sigma=1.4, tau~9, ce=0.02: sd > the whole budget
  CHECK(winOk < total / 2);
}

// ------------------------------------------------------------------ part (b)
TEST_CASE("instrumented replica is RNG-identical to volume_cooling_gaussians")
{
  for (const Body& B : {cube(3), simplex(4)})
    for (unsigned seed : {11u, 12u})
    {
      HPoly P1 = toHPoly(B), P2 = toHPoly(B);
      volestiBall(P1);
      volestiBall(P2);
      const int n = static_cast<int>(B.A.cols());
      double lib = gaussianOracle(P1, seed, 0.8 / 12, ttcWalkLength(n));
      TtcSampleRNG rng(n);
      rng.set_seed(seed);
      Trace T;
      double rep = tracedVolume(P2, rng, 0.8 / 12, ttcWalkLength(n), T);
      // The replica encodes the HEAD stopping rule; once the library's rule is
      // fixed this only warns (the replica-based diagnostics then describe the
      // old rule, which is what they are meant to audit).
      WARN_MESSAGE(rep == lib, "replica differs from library (rule changed?): " << rep << " vs " << lib);
    }
}

TEST_CASE("cube [-1,1]^3: exact per-phase decomposition of one ttc run")
{
  const int n = 3;
  Body B = cube(n);
  HPoly P = toHPoly(B);
  volestiBall(P);
  TtcSampleRNG rng(n);
  rng.set_seed(5);
  Trace T;
  const double err = 0.8 / 12;
  double v = tracedVolume(P, rng, err, ttcWalkLength(n), T);
  const double ce = err / std::sqrt(double(T.phases.size()));
  std::printf("[cube3 seed 5] vol=%.4f exact=%.0f rel=%+.4f mm=%zu curr_eps=%.4f\n", v,
              B.exactVolume, v / B.exactVolume - 1, T.phases.size(), ce);
  std::printf("  first Gaussian a0=%.3f log-err=%+.2e\n", T.a_vals[0],
              std::log(T.firstGaussVol) - cubeGaussLog(n, T.a_vals[0]));
  for (size_t i = 0; i < T.phases.size(); ++i)
  {
    const auto& ph = T.phases[i];
    double ex = std::exp(cubeGaussLog(n, ph.aNext) - cubeGaussLog(n, ph.a));
    std::printf("  phase %2zu a=%9.4f->%9.4f N=%6.0f ratio=%.5f exact=%.5f rel=%+.4f relVar=%.3f "
                "tau=%.1f\n",
                i, ph.a, ph.aNext, ph.N, ph.ratio, ex, ph.ratio / ex - 1, ph.relVar, ph.tauInt);
    // the stopping rule hides a min_steps == W
    CHECK(ph.N >= 6 * n * n + 800);
  }
  // first-Gaussian truncation must consume <= frac(=0.1) of the error budget
  CHECK(std::abs(std::log(T.firstGaussVol) - cubeGaussLog(n, T.a_vals[0])) <= 0.1 * err);
}

// Seed sweep through the *library* oracle (so it picks up a fix in src/volesti).
inline Stats libSweep(const Body& B, double error, int seeds, unsigned seed0)
{
  std::vector<double> est;
  for (int s = 0; s < seeds; ++s) est.push_back(runSeed(B, seed0 + s, error));
  return relStats(est, B.exactVolume, error);
}

// Correct behaviour of an (error)-volume oracle whose budget is split as
// sqrt(mm) * curr_eps = error: rel. sd <= error/2 (=> ~95% within `error`).
// With 20-40 seeds the sample sd has ~12-15% relative noise, so we assert
// sd <= 0.6*error, >= 80% of runs within `error`, and a bias bound of
// frac*error (the deterministic first-Gaussian truncation, frac = 0.1) plus
// 2.5 standard errors of the seed mean.
inline void checkOracle(const Body& B, double err, const Stats& st)
{
  std::printf("  %-14s err=%.4f %s sd/err=%.2f\n", B.name.substr(0, 14).c_str(), err, fmt(st).c_str(),
              st.sd / err);
  CHECK_MESSAGE(st.sd <= 0.6 * err, B.name << " err=" << err << " sd=" << st.sd);
  const double biasTol = 0.1 * err + 2.5 * st.sd / std::sqrt(double(st.n - st.nonFinite));
  CHECK_MESSAGE(std::abs(st.meanRel) <= biasTol,
                B.name << " err=" << err << " bias=" << st.meanRel << " tol=" << biasTol);
  CHECK_MESSAGE(st.fracWithin >= 0.8, B.name << " err=" << err << " within=" << st.fracWithin);
  CHECK(st.nonFinite == 0);
}

TEST_CASE("full oracle over seeds at ttc's error 0.8/12: simplex4 and cube3")
{
  const double err = 0.8 / 12;
  std::printf("[seeds] library oracle, ttc config, error=%.4f\n", err);
  for (const Body& B : {simplex(4), cube(3)})
    checkOracle(B, err, libSweep(B, err, 20, 100));
}

TEST_CASE("full oracle at error 0.03: requested accuracy is not delivered")
{
  // Measured on HEAD: cube3 sd~0.031 (sd/err 1.05, 65% within).
  const double err = 0.03;
  std::printf("[seeds] library oracle, error=%.4f\n", err);
  Body B = cube(3);
  checkOracle(B, err, libSweep(B, err, 20, 200));
}

TEST_CASE("cube5/cube8/simplex4 over seeds, error in {0.2,0.1,0.0667,0.03}" *
          doctest::test_suite("slow") * doctest::skip())
{
  for (const Body& B : {cube(5), simplex(4), cube(8)})
    for (double err : {0.2, 0.1, 0.8 / 12, 0.03})
      checkOracle(B, err, libSweep(B, err, 40, 1000));
}

TEST_CASE("diagnostic: replica sweep with exact cube decomposition" * doctest::test_suite("slow") * doctest::skip())
{
  // firstGauss = log((pi/a0)^{n/2} / I(a0)) is the deterministic truncation bias
  // (<= frac*error by construction); ratio mean/rms are log-errors of the
  // individual ratios against closed-form cube Gaussian integrals.  iid=1 runs
  // the same rule on exact i.i.d. samples (no MCMC error).
  for (int n : {3, 5, 8})
    for (double err : {0.2, 0.8 / 12, 0.03})
      for (bool iid : {false, true})
      {
        Body B = cube(n);
        SeedSweep S = sweep(B, err, 20, 700, true, 1.0, iid);
        printSweep(iid ? "iid" : "mcmc", B, err, S);
      }
}
