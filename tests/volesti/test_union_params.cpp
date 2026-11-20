// T4 of text-for-agents/ttc-union-test-plan.md: the parameters of Alg. 1
// (KR25 manuscript, sections/algo.tex) and the phase-1b volume clipping.
//
//   [thresh]   lraSampleThreshold == Alg. 1 line 5,
//              max(24 ln(24/delta) / ((1-e') e'^2), 15 (ln(6/delta) + ln m)), e' = eps/12,
//              recomputed independently in long double on a grid of (eps, delta, m).
//   [proof]    Thresh meets the two inequalities the proof of lem:numlattice needs
//              (4 exp(-(1-e')e'^2 T/24) <= delta/6, exp(-T/15) <= delta/(6m)), and the
//              constant inequalities of that proof and of thm:volumettc hold for eps in (0,1].
//   [wiring]   ./build/ttc reports e' = eps/12, delta' = delta/(4m), Thresh(eps, delta, m)
//              and b = latticePrecision(n, r, m, gamma_min, eps/96) with m = all cubes of the
//              DNF (including the zero-volume one it drops), r = facets of the kept ones.
//   [clip]     the clip of src/volume/vol_est.cpp (ballVolumeBounds),
//              volume <- clamp(volume, (1-e')(2 gamma/sqrt n)^n, (1+e')(2R)^n),
//              is sound: certifiedInnerBall's ball lies in P, gamma <= the true inradius,
//              outerRadius R >= max_v |v - centre| over the true vertices, and
//              (2 gamma/sqrt n)^n <= vol <= (2R)^n, for boxes, simplices, cross-polytopes,
//              rotated, skinny, translated, tiny and huge bodies in dims 1..10; a clamped
//              (1 +- e') oracle value stays in (1 +- e') vol.
//   [scale]    dims up to 40 and sides 1e-7..1e7: the clip bounds computed in double never
//              produce NaN / lower > upper / a wrong finite value while vol is a normal
//              double, and estimateUnionVolume's log2Scale = log2 T + b n log2 10 bookkeeping
//              returns vol(P) (Poisson z-test, |z| <= 5) for volumes 1e-300..1e300.
#include "doctest.h"
#include "common.hpp"

#include "volume/get_precision.hpp"
#include "volume/union_estimator.hpp"
#include "volume/union_parameters.hpp"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>

#include <unistd.h>

using namespace vt;

namespace
{

// ---------------------------------------------------------------- Alg. 1 l.5
long double paperThresh(long double eps, long double delta, long double m)
{
  const long double e = eps / 12.0L;
  const long double a = 24.0L * std::log(24.0L / delta) / ((1.0L - e) * e * e);
  const long double c = 15.0L * (std::log(6.0L / delta) + std::log(m));
  return std::max(a, c);
}

const std::vector<double> kEps = {0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 0.95, 1.0};
const std::vector<double> kDelta = {1e-9, 1e-4, 0.01, 0.05, 0.1, 0.2, 0.5, 0.9, 0.999};
const std::vector<double> kM = {1, 2, 3, 10, 100, 1e4, 1e6, 1e9, 1e12, 1e15};

// ------------------------------------------------------- clip of vol_est.cpp
// Mirror of src/volume/vol_est.cpp (VolEst clipping): returns false where
// ttc throws.
bool ttcClip(double volume, double eps, int n, double gamma, double R, double& out,
             double* lowerOut = nullptr, double* upperOut = nullptr)
{
  const double mvcEps = eps / 12.0;
  const double lower = (1.0 - mvcEps) * std::pow(2.0 * gamma / std::sqrt(double(n)), double(n));
  const double upper = (1.0 + mvcEps) * std::pow(2.0 * R, double(n));
  if (lowerOut) *lowerOut = lower;
  if (upperOut) *upperOut = upper;
  if (std::isnan(lower) || std::isnan(upper) || lower > upper) return false;
  out = std::clamp(volume, lower, upper);
  return std::isfinite(out) && out > 0.0;
}

// ------------------------------------------------------- bodies with vertices
struct VBody
{
  std::string name;
  MT A;
  VT b;
  double logVol;    // ln vol
  double inradius;  // exact; < 0 if unknown
  std::vector<VT> verts;
};

VBody vbox(const std::vector<double>& lo, const std::vector<double>& hi)
{
  const int n = int(lo.size());
  Body B = box(lo, hi);
  VBody V{"box", B.A, B.b, 0.0, 0.0, {}};
  double r = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i)
  {
    V.logVol += std::log(hi[i] - lo[i]);
    r = std::min(r, 0.5 * (hi[i] - lo[i]));
  }
  V.inradius = r;
  for (int k = 0; k < (1 << n); ++k)
  {
    VT v(n);
    for (int i = 0; i < n; ++i) v(i) = (k >> i) & 1 ? hi[i] : lo[i];
    V.verts.push_back(v);
  }
  return V;
}

VBody vsimplex(int n)
{
  Body B = simplex(n);
  VBody V{"simplex", B.A, B.b, -std::lgamma(n + 1.0), 1.0 / (n + std::sqrt(double(n))), {}};
  V.verts.push_back(VT::Zero(n));
  for (int i = 0; i < n; ++i) V.verts.push_back(VT::Unit(n, i));
  return V;
}

VBody vcross(int n)
{
  Body B = crossPolytope(n);
  VBody V{"cross", B.A, B.b, n * std::log(2.0) - std::lgamma(n + 1.0), 1.0 / std::sqrt(double(n)),
          {}};
  for (int i = 0; i < n; ++i)
  {
    V.verts.push_back(VT::Unit(n, i));
    V.verts.push_back(-VT::Unit(n, i));
  }
  return V;
}

// x -> s Q x + t, Q orthogonal: inradius scales by s.
VBody vmap(const VBody& V, const MT& Q, double s, const VT& t, const std::string& tag)
{
  const int n = int(V.A.cols());
  Body B{V.name, V.A, V.b, 1.0};
  MT M = s * Q;
  Body R = linearImage(B, M, t, tag);
  VBody W{V.name + "+" + tag, R.A, R.b, V.logVol + n * std::log(s), V.inradius * s, {}};
  for (const auto& v : V.verts) W.verts.push_back(M * v + t);
  return W;
}

std::vector<VBody> clipBodies(int n, unsigned seed)
{
  std::vector<VBody> out;
  std::vector<double> lo(n), hi(n);
  std::mt19937 g(seed);
  std::uniform_real_distribution<double> U(0.2, 3.0);
  for (int i = 0; i < n; ++i)
  {
    lo[i] = -U(g);
    hi[i] = U(g);
  }
  out.push_back(vbox(lo, hi));
  out.push_back(vbox(std::vector<double>(n, 0.0), std::vector<double>(n, 1.0)));
  // skinny: aspect 1e3 and 1e6
  for (double w : {1e-3, 1e-6})
  {
    std::vector<double> h(n, 1.0);
    h[0] = w;
    VBody S = vbox(std::vector<double>(n, 0.0), h);
    S.name = "skinny" + std::to_string(int(std::round(-std::log10(w))));
    out.push_back(S);
  }
  out.push_back(vsimplex(n));
  if (n <= 10) out.push_back(vcross(n));
  const std::size_t base = out.size();
  for (std::size_t k = 0; k < base; ++k)
  {
    const MT Q = randomRotation(n, seed * 31u + unsigned(k));
    VT t = VT::Constant(n, 0.37);
    out.push_back(vmap(out[k], Q, 1.0, t, "rot"));
  }
  // translated far, tiny, huge (rotated simplex and skinny box)
  const MT Q = randomRotation(n, seed + 7u);
  out.push_back(vmap(vsimplex(n), Q, 1.0, VT::Constant(n, 1e6), "far1e6"));
  out.push_back(vmap(vsimplex(n), Q, 1e-6, VT::Constant(n, 1e-3), "tiny1e-6"));
  out.push_back(vmap(vsimplex(n), Q, 1e6, VT::Constant(n, -1e5), "huge1e6"));
  out.push_back(vmap(out[2], Q, 1.0, VT::Constant(n, 5.0), "rot"));
  return out;
}

// ---------------------------------------------------------------- ttc runner
std::string findTtc()
{
  if (const char* e = std::getenv("TTC_BIN")) return e;
  for (const char* p : {"../../build/ttc", "build/ttc", "./ttc"})
    if (std::filesystem::exists(p)) return p;
  return {};
}

std::string runCapture(const std::string& cmd)
{
  std::string out;
  FILE* f = popen(cmd.c_str(), "r");
  if (!f) return out;
  char buf[4096];
  while (std::fgets(buf, sizeof buf, f)) out += buf;
  pclose(f);
  return out;
}

double grabNumber(const std::string& text, const std::string& pattern)
{
  std::smatch m;
  if (std::regex_search(text, m, std::regex(pattern))) return std::stod(m[1].str());
  return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace

// ============================================================================
TEST_CASE("[T4][thresh] lraSampleThreshold equals Alg. 1 line 5 on an (eps, delta, m) grid")
{
  int checked = 0, bad = 0;
  for (double eps : kEps)
    for (double delta : kDelta)
      for (double m : kM)
      {
        const double got = ttc::lraSampleThreshold(eps, delta, std::size_t(m));
        const long double want = paperThresh(eps, delta, m);
        ++checked;
        if (!(std::abs(got / double(want) - 1.0) <= 1e-12))
        {
          ++bad;
          MESSAGE("eps=" << eps << " delta=" << delta << " m=" << m << " got=" << got
                         << " want=" << double(want));
        }
      }
  MESSAGE("lraSampleThreshold: " << checked << " grid points, " << bad << " mismatches");
  CHECK(bad == 0);
  // m = 0 (no polytope) is treated as m = 1.
  CHECK(ttc::lraSampleThreshold(0.8, 0.2, 0) == ttc::lraSampleThreshold(0.8, 0.2, 1));
  // Monotone: more samples for smaller eps, smaller delta, larger m.
  for (double m : kM)
    for (double delta : kDelta)
      for (std::size_t i = 1; i < kEps.size(); ++i)
        CHECK(ttc::lraSampleThreshold(kEps[i], delta, std::size_t(m))
              <= ttc::lraSampleThreshold(kEps[i - 1], delta, std::size_t(m)));
  // Which branch is active at the defaults: the 15 ln(6m/delta) branch only
  // dominates for astronomically many polytopes.
  const double first = 24.0 * std::log(24.0 / 0.2) / ((1 - 0.8 / 12) * std::pow(0.8 / 12, 2));
  MESSAGE("eps=0.8 delta=0.2: Thresh = " << ttc::lraSampleThreshold(0.8, 0.2, 1)
          << " (first branch " << first << "); the m-branch needs ln m > "
          << first / 15.0 - std::log(6.0 / 0.2));
}

TEST_CASE("[T4][proof] Thresh and eps' = eps/12 satisfy the inequalities of lem:numlattice")
{
  int bad = 0;
  for (double eps : kEps)
    for (double delta : kDelta)
      for (double m : kM)
      {
        const long double T = ttc::lraSampleThreshold(eps, delta, std::size_t(m));
        const long double e = eps / 12.0L;
        // Bound 2 (Err_j): 4 exp(-(1-e')e'^2 T/24) <= delta/6.
        const long double errBound = 4.0L * std::exp(-(1.0L - e) * e * e * T / 24.0L);
        // Bound 1 (Noise_j): exp(-T/15) <= delta/(6m).
        const long double noiseBound = std::exp(-T / 15.0L);
        if (!(errBound <= (delta / 6.0L) * (1 + 1e-12L)) ||
            !(noiseBound <= delta / (6.0L * m) * (1 + 1e-12L)))
        {
          ++bad;
          MESSAGE("eps=" << eps << " delta=" << delta << " m=" << m << " err=" << double(errBound)
                         << " noise=" << double(noiseBound));
        }
        // Total failure probability delta/6 + delta/6 + delta/4 < delta.
      }
  CHECK(bad == 0);

  // Constant inequalities used in the proof of lem:numlattice (eps' = eps/12,
  // lattice error eps/8) and of thm:volumettc, on a fine grid of eps in (0,1].
  int badConst = 0;
  for (int k = 1; k <= 1000; ++k)
  {
    const double eps = k / 1000.0, e = eps / 12.0, eta = eps / 96.0;
    bool ok = true;
    ok &= std::pow(1 - e, 3) * (1 - eps / 8) >= 1 - eps / 2;
    ok &= std::pow(1 + e, 3) * (1 + eps / 8) <= 1 + eps / 2;
    ok &= (1 + e) * (1 + eps / 8) <= 4.0 / 3.0;
    ok &= (1 - e) * (1 - eps / 8) >= 0.5;
    // thm:volumettc with the eps/8 lattice error used in the invariant ...
    ok &= (1 - eps / 2) / (1 + eps / 8) >= 1 - eps;
    ok &= (1 + eps / 2) / (1 - eps / 8) <= 1 + eps;
    // ... and with eta = eps/96 of GetPrecision (eta <= eps/8).
    ok &= eta <= eps / 8;
    if (!ok)
    {
      ++badConst;
      MESSAGE("constant inequality fails at eps=" << eps);
    }
  }
  CHECK(badConst == 0);
}

TEST_CASE("[T4][wiring] ttc uses eps' = eps/12, delta' = delta/(4m), Thresh(eps,delta,m) and b")
{
  const std::string ttcBin = findTtc();
  if (ttcBin.empty())
  {
    MESSAGE("ttc binary not found (set TTC_BIN); skipping");
    return;
  }
  // Three disjoint unit squares plus a zero-volume cube (x = 7): m = 4 cubes,
  // 3 kept, r = 12 facets kept, gamma_min = 0.5, n = 2.
  const auto path = std::filesystem::temp_directory_path()
                    / ("ttc_union_params_" + std::to_string(::getpid()) + ".smt2");
  {
    std::ofstream f(path);
    f << "(set-logic QF_LRA)\n(declare-const x Real)\n(declare-const y Real)\n"
         "(assert (or (and (>= x 0) (<= x 1) (>= y 0) (<= y 1))\n"
         "            (and (>= x 2) (<= x 3) (>= y 0) (<= y 1))\n"
         "            (and (>= x 4) (<= x 5) (>= y 0) (<= y 1))\n"
         "            (and (= x 7) (>= y 0) (<= y 1))))\n(check-sat)\n";
  }
  // eps = 0.3 separates m = 4 (b = 7) from m = 3 (b = 6).
  REQUIRE(ttc::latticePrecision(2, 12, 4, 0.5, 0.3 / 96.0) != ttc::latticePrecision(2, 12, 3, 0.5, 0.3 / 96.0));
  for (auto [eps, delta] : std::vector<std::pair<double, double>>{{0.3, 0.2}, {0.8, 0.05}})
  {
    const std::string out = runCapture(ttcBin + " " + path.string() + " -v 2 --seed 1 -e "
                                       + std::to_string(eps) + " -d " + std::to_string(delta)
                                       + " 2>&1");
    const double m = grabNumber(out, R"(c polytopes:\s*(\S+))");
    const double thresh = grabNumber(out, R"(Starting union algorithm, threshold:\s*(\S+))");
    const double ePrime = grabNumber(out, R"(oracle targets: epsilon=([0-9.eE+-]+))");
    const double dPrime = grabNumber(out, R"(oracle targets: epsilon=[0-9.eE+-]+, delta=([0-9.eE+-]+))");
    const double b = grabNumber(out, R"(lattice precision \(b\):\s*(\S+))");
    const double digits = grabNumber(out, R"(gmp sampling precision \(digits\):\s*(\S+))");
    const double vol = grabNumber(out, R"(\ns volume\s+(\S+))");
    MESSAGE("eps=" << eps << " delta=" << delta << ": m=" << m << " Thresh=" << thresh
                   << " eps'=" << ePrime << " delta'=" << dPrime << " b=" << b
                   << " digits=" << digits << " volume=" << vol << " (exact 3)");
    REQUIRE(m == 4);
    CHECK(std::abs(thresh / ttc::lraSampleThreshold(eps, delta, 4) - 1) < 1e-5);
    CHECK(std::abs(ePrime / (eps / 12) - 1) < 1e-5);
    CHECK(std::abs(dPrime / (delta / (4 * 4)) - 1) < 1e-5);
    CHECK(b == ttc::latticePrecision(2, 12, 4, 0.5, eps / 96.0));
    CHECK(digits >= b + 32);
  }
  std::filesystem::remove(path);
}

TEST_CASE("[T4][clip] inner ball, outer radius and clip bounds are sound in dims 1..10")
{
  int bodies = 0, badBall = 0, badGamma = 0, weakGamma = 0, badR = 0, badBounds = 0, badClamp = 0;
  double worstGammaRatio = 1.0;
  for (int n = 1; n <= 10; ++n)
    for (const VBody& V : clipBodies(n, 1000u + unsigned(n)))
    {
      ++bodies;
      const std::string tag = V.name + " n=" + std::to_string(n);
      const ttc::InnerBall ball = ttc::certifiedInnerBall(V.A, V.b);
      REQUIRE_MESSAGE(ball.radius > 0, tag);
      const double gamma = ball.radius;
      const VT& c = ball.center;
      double R = 0;
      try
      {
        R = ttc::outerRadius(V.A, V.b, c);
      }
      catch (const std::exception& e)
      {
        FAIL_CHECK(tag << ": outerRadius threw " << e.what());
        continue;
      }
      // (a) B(c, gamma) subset P: slack_i >= gamma |A_i|.
      for (Eigen::Index i = 0; i < V.A.rows(); ++i)
      {
        const double an = V.A.row(i).norm();
        const double slack = V.b(i) - V.A.row(i).dot(c);
        const double tol = 1e-12 * (std::abs(V.b(i)) + an * c.norm()) + 1e-9 * gamma * an;
        if (slack < gamma * an - tol)
        {
          ++badBall;
          MESSAGE(tag << ": ball leaves P at row " << i << " slack=" << slack
                      << " gamma|a|=" << gamma * an);
          break;
        }
      }
      // (b) gamma <= exact inradius (soundness) and close to it (tightness).
      if (V.inradius > 0)
      {
        worstGammaRatio = std::min(worstGammaRatio, gamma / V.inradius);
        if (gamma > V.inradius * (1 + 1e-7))
        {
          ++badGamma;
          MESSAGE(tag << ": gamma=" << gamma << " > inradius " << V.inradius);
        }
        if (gamma < V.inradius * (1 - 1e-3))
        {
          ++weakGamma;
          MESSAGE(tag << ": gamma=" << gamma << " < 0.999 inradius " << V.inradius);
        }
      }
      // (c) R >= max_v |v - c|.
      double maxDist = 0;
      for (const auto& v : V.verts) maxDist = std::max(maxDist, (v - c).norm());
      if (R < maxDist * (1 - 1e-9))
      {
        ++badR;
        MESSAGE(tag << ": outerRadius " << R << " < max vertex distance " << maxDist);
      }
      // (d) (2 gamma/sqrt n)^n <= vol <= (2R)^n, in logs.
      const double logLo = n * std::log(2 * gamma / std::sqrt(double(n)));
      const double logHi = n * std::log(2 * R);
      if (!(logLo <= V.logVol + 1e-9) || !(V.logVol <= logHi + 1e-9))
      {
        ++badBounds;
        MESSAGE(tag << ": ln bounds [" << logLo << ", " << logHi << "] miss ln vol " << V.logVol);
      }
      // (e) ttc's clamp keeps a (1 +- e') oracle value inside (1 +- e') vol.
      for (double eps : {0.8, 0.3, 0.05})
      {
        const double e = eps / 12, vol = std::exp(V.logVol);
        for (double f : {1 - e, 1.0, 1 + e})
        {
          double out = 0;
          const bool ok = ttcClip(f * vol, eps, n, gamma, R, out);
          if (!ok || out < (1 - e) * vol * (1 - 1e-12) || out > (1 + e) * vol * (1 + 1e-12))
          {
            ++badClamp;
            MESSAGE(tag << ": clamp(" << f << " vol) = " << out << ", vol=" << vol);
          }
        }
      }
    }
  MESSAGE(bodies << " bodies; worst gamma/inradius = " << worstGammaRatio);
  CHECK(badBall == 0);
  CHECK(badGamma == 0);
  CHECK(weakGamma == 0);
  CHECK(badR == 0);
  CHECK(badBounds == 0);
  CHECK(badClamp == 0);
}

TEST_CASE("[T4][scale] clip bounds in double for dims up to 40 and extreme sides")
{
  // Unit-free: for a cube of side s in R^n, gamma = s/2 and the computed R (LP
  // extents) is s sqrt(n)/2; vol = s^n.  Whenever vol is a normal double ttc
  // must neither throw nor move a (1 +- e') oracle value out of (1 +- e') vol.
  int representable = 0, bad = 0, lpBad = 0;
  for (int n : {1, 2, 5, 10, 20, 30, 40})
    for (double s : {1e-9, 1e-7, 1e-3, 1.0, 1e3, 1e7, 1e9})
    {
      const double logVol = n * std::log(s);
      Body B = cube(n, 0.0, s);
      double gamma = s / 2, R = s * std::sqrt(double(n)) / 2;
      {
        // The real ttc components on the real body.
        const ttc::InnerBall ball = ttc::certifiedInnerBall(B.A, B.b);
        const double Rlp = ttc::outerRadius(B.A, B.b, ball.center);
        // gamma: sound (<= exact) and within 0.1% (certifiedInnerBall's contract).
        if (ball.radius > gamma * (1 + 1e-7) || ball.radius < gamma * (1 - 1e-3)
            || Rlp < R * (1 - 1e-9))
        {
          ++lpBad;
          MESSAGE("n=" << n << " s=" << s << ": gamma=" << ball.radius << " (exact " << gamma
                       << "), R=" << Rlp << " (circumradius " << R << ")");
        }
        gamma = ball.radius;
        R = Rlp;
      }
      double lower = 0, upper = 0;
      const bool normal = logVol > std::log(std::numeric_limits<double>::min())
                          && logVol < std::log(std::numeric_limits<double>::max() / 2);
      if (!normal)
      {
        double out = 0;
        const bool ok = ttcClip(std::exp(logVol), 0.8, n, gamma, R, out, &lower, &upper);
        MESSAGE("n=" << n << " s=" << s << ": vol=e^" << logVol << " not a normal double;"
                     << " ttc clip " << (ok ? "returns " + std::to_string(out) : "throws")
                     << " (lower=" << lower << ", upper=" << upper << ")");
        continue;
      }
      ++representable;
      const double vol = std::exp(logVol);
      for (double eps : {0.8, 0.05})
      {
        const double e = eps / 12;
        for (double f : {1 - e, 1.0, 1 + e})
        {
          double out = 0;
          const bool ok = ttcClip(f * vol, eps, n, gamma, R, out, &lower, &upper);
          if (!ok || out < (1 - e) * vol * (1 - 1e-12) || out > (1 + e) * vol * (1 + 1e-12))
          {
            ++bad;
            MESSAGE("n=" << n << " s=" << s << " eps=" << eps << ": clip(" << f * vol << ") -> "
                         << (ok ? std::to_string(out) : std::string("throw")) << " lower=" << lower
                         << " upper=" << upper);
          }
        }
      }
    }
  MESSAGE(representable << " representable (n, s) pairs checked");
  CHECK(lpBad == 0);
  CHECK(bad == 0);
}

TEST_CASE("[T4][scale] log2Scale bookkeeping of estimateUnionVolume for vol 1e-300..1e300, n <= 40")
{
  // One polytope (or two disjoint ones), exact volumes, a generator that
  // returns the polytope's centre: the estimator reduces to N ~ Pois(t p),
  // est = N / (p C), so est/vol - 1 has sd 1/sqrt(t p) <= sqrt(2/T).
  const double eps = 0.8, delta = 0.2;
  int runs = 0, bad = 0;
  for (int n : {1, 2, 10, 40})
    for (double logVol10 : {-300.0, -100.0, -10.0, 0.0, 10.0, 100.0, 300.0})
    {
      // Side of the cube with this volume, and ttc's b for it (gamma = s/2).
      const double s = std::pow(10.0, logVol10 / n);
      if (!std::isfinite(s) || s <= 0) continue;
      const int b = ttc::latticePrecision(n, 2 * n, 1, s / 2, eps / 96.0);
      const double thresh = ttc::lraSampleThreshold(eps, delta, 2);
      const double log2Scale = std::log2(thresh) + double(b) * n * std::log2(10.0);
      REQUIRE(std::isfinite(log2Scale));
      // Two disjoint polytopes: P1 of volume 10^logVol10 (a "cube" of side s at
      // the origin, membership only matters for the second), P2 a unit-free
      // tiny or equal copy far away.
      for (double ratio : {1.0, 1e-30})
      {
        if (logVol10 - 30 < -290 && ratio < 1) continue;  // 2nd volume not a normal double
        std::vector<ttc::PolytopeVolume> kept(2);
        for (int k = 0; k < 2; ++k)
        {
          Body B = cube(n, 10.0 * k * std::max(1.0, s), 10.0 * k * std::max(1.0, s) + s);
          kept[k] = {B.A, B.b, std::pow(10.0, logVol10) * (k ? ratio : 1.0)};
        }
        auto generate = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& bb, long long cnt,
                            double& t) {
          t = 0;
          // centre of the box: (hi + lo)/2 per coordinate, rows (+e_i, -e_i).
          Eigen::VectorXd c(A.cols());
          for (Eigen::Index i = 0; i < A.cols(); ++i) c(i) = 0.5 * (bb(2 * i) - bb(2 * i + 1));
          return std::vector<Eigen::VectorXd>(std::size_t(cnt), c);
        };
        ttc::VolumeComputationResult res;
        ttc::estimateUnionVolume<double>(kept, thresh, log2Scale, 17u + unsigned(runs), ttc::storedVolumes(kept), generate,
                                       res, {});
        const double exact = kept[0].volume + kept[1].volume;
        const double rel = res.volumeEstimate / exact - 1;
        const double sd = std::sqrt(2.0 / thresh);
        ++runs;
        if (!std::isfinite(res.volumeEstimate) || std::abs(rel) > 5 * sd)
        {
          ++bad;
          MESSAGE("n=" << n << " vol=1e" << logVol10 << " ratio=" << ratio << " b=" << b
                       << " log2Scale=" << log2Scale << ": est=" << res.volumeEstimate
                       << " rel=" << rel << " (5 sd = " << 5 * sd << ")");
        }
      }
    }
  MESSAGE(runs << " runs");
  CHECK(bad == 0);
}
