// Uniformity of ttc's single-polytope sampler (src/volume/polytope_sampler.cpp,
// sampleRestartWalk: every sample is the endpoint of a fresh walk of
// 10 + n^2 steps started at the certified Chebyshev centre) on shapes whose
// uniform distribution is known exactly, dims 1..10, for the default Billiard
// walk and for CDHR and RDHR.
//
// Statistics.  Every check reduces to a p-value of an exact null:
//   * KS: a statistic s with known CDF F is mapped to u = F(s) ~ U(0,1) and
//     tested with the one-sample Kolmogorov-Smirnov test (asymptotic p-value
//     with Stephens' correction; N = 4000 >> 35, so the asymptotics are exact
//     to well below the significance level).
//   * chi-square goodness of fit on cells of exactly known probability (all
//     expected counts >= 15).
//   * vertex occupancy: for a vertex v and 0 < c < 1 the homothet
//     v + c (P - v) lies in P and has volume fraction exactly c^n, so its count
//     is Binomial(N, c^n); two-sided normal-approximation p-value (expected
//     count >= 40).
//   * radial ("gauge") test: for an interior point z the gauge
//     g(x) = max_j A_j (x - z) / (b_j - A_j z) of a uniform point satisfies
//     P(g <= t) = t^n (star-shaped homothety about z), for ANY polytope.
// Samples of sampleRestartWalk are independent (every one restarts at the
// centre on a continuing RNG stream), so these i.i.d. tests apply.
// Each individual check uses alpha = 1e-4 (the default suite runs a few
// hundred checks: family-wise false-alarm rate <~ 5%; seeds are fixed, so the
// outcome is deterministic).  Failures here are sampler bugs, not noise.
//
// Default run ~45 s.  The walk-length / rounding sweeps are tagged "slow":
//   test_sampling_shapes -ts=slow --no-skip
// Findings (2026-11-14, N = 40000 corner test unless noted):
//   * Billiard is too central at EVERY walk length (gauge D+ ~ 0.009-0.015,
//     D- ~ 0, walk 46..1760, rounded or not): volesti reflects at 0.995 of the
//     distance to the facet (uniform_billiard_walk.hpp `dl = 0.995`), so each
//     bounce skips a skin next to the facet.  dl = 0.999999 removes it.
//   * RDHR under-samples the acute simplex vertices e_i at 10 + n^2 steps
//     (obs/exp 0.30-0.65 for the 0.2%-1% homothets); x4 the walk fixes it.
//   * CDHR in rounded coordinates (ttc's default since rounding) loses its
//     axis alignment and is too central on the simplex; x4 the walk fixes it.
//   * Unrounded (e8a8b0a), billiard and RDHR fail on boxes with aspect >= 4 and
//     on simplex2 x simplex8: the Chebyshev centre is not unique and the walks
//     never leave its neighbourhood; rounding (or x4-x16 the walk) fixes it.
#include "doctest.h"
#include "common.hpp"

#include <boost/math/special_functions/gamma.hpp>
#include <chrono>
#include <cstdio>
#include <functional>
#include <map>


using namespace vt;

namespace
{

constexpr double kAlpha = 1e-4;
constexpr long long kN = 4000;

const char* walkName(ttc::SamplerWalk w)
{
  switch (w)
  {
    case ttc::SamplerWalk::Billiard: return "billiard";
    case ttc::SamplerWalk::AcceleratedBilliard: return "accel-billiard";
    case ttc::SamplerWalk::Ball: return "ball";
    case ttc::SamplerWalk::RDHR: return "rdhr";
    case ttc::SamplerWalk::CDHR: return "cdhr";
  }
  return "?";
}

const ttc::SamplerWalk kWalks[] = {ttc::SamplerWalk::Billiard, ttc::SamplerWalk::CDHR,
                                   ttc::SamplerWalk::RDHR};

// ------------------------------------------------------------- statistics
// Kolmogorov distribution tail Q(lambda) = 2 sum (-1)^{k-1} exp(-2 k^2 lambda^2).
double kolmogorovQ(double lambda)
{
  if (lambda < 0.2) return 1.0;
  double s = 0.0;
  for (int k = 1; k <= 100; ++k)
  {
    const double term = std::exp(-2.0 * k * k * lambda * lambda);
    s += (k % 2 ? 2.0 : -2.0) * term;
    if (term < 1e-16) break;
  }
  return std::min(1.0, std::max(0.0, s));
}

struct KS
{
  double D = 0, p = 1, dPlus = 0, dMinus = 0;
};

// u: values that are U(0,1) under the null.
KS ksUniform(std::vector<double> u)
{
  std::sort(u.begin(), u.end());
  const double N = u.size();
  KS r;
  for (size_t i = 0; i < u.size(); ++i)
  {
    r.dPlus = std::max(r.dPlus, (i + 1) / N - u[i]);   // sample CDF above F
    r.dMinus = std::max(r.dMinus, u[i] - i / N);       // sample CDF below F
  }
  r.D = std::max(r.dPlus, r.dMinus);
  const double en = std::sqrt(N);
  r.p = kolmogorovQ((en + 0.12 + 0.11 / en) * r.D);
  return r;
}

double chiSquareP(const std::vector<long long>& counts, const std::vector<double>& probs,
                  double* stat = nullptr)
{
  long long N = 0;
  for (auto c : counts) N += c;
  double x2 = 0;
  for (size_t i = 0; i < counts.size(); ++i)
  {
    const double e = N * probs[i];
    x2 += (counts[i] - e) * (counts[i] - e) / e;
  }
  if (stat) *stat = x2;
  const double df = counts.size() - 1.0;
  return boost::math::gamma_q(df / 2.0, x2 / 2.0);
}

double binomialP(long long k, long long N, double p)
{
  const double z = (k - N * p) / std::sqrt(N * p * (1 - p));
  return std::erfc(std::abs(z) / std::sqrt(2.0));
}

// Collects the checks of one (shape, dim, walk) run.
struct Checker
{
  std::string tag;
  int checks = 0, fails = 0;
  double worstP = 1;
  std::string worst;
  bool assert_ = true;

  void check(const std::string& what, double p, const std::string& detail)
  {
    ++checks;
    const bool ok = p > kAlpha;
    if (!ok) ++fails;
    if (p < worstP)
    {
      worstP = p;
      worst = what + " " + detail;
    }
    if (!ok && assert_) std::printf("    FAIL %-40s p=%.2e  %s\n", (tag + " " + what).c_str(), p, detail.c_str());
    if (assert_) CHECK_MESSAGE(ok, tag << " " << what << " p=" << p << " " << detail);
  }

  void ks(const std::string& what, const std::vector<double>& u)
  {
    const KS r = ksUniform(u);
    char buf[160];
    std::snprintf(buf, sizeof buf, "D=%.4f (D+=%.4f D-=%.4f)", r.D, r.dPlus, r.dMinus);
    check("KS " + what, r.p, buf);
  }

  void summary() const
  {
    std::printf("  %-34s checks=%3d fails=%3d worst p=%.1e [%s]\n", tag.c_str(), checks, fails,
                worstP, worst.c_str());
  }
};

// ------------------------------------------------------------- generic tests
// Gauge of x about the interior point z.
double gauge(const MT& A, const VT& b, const VT& z, const VT& x)
{
  const VT num = A * (x - z), den = b - A * z;
  double g = 0;
  for (int j = 0; j < A.rows(); ++j) g = std::max(g, num(j) / den(j));
  return g;
}

void gaugeTest(Checker& C, const std::string& what, const MT& A, const VT& b, const VT& z,
               const std::vector<VT>& pts)
{
  const int n = A.cols();
  std::vector<double> u;
  u.reserve(pts.size());
  for (const auto& x : pts) u.push_back(std::pow(std::min(1.0, gauge(A, b, z, x)), n));
  const KS r = ksUniform(u);
  // D+ > 0 with D- ~ 0: too few samples near the boundary (too central).
  std::printf("    gauge^n about %s: D+=%.4f D-=%.4f p=%.1e\n", what.c_str(), r.dPlus,
              r.dMinus, r.p);
  C.ks("gauge^n about " + what, u);
}

bool inside(const MT& A, const VT& b, const VT& x, double tol = 1e-12)
{
  return ((A * x - b).array() <= tol * (1.0 + b.array().abs())).all();
}

// Counts in the vertex homothets v + c (P - v), c = frac^(1/n); each is
// Binomial(N, frac).  Prints the worst observed/expected ratio.
void vertexTest(Checker& C, const MT& A, const VT& b, const std::vector<VT>& vertices,
                const std::vector<VT>& pts, double frac)
{
  const int n = A.cols();
  const double c = std::pow(frac, 1.0 / n);
  double worstRatio = 1;
  for (size_t k = 0; k < vertices.size(); ++k)
  {
    const VT& v = vertices[k];
    long long cnt = 0;
    for (const auto& x : pts) cnt += inside(A, b, v + (x - v) / c);
    const double ratio = cnt / (frac * pts.size());
    if (std::abs(ratio - 1) > std::abs(worstRatio - 1)) worstRatio = ratio;
    char buf[160];
    std::snprintf(buf, sizeof buf, "count=%lld expected=%.0f ratio=%.2f", cnt, frac * pts.size(),
                  ratio);
    char what[64];
    std::snprintf(what, sizeof what, "vertex%zu frac=%.3f", k, frac);
    C.check(what, binomialP(cnt, pts.size(), frac), buf);
  }
  std::printf("    vertex homothets frac=%.3f (c=%.3f): worst obs/exp = %.2f over %zu vertices\n",
              frac, c, worstRatio, vertices.size());
}

// ------------------------------------------------------------- shapes
// A shape is a body in its own coordinates (A, b, vertices, centre z) plus the
// affine map x = Q y + t to the coordinates the sampler sees.
struct Shape
{
  std::string name;
  MT A;
  VT b;
  VT z;  // interior point used for the gauge test (Chebyshev centre)
  std::vector<VT> vertices;
  MT Q;
  VT t;
  std::function<void(Checker&, const std::vector<VT>&)> extra;  // shape-specific tests
};

using Sampler = std::function<std::vector<VT>(const MT&, const VT&, long long, unsigned)>;

Sampler ttcSampler(ttc::SamplerWalk w, int walkLength = -1)
{
  return [=](const MT& A, const VT& b, long long N, unsigned seed) {
    return ttcProposals(A, b, N, seed, w, walkLength);
  };
}

// Runs all tests of the shape on N samples; returns number of failed checks.
int runShape(const Shape& S, const std::string& samplerTag, const Sampler& sampler,
             unsigned seed, bool assertChecks = true, long long N = kN)
{
  const int n = S.A.cols();
  // Body as the sampler sees it: x = Q y + t  <=>  y = Q^{-1}(x - t).
  MT Qinv = S.Q.inverse();
  MT Ax = S.A * Qinv;
  VT bx = S.b + Ax * S.t;
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<VT> xs = sampler(Ax, bx, N, seed);
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  Checker C;
  C.tag = S.name + " " + samplerTag;
  C.assert_ = assertChecks;
  REQUIRE(static_cast<long long>(xs.size()) == N);
  if (assertChecks)
  {
    // Where every walk restarts: ttc's certified Chebyshev centre (own coords).
    const ttc::InnerBall ball = ttc::certifiedInnerBall(Ax, bx);
    const VT c = Qinv * (ball.center - S.t);
    std::printf("    %s: Chebyshev centre (own coords) [", S.name.c_str());
    for (int i = 0; i < n; ++i) std::printf("%s%.3f", i ? " " : "", c(i));
    std::printf("] r=%.3f\n", ball.radius);
  }
  std::vector<VT> ys;
  ys.reserve(xs.size());
  long long outside = 0;
  for (const auto& x : xs)
  {
    ys.push_back(Qinv * (x - S.t));
    outside += !inside(S.A, S.b, ys.back(), 1e-9);
  }
  C.check("all samples inside P", outside == 0 ? 1.0 : 0.0,
          "outside=" + std::to_string(outside));
  gaugeTest(C, "centre", S.A, S.b, S.z, ys);
  if (!S.vertices.empty())
  {
    // A point 90% of the way from the centre to vertex 0: gauge about an
    // off-centre point probes the far side of the body.
    const VT z1 = S.z + 0.9 * (S.vertices[0] - S.z);
    gaugeTest(C, "near-vertex0", S.A, S.b, z1, ys);
    // Expected counts >= 40 per vertex homothet.
    for (double frac : {0.05, 0.01, 0.002})
      if (frac * N >= 40) vertexTest(C, S.A, S.b, S.vertices, ys, frac);
  }
  if (S.extra) S.extra(C, ys);
  C.summary();
  std::printf("    (%.2fs)\n", secs);
  return C.fails;
}

int runAllWalks(const Shape& S, unsigned seed)
{
  int fails = 0;
  for (auto w : kWalks) fails += runShape(S, walkName(w), ttcSampler(w), seed);
  return fails;
}

void identityMap(Shape& S)
{
  const int n = S.A.cols();
  S.Q = MT::Identity(n, n);
  S.t = VT::Zero(n);
}

// Standard simplex {y >= 0, sum y <= 1}: (y, 1 - sum y) ~ Dirichlet(1,...,1),
// so y_i ~ Beta(1, n) (CDF 1 - (1-t)^n) and sum y has CDF t^n.
Shape simplexShape(int n)
{
  Body B = simplex(n);
  Shape S;
  S.name = "simplex" + std::to_string(n);
  S.A = B.A;
  S.b = B.b;
  // Chebyshev centre: r (1,...,1) with r = 1 / (n + sqrt n).
  S.z = VT::Constant(n, 1.0 / (n + std::sqrt(double(n))));
  S.vertices.push_back(VT::Zero(n));
  for (int i = 0; i < n; ++i) S.vertices.push_back(VT::Unit(n, i));
  identityMap(S);
  S.extra = [n](Checker& C, const std::vector<VT>& ys) {
    std::vector<double> us;
    for (const auto& y : ys) us.push_back(std::pow(std::min(1.0, y.sum()), n));
    C.ks("sum ~ Beta(n,1)", us);
    for (int i = 0; i < n; ++i)
    {
      std::vector<double> u;
      for (const auto& y : ys) u.push_back(1 - std::pow(1 - y(i), n));
      C.ks("y" + std::to_string(i) + " ~ Beta(1,n)", u);
    }
  };
  return S;
}

// Cross-polytope {|y|_1 <= 1}: |y|_1 has CDF t^n, |y_i| has CDF 1 - (1-t)^n,
// and the sign pattern is uniform over the 2^n orthants.
Shape crossShape(int n)
{
  Body B = crossPolytope(n);
  Shape S;
  S.name = "cross" + std::to_string(n);
  S.A = B.A;
  S.b = B.b;
  S.z = VT::Zero(n);
  for (int i = 0; i < n; ++i)
  {
    S.vertices.push_back(VT::Unit(n, i));
    S.vertices.push_back(-VT::Unit(n, i));
  }
  identityMap(S);
  S.extra = [n](Checker& C, const std::vector<VT>& ys) {
    for (int i = 0; i < n; ++i)
    {
      std::vector<double> u;
      for (const auto& y : ys) u.push_back(1 - std::pow(1 - std::min(1.0, std::abs(y(i))), n));
      C.ks("|y" + std::to_string(i) + "| ~ Beta(1,n)", u);
    }
    const int cells = 1 << n;
    if (double(ys.size()) / cells >= 15)
    {
      std::vector<long long> cnt(cells, 0);
      for (const auto& y : ys)
      {
        int k = 0;
        for (int i = 0; i < n; ++i) k |= (y(i) > 0) << i;
        ++cnt[k];
      }
      double x2;
      const double p = chiSquareP(cnt, std::vector<double>(cells, 1.0 / cells), &x2);
      char buf[96];
      std::snprintf(buf, sizeof buf, "X2=%.1f df=%d", x2, cells - 1);
      C.check("orthant chi-square", p, buf);
    }
  };
  return S;
}

// Cube [-1,1]^n under a random rotation; tested in the cube's own
// coordinates: y_i ~ U(-1,1) and max |y_i| has CDF t^n (= gauge about 0).
Shape rotatedCubeShape(int n, bool rotate = true)
{
  Body B = cube(n, -1, 1);
  Shape S;
  S.name = (rotate ? "rotcube" : "cube") + std::to_string(n);
  S.A = B.A;
  S.b = B.b;
  S.z = VT::Zero(n);
  // A few vertices: all-plus, all-minus, alternating signs and its negation.
  VT alt(n);
  for (int i = 0; i < n; ++i) alt(i) = i % 2 ? -1.0 : 1.0;
  S.vertices = {VT::Ones(n), -VT::Ones(n), alt, -alt};
  if (n == 1) S.vertices.resize(2);
  S.Q = rotate ? randomRotation(n, 1000 + n) : MT::Identity(n, n);
  S.t = rotate ? VT::Constant(n, 0.5) : VT::Zero(n);
  S.extra = [n](Checker& C, const std::vector<VT>& ys) {
    for (int i = 0; i < n; ++i)
    {
      std::vector<double> u;
      for (const auto& y : ys) u.push_back(std::clamp((y(i) + 1) / 2, 0.0, 1.0));
      C.ks("y" + std::to_string(i) + " ~ U(-1,1)", u);
    }
  };
  return S;
}

// Axis-aligned box prod [0, L_i], L_i = 1 + i (aspect ratio up to n):
// x_i / L_i ~ U(0,1) independently.
Shape boxShape(int n)
{
  std::vector<double> lo(n, 0.0), hi(n);
  for (int i = 0; i < n; ++i) hi[i] = 1.0 + i;
  Body B = box(lo, hi);
  Shape S;
  S.name = "box1.." + std::to_string(n);
  S.A = B.A;
  S.b = B.b;
  S.z = VT::Zero(n);
  for (int i = 0; i < n; ++i) S.z(i) = hi[i] / 2;
  VT top(n);
  for (int i = 0; i < n; ++i) top(i) = hi[i];
  S.vertices = {VT::Zero(n), top};
  identityMap(S);
  S.extra = [n, hi](Checker& C, const std::vector<VT>& ys) {
    for (int i = 0; i < n; ++i)
    {
      std::vector<double> u;
      for (const auto& y : ys) u.push_back(std::clamp(y(i) / hi[i], 0.0, 1.0));
      C.ks("x" + std::to_string(i) + "/L ~ U(0,1)", u);
    }
  };
  return S;
}

// Regular k-gon with circumradius 1, vertex j at angle 2 pi j / k.  By
// symmetry the 2k half-sectors (vertex ray to edge-midpoint ray) have equal
// area; gauge about 0 has CDF t^2.
Shape polygonShape(int k)
{
  Shape S;
  S.name = std::to_string(k) + "-gon";
  S.A = MT(k, 2);
  S.b = VT::Constant(k, std::cos(M_PI / k));
  for (int j = 0; j < k; ++j)
  {
    const double a = 2 * M_PI * (j + 0.5) / k;
    S.A(j, 0) = std::cos(a);
    S.A(j, 1) = std::sin(a);
    VT v(2);
    v << std::cos(2 * M_PI * j / k), std::sin(2 * M_PI * j / k);
    S.vertices.push_back(v);
  }
  S.z = VT::Zero(2);
  identityMap(S);
  S.extra = [k](Checker& C, const std::vector<VT>& ys) {
    std::vector<long long> cnt(2 * k, 0);
    for (const auto& y : ys)
    {
      double phi = std::atan2(y(1), y(0));
      if (phi < 0) phi += 2 * M_PI;
      ++cnt[std::min(2 * k - 1, int(phi / (M_PI / k)))];
    }
    double x2;
    const double p = chiSquareP(cnt, std::vector<double>(2 * k, 0.5 / k), &x2);
    char buf[96];
    std::snprintf(buf, sizeof buf, "X2=%.1f df=%d", x2, 2 * k - 1);
    C.check("half-sector chi-square", p, buf);
  };
  return S;
}

// Delta_a x Delta_b: the two factors are independent uniform simplices, so
// U1 = (sum u)^a, U2 = (sum w)^b are independent U(0,1): KS on each and a
// chi-square on the 4x4 grid of (U1, U2) (each cell 1/16).
Shape simplexProductShape(int a, int bdim)
{
  const int n = a + bdim;
  Body Sa = simplex(a), Sb = simplex(bdim);
  Shape S;
  S.name = "simplex" + std::to_string(a) + "xsimplex" + std::to_string(bdim);
  S.A = MT::Zero(a + 1 + bdim + 1, n);
  S.b = VT::Zero(a + 1 + bdim + 1);
  S.A.block(0, 0, a + 1, a) = Sa.A;
  S.b.head(a + 1) = Sa.b;
  S.A.block(a + 1, a, bdim + 1, bdim) = Sb.A;
  S.b.tail(bdim + 1) = Sb.b;
  // Chebyshev centre: the inner radius is min(r_a, r_b); the factor with the
  // larger inradius can sit anywhere, use its own centre (still interior).
  S.z = VT(n);
  S.z.head(a) = VT::Constant(a, 1.0 / (a + std::sqrt(double(a))));
  S.z.tail(bdim) = VT::Constant(bdim, 1.0 / (bdim + std::sqrt(double(bdim))));
  for (int i = -1; i < a; ++i)
    for (int j = -1; j < bdim; ++j)
    {
      VT v = VT::Zero(n);
      if (i >= 0) v(i) = 1;
      if (j >= 0) v(a + j) = 1;
      S.vertices.push_back(v);
    }
  identityMap(S);
  S.extra = [a, bdim](Checker& C, const std::vector<VT>& ys) {
    std::vector<double> u1, u2;
    std::vector<long long> grid(16, 0);
    for (const auto& y : ys)
    {
      const double s1 = std::pow(std::min(1.0, y.head(a).sum()), a);
      const double s2 = std::pow(std::min(1.0, y.tail(bdim).sum()), bdim);
      u1.push_back(s1);
      u2.push_back(s2);
      ++grid[std::min(3, int(4 * s1)) * 4 + std::min(3, int(4 * s2))];
    }
    C.ks("factor1 sum ~ Beta(a,1)", u1);
    C.ks("factor2 sum ~ Beta(b,1)", u2);
    double x2;
    const double p = chiSquareP(grid, std::vector<double>(16, 1.0 / 16), &x2);
    char buf[96];
    std::snprintf(buf, sizeof buf, "X2=%.1f df=15", x2);
    C.check("4x4 joint (independence) chi-square", p, buf);
  };
  return S;
}

// ttc's sampler with rounding explicitly on (ttc's default: the walk runs in
// the max-inscribed-ellipsoid-rounded body, points mapped back x = shift + T y)
// or off (the walk runs in the original coordinates, as at e8a8b0a).
Sampler roundingSampler(ttc::SamplerWalk w, int walkLength, bool round)
{
  return [=](const MT& A, const VT& b, long long N, unsigned seed) {
    double ballSeconds = 0;
    return ttc::sampleRestartWalk(A, b, N, w, unsigned(walkLength), seed, ballSeconds, round);
  };
}

}  // namespace

// ================================================================== tests
TEST_CASE("sampling shapes: axis-aligned cube and box, dims 1..10")
{
  int fails = 0;
  for (int n = 1; n <= 10; ++n)
  {
    fails += runAllWalks(rotatedCubeShape(n, false), 11 + n);
    fails += runAllWalks(boxShape(n), 21 + n);
  }
  std::printf("cube/box: total failed checks %d\n", fails);
}

TEST_CASE("sampling shapes: standard simplex, dims 2..10")
{
  int fails = 0;
  for (int n : {2, 3, 4, 6, 8, 10}) fails += runAllWalks(simplexShape(n), 31 + n);
  std::printf("simplex: total failed checks %d\n", fails);
}

TEST_CASE("sampling shapes: cross-polytope, dims 2..8")
{
  int fails = 0;
  for (int n : {2, 3, 4, 6, 8}) fails += runAllWalks(crossShape(n), 41 + n);
  std::printf("cross: total failed checks %d\n", fails);
}

TEST_CASE("sampling shapes: rotated cube, dims 2..10")
{
  int fails = 0;
  for (int n : {2, 3, 5, 10}) fails += runAllWalks(rotatedCubeShape(n), 51 + n);
  std::printf("rotated cube: total failed checks %d\n", fails);
}

TEST_CASE("sampling shapes: regular hexagon and octagon")
{
  int fails = 0;
  for (int k : {6, 8}) fails += runAllWalks(polygonShape(k), 61 + k);
  std::printf("polygons: total failed checks %d\n", fails);
}

TEST_CASE("sampling shapes: product of two simplices")
{
  int fails = 0;
  for (auto ab : {std::pair{2, 3}, std::pair{3, 5}, std::pair{2, 8}})
    fails += runAllWalks(simplexProductShape(ab.first, ab.second), 71 + ab.first * ab.second);
  std::printf("simplex products: total failed checks %d\n", fails);
}

TEST_CASE("sampling shapes: sharp corners with N = 40000 (vertex homothets of mass 0.2%)")
{
  // At N = 40000 the 0.2% vertex homothets (c = 0.002^(1/n): 0.35 in 6-D,
  // 0.54 in 10-D) expect 80 samples each, and the KS tests resolve
  // deviations of ~0.011 in CDF.
  int fails = 0;
  for (const Shape& S : {simplexShape(6), simplexShape(10), crossShape(6)})
    for (auto w : kWalks) fails += runShape(S, walkName(w), ttcSampler(w), 81, true, 40000);
  std::printf("sharp corners: total failed checks %d\n", fails);
}

// ------------------------------------------------------------------- slow
// Diagnostics (no assertions on ttc's defaults; each run prints its number of
// failed checks): how the failures above react to longer walks (x4, x16 the
// default 10 + n^2), with and without rounding (roundingSampler).
namespace
{
void sweep(const std::vector<Shape>& shapes, long long N, std::vector<int> roundedMults)
{
  for (const auto& S : shapes)
  {
    const int n = S.A.cols();
    for (auto w : kWalks)
      for (int mult : {1, 4, 16})
      {
        const int wl = mult * int(ttcWalkLength(n));
        const std::string tag = std::string(walkName(w)) + " wl=" + std::to_string(wl);
        runShape(S, "unrounded " + tag, roundingSampler(w, wl, false), 91, false, N);
        if (std::count(roundedMults.begin(), roundedMults.end(), mult))
          runShape(S, "rounded " + tag, roundingSampler(w, wl, true), 91, false, N);
      }
  }
}
}  // namespace

// Run with: test_sampling_shapes -ts=slow --no-skip
TEST_CASE("sampling shapes sweep: elongated bodies, walk length and rounding (N = 4000)"
          * doctest::test_suite("slow") * doctest::skip())
{
  sweep({boxShape(10), boxShape(6), simplexProductShape(2, 8), simplexShape(10)}, kN,
        {1, 4, 16});
}

TEST_CASE("sampling shapes sweep: sharp corners, walk length and rounding (N = 40000)"
          * doctest::test_suite("slow") * doctest::skip())
{
  sweep({simplexShape(6), simplexShape(10), crossShape(6)}, 40000, {1, 4});
}
