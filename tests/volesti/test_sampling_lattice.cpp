// Audit of ttc::sampleLattice (src/volume/polytope_sampler.hpp): GenerateSamples
// of the manuscript -- inflate every facet by 10^-b * kappa * |A_j|, draw
// proposals from the inflated body, round each coordinate stochastically to
// the lattice h Z (h = 10^-b), reject against P.
//
// Bodies are given with integer facet normals and right-hand sides b = B/10^q,
// so membership of a lattice point k*h is decided exactly in integers.
//
//   [rounding]   a fixed proposal is rounded to one of its two neighbouring
//                lattice values, unbiasedly (z-test, |z| < 4.5 per coordinate).
//   [exact]      with EXACT uniform proposals (rejection from a bounding box)
//                the output must be uniform over P cap hZ^n: chi-square over the
//                lattice points (p > 0.001), boundary/corner points neither
//                over- nor under-represented (binomial z, |z| < 4), no point
//                outside P.
//   [budget]     rejection rate vs the paper's bound eps/3 and the retry budget.
//   [walks]      composed with ttc's restart walks (ttcLatticeSamples): uniform
//                over the lattice points?  (chi-square, p > 0.001)
//   [seeding]    rounding RNG vs walk RNG.
#include "doctest.h"
#include "common.hpp"

#include <boost/math/distributions/chi_squared.hpp>
#include <cstdio>
#include <map>
#include <stdexcept>

using namespace vt;

namespace
{

// ------------------------------------------------------------------ bodies
// {x : A x <= B / 10^q}, A integer.
struct IBody
{
  std::string name;
  std::vector<std::vector<long long>> A;
  std::vector<long long> B;
  int q = 0;
  int dim() const { return int(A[0].size()); }
  MT Ad() const
  {
    MT M(A.size(), dim());
    for (size_t i = 0; i < A.size(); ++i)
      for (int j = 0; j < dim(); ++j) M(i, j) = double(A[i][j]);
    return M;
  }
  VT bd() const
  {
    VT v(B.size());
    for (size_t i = 0; i < B.size(); ++i) v(i) = double(B[i]) / std::pow(10.0, q);
    return v;
  }
};

long long ipow10(int e)
{
  long long r = 1;
  while (e-- > 0) r *= 10;
  return r;
}

// Exact membership of the lattice point k * 10^-p (requires q >= p).
bool memberExact(const IBody& P, const std::vector<long long>& k, int p)
{
  const long long s = ipow10(P.q - p);
  for (size_t i = 0; i < P.A.size(); ++i)
  {
    long long dot = 0;
    for (int j = 0; j < P.dim(); ++j) dot += P.A[i][j] * k[j];
    if (dot * s > P.B[i]) return false;
  }
  return true;
}

// Integer slack of facet i at k, in units of 10^-q.
long long slack(const IBody& P, size_t i, const std::vector<long long>& k, int p)
{
  long long dot = 0;
  for (int j = 0; j < P.dim(); ++j) dot += P.A[i][j] * k[j];
  return P.B[i] - dot * ipow10(P.q - p);
}

// Box prod [lo_j, hi_j], ends given in units of 10^-q.
IBody ibox(const std::vector<long long>& lo, const std::vector<long long>& hi, int q,
           std::string name)
{
  IBody P;
  P.name = std::move(name);
  P.q = q;
  const int n = int(lo.size());
  for (int j = 0; j < n; ++j)
  {
    std::vector<long long> r(n, 0);
    r[j] = 1;
    P.A.push_back(r);
    P.B.push_back(hi[j]);
    r[j] = -1;
    P.A.push_back(r);
    P.B.push_back(-lo[j]);
  }
  return P;
}

IBody addFacet(IBody P, std::vector<long long> a, long long B)
{
  P.A.push_back(std::move(a));
  P.B.push_back(B);
  return P;
}

// Bounding box of {A x <= b} by enumerating vertices (n <= 3).
std::pair<VT, VT> bbox(const MT& A, const VT& b)
{
  const int n = int(A.cols()), m = int(A.rows());
  VT lo = VT::Constant(n, 1e300), hi = VT::Constant(n, -1e300);
  std::vector<int> idx(n);
  std::function<void(int, int)> rec = [&](int start, int depth) {
    if (depth == n)
    {
      MT M(n, n);
      VT r(n);
      for (int i = 0; i < n; ++i) M.row(i) = A.row(idx[i]), r(i) = b(idx[i]);
      Eigen::FullPivLU<MT> lu(M);
      if (lu.rank() < n) return;
      VT x = lu.solve(r);
      if (((A * x - b).array() <= 1e-9 * (1 + b.cwiseAbs().maxCoeff())).all())
        lo = lo.cwiseMin(x), hi = hi.cwiseMax(x);
      return;
    }
    for (int i = start; i < m; ++i) idx[depth] = i, rec(i + 1, depth + 1);
  };
  rec(0, 0);
  return {lo, hi};
}

// Exact uniform sampler for the body it is handed (the inflated body inside
// sampleLattice), by rejection from its bounding box.  Records the volume of
// the body it was asked for (estimated from the acceptance rate).
struct ExactUniform
{
  std::mt19937_64 g;
  long long calls = 0, proposals = 0;
  double bodyVolume = 0;  // MC estimate of the requested body's volume
  explicit ExactUniform(unsigned seed) : g(seed) {}
  std::vector<VT> operator()(const MT& A, const VT& b, long long k, double& t)
  {
    t = 0;
    ++calls;
    auto [lo, hi] = bbox(A, b);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<VT> out;
    long long tries = 0;
    while ((long long)out.size() < k)
    {
      VT x = lo + (hi - lo).cwiseProduct(VT::NullaryExpr(lo.size(), [&] { return U(g); }));
      ++tries;
      if (((A * x - b).array() <= 0).all()) out.push_back(x);
    }
    proposals += k;
    if (calls == 1 && k > 0) bodyVolume = (hi - lo).prod() * double(k) / double(tries);
    return out;
  }
};

// Lattice points of P at precision p (n <= 3), via the bounding box.
std::vector<std::vector<long long>> latticePoints(const IBody& P, int p)
{
  auto [lo, hi] = bbox(P.Ad(), P.bd());
  const double s = std::pow(10.0, p);
  const int n = P.dim();
  std::vector<long long> klo(n), khi(n);
  for (int j = 0; j < n; ++j)
    klo[j] = (long long)std::floor(lo(j) * s) - 1, khi[j] = (long long)std::ceil(hi(j) * s) + 1;
  std::vector<std::vector<long long>> out;
  std::vector<long long> k(n);
  std::function<void(int)> rec = [&](int j) {
    if (j == n)
    {
      if (memberExact(P, k, p)) out.push_back(k);
      return;
    }
    for (k[j] = klo[j]; k[j] <= khi[j]; ++k[j]) rec(j + 1);
  };
  rec(0);
  return out;
}

std::vector<long long> toLattice(const VT& z, int p)
{
  std::vector<long long> k(z.size());
  for (Eigen::Index j = 0; j < z.size(); ++j) k[j] = std::llround(z(j) * std::pow(10.0, p));
  return k;
}

struct Uniformity
{
  int cells = 0, outside = 0, offLattice = 0, missing = 0;
  double chi2 = 0, pValue = 0;
  double zBoundary = 0, zCorner = 0, zFacet = 0;  // binomial z of boundary/corner/on-facet
  long long boundaryCount = 0, cornerCount = 0, facetCount = 0, N = 0;
  double boundaryExpect = 0, cornerExpect = 0, facetExpect = 0;
};

// Chi-square goodness of fit of samples to the uniform law on P cap 10^-p Z^n.
// A lattice point is "boundary" if its rounding cell [z-h, z+h]^n meets a
// facet (slack < h |A_i|_1), "corner" if it meets two. "On-facet": slack exactly 0 on some facet.
Uniformity uniformity(const IBody& P, int p, const std::vector<VT>& samples)
{
  Uniformity U;
  auto pts = latticePoints(P, p);
  std::map<std::vector<long long>, long long> count;
  for (auto& k : pts) count[k] = 0;
  for (const VT& z : samples)
  {
    auto k = toLattice(z, p);
    if ((z * std::pow(10.0, p) - Eigen::Map<Eigen::Matrix<long long, -1, 1>>(k.data(), k.size())
                                       .cast<double>())
            .cwiseAbs()
            .maxCoeff()
        > 1e-6)
      ++U.offLattice;
    auto it = count.find(k);
    if (it == count.end()) { ++U.outside; continue; }
    ++it->second;
  }
  U.cells = int(pts.size());
  U.N = (long long)samples.size();
  const double e = double(U.N) / U.cells;
  long long nb = 0, nc = 0, nf = 0;
  for (auto& [k, c] : count)
  {
    U.chi2 += (c - e) * (c - e) / e;
    if (c == 0 && U.missing < 6)
    {
      std::printf("  missing lattice point (");
      for (auto v : k) std::printf(" %lld", v);
      std::printf(" ) slacks:");
      for (size_t i = 0; i < P.A.size(); ++i) std::printf(" %lld", slack(P, i, k, p));
      std::printf("\n");
    }
    U.missing += (c == 0);
    int touching = 0;
    bool onFacet = false;
    for (size_t i = 0; i < P.A.size(); ++i)
    {
      onFacet |= slack(P, i, k, p) == 0;
      long long l1 = 0;
      for (auto a : P.A[i]) l1 += std::llabs(a);
      // slack (units 10^-q) < h |A_i|_1  (h = 10^{q-p} units)
      if (slack(P, i, k, p) < l1 * ipow10(P.q - p)) ++touching;
    }
    if (touching >= 1) ++nb, U.boundaryCount += c;
    if (touching >= 2) ++nc, U.cornerCount += c;
    if (onFacet) ++nf, U.facetCount += c;
  }
  boost::math::chi_squared dist(std::max(1, U.cells - 1));
  U.pValue = boost::math::cdf(boost::math::complement(dist, U.chi2));
  auto zgroup = [&](long long obs, long long cells) {
    const double pr = double(cells) / U.cells;
    if (cells == 0 || pr >= 1) return 0.0;
    return (obs - U.N * pr) / std::sqrt(U.N * pr * (1 - pr));
  };
  U.boundaryExpect = U.N * double(nb) / U.cells;
  U.cornerExpect = U.N * double(nc) / U.cells;
  U.zBoundary = zgroup(U.boundaryCount, nb);
  U.zCorner = zgroup(U.cornerCount, nc);
  U.facetExpect = U.N * double(nf) / U.cells;
  U.zFacet = zgroup(U.facetCount, nf);
  return U;
}

void print(const char* tag, const IBody& P, int p, const Uniformity& U)
{
  std::printf("%-28s %-22s b=%d cells=%5d N=%7lld chi2/df=%.3f p=%.3g missing=%d "
              "boundary %lld/%.0f (z=%+.2f) corner %lld/%.0f (z=%+.2f) on-facet %lld/%.0f (z=%+.2f) "
              "outside=%d offLattice=%d\n",
              tag, P.name.c_str(), p, U.cells, U.N, U.chi2 / std::max(1, U.cells - 1), U.pValue,
              U.missing, U.boundaryCount, U.boundaryExpect, U.zBoundary, U.cornerCount,
              U.cornerExpect, U.zCorner, U.facetCount, U.facetExpect, U.zFacet, U.outside,
              U.offLattice);
}

void checkUniform(const Uniformity& U)
{
  CHECK(U.outside == 0);
  CHECK(U.offLattice == 0);
  CHECK(U.pValue > 1e-3);           // chi-square, alpha = 0.001
  CHECK(std::abs(U.zBoundary) < 4);  // two-sided alpha ~ 6e-5
  CHECK(std::abs(U.zCorner) < 4);
  CHECK(std::abs(U.zFacet) < 4);  // lattice points lying exactly on a facet
}

constexpr double kEps = 0.8 / 12.0;  // what ttc passes (mvcEps)
constexpr double kDelta = 0.05;
// For the distribution tests the retry budget is made irrelevant (it is sized
// for the eps/3 rejection premise, which only holds for small h, see [budget]);
// 1e-300 gives maxBatches = 185, enough for 80% rejection.
constexpr double kNoBudget = 1e-300;

// Test bodies (2D, non-lattice and lattice-aligned facets).
IBody squareOffLattice() { return ibox({53, 21}, {917, 978}, 3, "square[.053,.917]x[.021,.978]"); }
IBody rect3() { return ibox({53, 21}, {2917, 978}, 3, "rect[.053,2.917]x[.021,.978]"); }
IBody triangle()
{
  // x >= 0.013, y >= 0.027, x + y <= 0.961
  IBody P;
  P.name = "triangle(.013,.027,.961)";
  P.q = 3;
  P.A = {{-1, 0}, {0, -1}, {1, 1}};
  P.B = {-13, -27, 961};
  return P;
}
IBody latticeTriangle()
{
  // x >= 0, y >= 0, 3x + 7y <= 2.1 : facets through lattice points at b=1,2
  IBody P;
  P.name = "triangle 3x+7y<=2.1";
  P.q = 2;
  P.A = {{-1, 0}, {0, -1}, {3, 7}};
  P.B = {0, 0, 210};
  return P;
}
IBody thinSlab()
{
  // 0.307 <= x + 2y <= 0.352, 0 <= x <= 0.3  (width 0.045/sqrt5 ~ 0.02)
  IBody P;
  P.name = "slab .307<=x+2y<=.352";
  P.q = 3;
  P.A = {{1, 2}, {-1, -2}, {1, 0}, {-1, 0}};
  P.B = {352, -307, 300, 0};
  return P;
}

}  // namespace

// ------------------------------------------------------------ (1) rounding
TEST_CASE("lattice [rounding]: fixed proposal -> two neighbours, unbiased")
{
  // P is a large box, so no rounded point is rejected; generate always
  // returns the same point x.
  Body B = box({-50, -50, -50}, {50, 50, 50});
  const int p = 2;
  const long long N = 200000;
  for (VT x : {VT((VT(3) << 0.123456, -7.891011, 0.5).finished()),
               VT((VT(3) << 12.345, -0.005, 0.25).finished()),
               VT((VT(3) << 0.37, -1.0, 49.999999).finished())})
  {
    auto gen = [&](const MT&, const VT&, long long k, double& t) {
      t = 0;
      return std::vector<VT>(size_t(k), x);
    };
    double time = 0;
    auto out = ttc::sampleLattice<double>(B.A, B.b, N, p, kEps, kDelta, 12345u, gen, time);
    REQUIRE((long long)out.size() == N);
    for (int j = 0; j < 3; ++j)
    {
      const double s = 100.0, raw = x(j) * s, lo = std::floor(raw);
      const double frac = raw - lo;
      long long nlo = 0, nhi = 0, other = 0;
      double mean = 0;
      for (const VT& z : out)
      {
        const long long k = std::llround(z(j) * s);
        if (k == (long long)lo) ++nlo;
        else if (k == (long long)lo + 1) ++nhi;
        else ++other;
        mean += z(j);
      }
      mean /= N;
      // E[z] = x; Var[z] = h^2 f (1-f)
      const double sd = 0.01 * std::sqrt(frac * (1 - frac) / N);
      const double zstat = sd > 0 ? (mean - x(j)) / sd : (mean - x(j)) / 1e-15;
      std::printf("rounding x=%.8g: lo=%lld hi=%lld other=%lld mean-x=%.3g z=%+.2f\n", x(j),
                  nlo, nhi, other, mean - x(j), zstat);
      CHECK(other == 0);
      // x on the lattice up to double rounding of x*10^b: always its own point
      if (frac < 1e-9) CHECK(nhi == 0);
      else if (frac > 1 - 1e-9) CHECK(nlo == 0);
      else CHECK(std::abs(zstat) < 4.5);
    }
  }
}

TEST_CASE("lattice [rounding]: exact lattice proposals are returned unchanged")
{
  Body B = box({-10, -10}, {10, 10});
  for (int p : {0, 1, 2, 3})
  {
    const double h = std::pow(10.0, -p);
    VT x(2);
    x << 3 * h, -7 * h;
    auto gen = [&](const MT&, const VT&, long long k, double& t) {
      t = 0;
      return std::vector<VT>(size_t(k), x);
    };
    double time = 0;
    auto out = ttc::sampleLattice<double>(B.A, B.b, 2000, p, kEps, kDelta, 7u, gen, time);
    long long moved = 0;
    for (const VT& z : out) moved += (toLattice(z, p) != toLattice(x, p));
    INFO("p=" << p);
    CHECK(moved == 0);
  }
}

// --------------------------------------------- (2) exact uniform proposals
TEST_CASE("lattice [exact]: exact-uniform proposals give the uniform law on P cap hZ^n")
{
  struct Case { IBody P; int p; long long N; };
  std::vector<Case> cases = {
      {squareOffLattice(), 1, 20000},
      {squareOffLattice(), 2, 200000},
      {triangle(), 1, 20000},
      {triangle(), 2, 200000},
      {latticeTriangle(), 1, 20000},   // facets pass through lattice points
      {latticeTriangle(), 2, 200000},
      {thinSlab(), 2, 50000},
      {thinSlab(), 3, 200000},
      {ibox({53, 21, 107}, {917, 978, 893}, 3, "box3[.053,.917]x..."), 1, 50000},
  };
  for (auto& c : cases)
  {
    ExactUniform gen(1000u + c.p);
    double time = 0;
    auto out = ttc::sampleLattice<double>(c.P.Ad(), c.P.bd(), c.N, c.p, kEps, kNoBudget, 99u,
                                          std::ref(gen), time);
    REQUIRE((long long)out.size() == c.N);
    auto U = uniformity(c.P, c.p, out);
    print("exact", c.P, c.p, U);
    INFO(c.P.name << " b=" << c.p);
    checkUniform(U);
  }
}

// ----------------------------------------- (3) inflation, rejection, budget
namespace
{
// GetPrecision of the revised manuscript, as ttc computes it
// (lra_volume.cpp precisionFromCubes): eta = eps/96.
int ttcPrecision(int dim, int facets, int cubes, double minRadius, double epsilon)
{
  const double eta = epsilon / 96.0;
  const double logGamma = std::log(64.0) + std::log(double(cubes)) + std::log(double(dim))
                          - std::log(eta)
                          + 0.5 * std::log(std::log(4.0) + std::log(double(facets)) - std::log(eta));
  return int(std::max(0.0, std::ceil((logGamma - std::log(minRadius)) / std::log(10.0))));
}

// Measured rejection rate (budget disabled) and theory 1 - |P cap hZ^n| h^n / vol(P').
struct Rejection { double measured, theory; long long batches; };
Rejection rejection(const IBody& P, int p, long long N, bool lattice)
{
  ExactUniform gen(5u);
  double time = 0;
  auto out = ttc::sampleLattice<double>(P.Ad(), P.bd(), N, p, kEps, kNoBudget, 3u,
                                        std::ref(gen), time);
  Rejection R;
  R.measured = 1.0 - double(out.size()) / double(gen.proposals);
  R.theory = lattice ? 1.0 - double(latticePoints(P, p).size()) * std::pow(10.0, -p * P.dim())
                                 / gen.bodyVolume
                     : NAN;
  R.batches = gen.calls;
  return R;
}
}  // namespace

TEST_CASE("lattice [budget]: kappa covers the rounding kernel; rejection rate vs eps/3")
{
  // Exactness: the rounding kernel of lattice point z is the tent on
  // [z-h, z+h]^n.  It lies inside the inflated body iff kappa |A_i|_2 >= |A_i|_1
  // for every facet (worst case |A|_1/|A|_2 = sqrt n).  Beyond that the paper's
  // Hoeffding argument bounds the lost mass per facet by exp(-kappa^2/2).
  for (int m : {3, 4, 6, 8})
  {
    const double kappa = std::sqrt(2 * std::log(16 / kEps)) + std::sqrt(2 * std::log(double(m)));
    std::printf("kappa(eps=%.4f, m=%d) = %.3f  -> kernel inside inflated body for n <= %d\n",
                kEps, m, kappa, int(kappa * kappa));
    CHECK(kappa * kappa >= 2);  // all 2D bodies tested here are exact
  }

  // (a) Small b: rejection far above eps/3 (diagnostic table, matches theory).
  struct Case { IBody P; int p; };
  for (auto& c : std::vector<Case>{{squareOffLattice(), 1}, {squareOffLattice(), 2},
                                   {squareOffLattice(), 3}, {triangle(), 2}, {triangle(), 3},
                                   {thinSlab(), 2}, {thinSlab(), 3}})
  {
    auto R = rejection(c.P, c.p, 20000, true);
    std::printf("rejection %-30s b=%d measured=%.4f theory=%.4f eps/3=%.4f batches=%lld\n",
                c.P.name.c_str(), c.p, R.measured, R.theory, kEps / 3, R.batches);
    INFO(c.P.name << " b=" << c.p);
    // binomial: N=20000 proposals, sd <= 0.0036
    CHECK(std::abs(R.measured - R.theory) < 4 * std::sqrt(0.25 / 20000) + 2e-3);
  }

  // (b) At ttc's own precision (GetPrecision with the body's inner radius) the
  // premise holds: rejection <= eps/3 and the default budget never throws.
  for (const IBody& P : {squareOffLattice(), rect3(), triangle(), latticeTriangle(), thinSlab()})
  {
    const double r = ttc::certifiedInnerBall(P.Ad(), P.bd()).radius;
    const int p = ttcPrecision(P.dim(), int(P.A.size()), 1, r, 0.8);
    auto R = rejection(P, p, 20000, false);
    std::printf("ttc precision %-30s r=%.4f b=%d rejection=%.2e (eps/3=%.4f)\n", P.name.c_str(),
                r, p, R.measured, kEps / 3);
    INFO(P.name << " b=" << p);
    CHECK(R.measured <= kEps / 3);
    int threw = 0;
    for (unsigned s = 0; s < 200; ++s)
    {
      ExactUniform gen(s);
      double t = 0;
      try
      {
        auto out = ttc::sampleLattice<double>(P.Ad(), P.bd(), 1 + s % 50, p, kEps, kDelta, s,
                                              std::ref(gen), t);
        CHECK(out.size() == 1 + s % 50);
      }
      catch (const std::runtime_error&) { ++threw; }
    }
    CHECK(threw == 0);
  }
}

TEST_CASE("lattice [budget]: n = 0, n = 1, and bodies of width ~ h")
{
  IBody P = squareOffLattice();
  {
    ExactUniform gen(1u);
    double t = 0;
    auto out = ttc::sampleLattice<double>(P.Ad(), P.bd(), 0, 2, kEps, kDelta, 1u, std::ref(gen), t);
    CHECK(out.empty());
    CHECK(gen.calls == 0);
  }
  {
    // n = 1, many seeds: the single point is uniform over the lattice points.
    std::vector<VT> all;
    for (unsigned s = 0; s < 20000; ++s)
    {
      ExactUniform gen(s + 77777u);
      double t = 0;
      auto out =
          ttc::sampleLattice<double>(P.Ad(), P.bd(), 1, 1, kEps, kNoBudget, s, std::ref(gen), t);
      REQUIRE(out.size() == 1);
      all.push_back(out[0]);
    }
    auto U = uniformity(P, 1, all);
    print("exact n=1 x20000 seeds", P, 1, U);
    checkUniform(U);
  }
  {
    // With the default budget and b too coarse for the body (premise violated):
    // the sampler must throw rather than return fewer points.  n=1 gets
    // maxBatches = ceil(ln(delta)/ln(eps/3)) = 1 batch, so it throws with
    // probability = rejection rate.
    int threw = 0, shortRuns = 0;
    for (unsigned s = 0; s < 500; ++s)
    {
      ExactUniform gen(s);
      double t = 0;
      try
      {
        auto out = ttc::sampleLattice<double>(P.Ad(), P.bd(), 1, 2, kEps, kDelta, s,
                                              std::ref(gen), t);
        shortRuns += out.size() != 1;
      }
      catch (const std::runtime_error&) { ++threw; }
    }
    std::printf("square b=2 n=1 default budget: threw %d/500 (rejection ~0.20)\n", threw);
    CHECK(shortRuns == 0);
  }
  {
    // Strip of width ~2h (y in {0.01, 0.02} at b=2), n = 1000.
    IBody S = ibox({0, 1}, {1000, 22}, 3, "strip[0,1]x[.001,.022]");
    int threw = 0, shortRuns = 0;
    for (unsigned s = 0; s < 20; ++s)
    {
      ExactUniform gen(s);
      double t = 0;
      try
      {
        auto out = ttc::sampleLattice<double>(S.Ad(), S.bd(), 1000, 2, kEps, kDelta, s,
                                              std::ref(gen), t);
        shortRuns += out.size() != 1000u;
      }
      catch (const std::runtime_error&) { ++threw; }
    }
    auto R = rejection(S, 2, 20000, true);
    std::printf("strip width ~2h (b=2): rejection %.3f (theory %.3f), threw %d/20 runs\n",
                R.measured, R.theory, threw);
    CHECK(shortRuns == 0);
    // Uniform over the 101 x 2 lattice points when the budget is lifted:
    ExactUniform gen(9u);
    double t = 0;
    auto out = ttc::sampleLattice<double>(S.Ad(), S.bd(), 20000, 2, kEps, kNoBudget, 9u,
                                          std::ref(gen), t);
    auto U = uniformity(S, 2, out);
    print("exact strip", S, 2, U);
    checkUniform(U);
  }
}
// --------------------------------------------- (4) composed with ttc walks
namespace
{
const char* walkName(ttc::SamplerWalk w)
{
  switch (w)
  {
    case ttc::SamplerWalk::Billiard: return "Billiard";
    case ttc::SamplerWalk::AcceleratedBilliard: return "AccelBilliard";
    case ttc::SamplerWalk::Ball: return "Ball";
    case ttc::SamplerWalk::RDHR: return "RDHR";
    case ttc::SamplerWalk::CDHR: return "CDHR";
  }
  return "?";
}

void composed(const IBody& P, int p, long long N, ttc::SamplerWalk w, unsigned seed)
{
  auto out = ttcLatticeSamples(P.Ad(), P.bd(), N, p, seed, w, 0.8, kNoBudget);
  REQUIRE((long long)out.size() == N);
  auto U = uniformity(P, p, out);
  char tag[64];
  std::snprintf(tag, sizeof tag, "walk=%s", walkName(w));
  print(tag, P, p, U);
  INFO(P.name << " b=" << p << " walk=" << std::string(walkName(w)));
  checkUniform(U);
}
}  // namespace

TEST_CASE("lattice [walks]: ttc's default sampler (billiard) is uniform on lattice points")
{
  composed(squareOffLattice(), 1, 20000, ttc::SamplerWalk::Billiard, 11u);
  composed(rect3(), 1, 30000, ttc::SamplerWalk::Billiard, 12u);
  composed(triangle(), 1, 20000, ttc::SamplerWalk::Billiard, 13u);
  composed(thinSlab(), 2, 20000, ttc::SamplerWalk::Billiard, 14u);
}

TEST_CASE("lattice [walks]: long box [0,2]x[0,1000] at b=0, y-deciles" * doctest::test_suite("slow"))
{
  IBody P = ibox({0, 0}, {2, 1000}, 0, "box[0,2]x[0,1000]");
  const long long N = 10000;
  auto out = ttcLatticeSamples(P.Ad(), P.bd(), N, 0, 5u, ttc::SamplerWalk::Billiard, 0.8, kNoBudget);
  // Lattice points: x in {0,1,2}, y in {0..1000}; decile d holds y in
  // [100d, 100d+99] (100 values), decile 9 holds 101 values.
  std::vector<long long> c(10, 0);
  for (auto& z : out) c[std::min(9, int(std::llround(z(1)) / 100))]++;
  double chi2 = 0;
  for (int d = 0; d < 10; ++d)
  {
    const double e = N * (d == 9 ? 101.0 : 100.0) / 1001.0;
    chi2 += (c[d] - e) * (c[d] - e) / e;
  }
  boost::math::chi_squared dist(9);
  const double pv = boost::math::cdf(boost::math::complement(dist, chi2));
  std::printf("long box b=0 billiard y-deciles:");
  for (auto v : c) std::printf(" %lld", v);
  std::printf("  chi2=%.1f p=%.3g\n", chi2, pv);
  CHECK(pv > 1e-3);
}

TEST_CASE("lattice [walks]: every ttc walk, square/rect/triangle/slab" * doctest::test_suite("slow"))
{
  for (auto w : {ttc::SamplerWalk::AcceleratedBilliard, ttc::SamplerWalk::Ball,
                 ttc::SamplerWalk::RDHR, ttc::SamplerWalk::CDHR})
  {
    composed(squareOffLattice(), 1, 20000, w, 21u);
    composed(rect3(), 1, 30000, w, 22u);
    composed(triangle(), 1, 20000, w, 23u);
    composed(thinSlab(), 2, 20000, w, 24u);
  }
}

// ------------------------------------------------------------- (5) seeding
TEST_CASE("lattice [seeding]: rounding seed equal to the walk seed")
{
  // ttc derives the rounding seed and every walk seed from distinct splitmix64
  // stream draws (lra_volume.cpp streamSeed); the test helper uses seed and
  // seed*2654435761+1.  Worst case: identical seeds, so std::mt19937 (rounding)
  // and boost::mt19937 (walk) emit the same raw stream.
  IBody P = squareOffLattice();
  const int p = 1;
  const long long N = 20000;
  const unsigned seed = 4242u;
  // First batch: walk seed == rounding seed.  Later batches advance the walk
  // seed (a constant one would replay the same rejected proposals forever).
  unsigned walkSeed = seed;
  auto gen = [&](const MT& A, const VT& b, long long k, double& t) {
    double bs = 0;
    t = 0;
    return ttc::sampleRestartWalk(A, b, k, ttc::SamplerWalk::Billiard,
                                  ttcWalkLength(int(A.cols())), walkSeed++, bs);
  };
  double time = 0;
  auto out = ttc::sampleLattice<double>(P.Ad(), P.bd(), N, p, kEps, kNoBudget, seed, gen, time);
  auto U = uniformity(P, p, out);
  print("same-seed billiard", P, p, U);
  checkUniform(U);
}
