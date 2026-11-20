// T3 of text-for-agents/ttc-union-test-plan.md: GetPrecision (KR25 manuscript,
// Alg. 2) and the lattice-count lemmas lem:nonconv / lem:lattice
// (sections/analysis.tex), against ttc::latticePrecision
// (src/volume/get_precision.cpp) and the (m, r, gamma_min) that
// getPrecision passes to it.
//
//   [alg2]     latticePrecision(n, r, m, gamma_min, eps/96) == ceil(log10(max(1, gamma/gamma_min))),
//              gamma = 64 m n / eta * sqrt(ln(4 r / eta)), eta = eps/96, recomputed in
//              long double on a grid.  The log must be natural: the Azuma step of
//              lem:nonconv (claim cl, part a) needs 2 exp(-kappa^2/2) <= eta/(4r) with
//              kappa = sqrt(2 ln(8/eta)) + sqrt(2 ln r).
//   [proof]    the step "(1 + kappa/gamma)^n <= (1 + eta/(32 m n))^n" of lem:nonconv, i.e.
//              kappa <= 2 sqrt(ln(4r/eta)), for Alg. 2's gamma (expected to FAIL: the
//              cited inequality sqrt u + sqrt v <= sqrt(2(u+v)) only gives
//              kappa <= 2 sqrt(ln(8r/eta))), and the lemma's conclusion
//              (1 - eta)|Q cap Z^n| <= vol Q <= (1 + eta)|Q cap Z^n| re-derived from the
//              proof's own chain with the actual kappa/gamma (must hold).
//   [counter]  the exact lattice counter used below agrees with closed forms (box
//              inclusion-exclusion, Ehrhart polynomials of simplices).
//   [lattice]  lem:lattice on small 2D/3D unions (boxes, L-shapes, frame, staircase,
//              skinny L, triangles tiling a square / a hexagon, overlapping triangles,
//              thin triangles, simplex, cube u simplex): with b from latticePrecision
//              and ttc's certified inner balls, |10^{-bn} |Q cap L| / vol Q - 1| <= eta.
//              Also reports the smallest b* that already gives eta-accuracy and what the
//              extra b_ttc - b* digits cost (C = Thresh 10^{bn}, halvings, GMP digits).
//   [cost]     GMP membership (SampleMultiset<MpFloat>::removeInside) time at b + 32 digits.
#include "doctest.h"
#include "common.hpp"

#include "volume/gmp_walks.hpp"
#include "volume/get_precision.hpp"
#include "volume/sample_multiset.hpp"
#include "volume/union_parameters.hpp"

#include <chrono>
#include <cstdio>
#include <functional>
#include <numeric>

using namespace vt;

namespace
{

using i64 = long long;
using i128 = __int128;
constexpr i64 kMilli = 1000;  // right-hand sides and vertices are integers / 1000

i64 pow10ll(int e)
{
  i64 r = 1;
  while (e-- > 0) r *= 10;
  return r;
}
i64 floorDiv(i64 a, i64 b)
{
  i64 q = a / b;
  if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
  return q;
}
i64 ceilDiv(i64 a, i64 b) { return -floorDiv(-a, b); }

// Polytope { x : a_j . x <= B_j / 1000 } with integer a_j, B_j; bounding box
// lo/hi in milli units.  Membership of a lattice point x = k 10^-b is exact:
// 1000 a_j . k <= B_j 10^b.
struct IPoly
{
  std::vector<std::vector<i64>> a;
  std::vector<i64> B;
  std::vector<i64> lo, hi;
  bool isBox = false;
};

IPoly boxM(const std::vector<i64>& lo, const std::vector<i64>& hi)
{
  const int n = int(lo.size());
  IPoly P;
  P.lo = lo;
  P.hi = hi;
  P.isBox = true;
  for (int i = 0; i < n; ++i)
  {
    std::vector<i64> e(n, 0);
    e[i] = 1;
    P.a.push_back(e);
    P.B.push_back(hi[i]);
    e[i] = -1;
    P.a.push_back(e);
    P.B.push_back(-lo[i]);
  }
  return P;
}

// Convex polygon from its vertices (milli units, either orientation).
IPoly polygonM(std::vector<std::array<i64, 2>> v)
{
  i64 area2 = 0;
  for (std::size_t i = 0; i < v.size(); ++i)
  {
    const auto& p = v[i];
    const auto& q = v[(i + 1) % v.size()];
    area2 += p[0] * q[1] - q[0] * p[1];
  }
  if (area2 < 0) std::reverse(v.begin(), v.end());
  IPoly P;
  P.lo = {v[0][0], v[0][1]};
  P.hi = P.lo;
  for (std::size_t i = 0; i < v.size(); ++i)
  {
    const auto& p = v[i];
    const auto& q = v[(i + 1) % v.size()];
    i64 n0 = q[1] - p[1], n1 = -(q[0] - p[0]);  // outward normal of a CCW edge
    const i64 g = std::gcd(std::abs(n0), std::abs(n1));
    n0 /= g;
    n1 /= g;
    P.a.push_back({n0, n1});
    P.B.push_back(n0 * p[0] + n1 * p[1]);
    for (int d = 0; d < 2; ++d)
    {
      P.lo[d] = std::min(P.lo[d], p[d]);
      P.hi[d] = std::max(P.hi[d], p[d]);
    }
  }
  return P;
}

// {x >= 0, sum x <= s}.
IPoly simplexM(int n, i64 s)
{
  IPoly P;
  for (int i = 0; i < n; ++i)
  {
    std::vector<i64> e(n, 0);
    e[i] = -1;
    P.a.push_back(e);
    P.B.push_back(0);
  }
  P.a.push_back(std::vector<i64>(n, 1));
  P.B.push_back(s);
  P.lo.assign(n, 0);
  P.hi.assign(n, s);
  return P;
}

struct UnionCase
{
  std::string name;
  int n;
  std::vector<IPoly> polys;
  long double vol;
  std::function<i128(int)> formula;  // exact |Q cap L| at precision b, if known
};

i128 binom(i128 N, int k)
{
  i128 r = 1;
  for (int i = 1; i <= k; ++i) r = r * (N - k + i) / i;
  return r;
}

std::vector<UnionCase> cases()
{
  std::vector<UnionCase> C;
  auto B2 = [](i64 x0, i64 y0, i64 x1, i64 y1) { return boxM({x0, y0}, {x1, y1}); };
  const i64 K = kMilli;
  C.push_back({"square", 2, {B2(0, 0, K, K)}, 1.0L, {}});
  C.push_back({"L-shape", 2, {B2(0, 0, 2 * K, K), B2(0, 0, K, 2 * K)}, 3.0L, {}});
  C.push_back({"frame", 2,
               {B2(0, 0, 3 * K, K), B2(0, 2 * K, 3 * K, 3 * K), B2(0, 0, K, 3 * K),
                B2(2 * K, 0, 3 * K, 3 * K)},
               8.0L, {}});
  C.push_back({"staircase", 2, {B2(0, 0, K, 3 * K), B2(0, 0, 2 * K, 2 * K), B2(0, 0, 3 * K, K)},
               6.0L, {}});
  C.push_back({"skinny-L", 2, {B2(0, 0, K, 10), B2(0, 0, 10, K)}, 0.0199L, {}});
  C.push_back({"2 triangles = square", 2,
               {polygonM({{{0, 0}, {K, 0}, {0, K}}}), polygonM({{{K, 0}, {K, K}, {0, K}}})},
               1.0L, {}});
  {
    const std::array<i64, 2> c{1500, 1000};
    const std::vector<std::array<i64, 2>> hv = {{1000, 0}, {2000, 0}, {3000, 1000},
                                                {2000, 2000}, {1000, 2000}, {0, 1000}};
    UnionCase H{"hexagon = 6 triangles", 2, {}, 4.0L, {}};
    for (std::size_t i = 0; i < hv.size(); ++i)
      H.polys.push_back(polygonM({c, hv[i], hv[(i + 1) % hv.size()]}));
    C.push_back(H);
  }
  C.push_back({"W = 2 overlapping triangles", 2,
               {polygonM({{{0, 0}, {2 * K, 0}, {0, 2 * K}}}),
                polygonM({{{K, 0}, {3 * K, 0}, {K, 2 * K}}})},
               3.5L, {}});
  C.push_back({"thin rectangle = 2 triangles", 2,
               {polygonM({{{0, 0}, {K, 0}, {0, 10}}}), polygonM({{{K, 0}, {K, 10}, {0, 10}}})},
               0.01L, {}});
  // 3D.
  C.push_back({"cube3", 3, {boxM({0, 0, 0}, {K, K, K})}, 1.0L, {}});
  C.push_back({"3D L", 3, {boxM({0, 0, 0}, {2 * K, K, K}), boxM({0, 0, 0}, {K, 2 * K, K})}, 3.0L,
               {}});
  C.push_back({"simplex3", 3, {simplexM(3, K)}, 1.0L / 6.0L,
               [](int b) { return binom(i128(pow10ll(b)) + 3, 3); }});
  // |cube u {x >= 0, sum x <= 2}| = C(2t+3, 3) + C(t+2, 3), t = 10^b; vol 4/3 + 1 - 5/6.
  C.push_back({"cube3 u simplex(2)", 3, {boxM({0, 0, 0}, {K, K, K}), simplexM(3, 2 * K)}, 1.5L,
               [](int b) {
                 const i128 t = pow10ll(b);
                 return binom(2 * t + 3, 3) + binom(t + 2, 3);
               }});
  return C;
}

// Inclusion-exclusion over boxes (exact for any b).
i128 countBoxes(const std::vector<IPoly>& polys, int n, int b)
{
  const i64 s = pow10ll(b);
  const int m = int(polys.size());
  i128 total = 0;
  for (int mask = 1; mask < (1 << m); ++mask)
  {
    std::vector<i64> lo(n, std::numeric_limits<i64>::min()), hi(n, std::numeric_limits<i64>::max());
    for (int p = 0; p < m; ++p)
      if (mask >> p & 1)
        for (int i = 0; i < n; ++i)
        {
          lo[i] = std::max(lo[i], polys[p].lo[i]);
          hi[i] = std::min(hi[i], polys[p].hi[i]);
        }
    i128 c = 1;
    for (int i = 0; i < n; ++i)
    {
      const i64 kl = ceilDiv(lo[i] * s, kMilli), kh = floorDiv(hi[i] * s, kMilli);
      c *= kh >= kl ? i128(kh - kl + 1) : i128(0);
    }
    total += (__builtin_popcount(unsigned(mask)) & 1) ? c : -c;
  }
  return total;
}

// Column enumeration: iterate over all coordinates but the longest one, take
// the union of the polytopes' integer intervals on it.  Returns false if the
// work (columns x polytopes) exceeds `budget`.
bool enumerate(const UnionCase& U, int b, double budget, i128& out)
{
  const int n = U.n;
  const i64 s = pow10ll(b);
  std::vector<i64> klo(n, std::numeric_limits<i64>::max()), khi(n, std::numeric_limits<i64>::min());
  for (const auto& P : U.polys)
    for (int i = 0; i < n; ++i)
    {
      klo[i] = std::min(klo[i], ceilDiv(P.lo[i] * s, kMilli));
      khi[i] = std::max(khi[i], floorDiv(P.hi[i] * s, kMilli));
    }
  int last = 0;
  for (int i = 1; i < n; ++i)
    if (khi[i] - klo[i] > khi[last] - klo[last]) last = i;
  double work = double(U.polys.size());
  for (int i = 0; i < n; ++i)
    if (i != last) work *= double(khi[i] - klo[i] + 1);
  if (work > budget) return false;
  // Overflow guard for 1000 a.k and B 10^b in 64 bits.
  double kAbs = 0;
  for (int i = 0; i < n; ++i)
    kAbs = std::max({kAbs, std::abs(double(klo[i])), std::abs(double(khi[i]))});
  for (const auto& P : U.polys)
    for (std::size_t j = 0; j < P.a.size(); ++j)
    {
      i64 amax = 0;
      for (auto x : P.a[j]) amax = std::max(amax, std::abs(x));
      REQUIRE(double(amax) * kMilli * n * kAbs < 1e18);
      REQUIRE(std::abs(double(P.B[j])) * double(s) < 1e18);
    }
  std::vector<int> others;
  for (int i = 0; i < n; ++i)
    if (i != last) others.push_back(i);
  std::vector<i64> k(n, 0);
  std::vector<std::pair<i64, i64>> iv;
  i128 total = 0;
  std::function<void(std::size_t)> rec = [&](std::size_t d) {
    if (d < others.size())
    {
      const int c = others[d];
      for (i64 v = klo[c]; v <= khi[c]; ++v)
      {
        k[c] = v;
        rec(d + 1);
      }
      return;
    }
    iv.clear();
    for (const auto& P : U.polys)
    {
      bool inBox = true;
      for (int c : others)
        if (k[c] * kMilli < P.lo[c] * s || k[c] * kMilli > P.hi[c] * s)
        {
          inBox = false;
          break;
        }
      if (!inBox) continue;
      i64 L = std::numeric_limits<i64>::min(), H = std::numeric_limits<i64>::max();
      for (std::size_t j = 0; j < P.a.size() && L <= H; ++j)
      {
        i64 rhs = P.B[j] * s;
        for (int c : others) rhs -= kMilli * P.a[j][c] * k[c];
        const i64 alpha = kMilli * P.a[j][last];
        if (alpha > 0) H = std::min(H, floorDiv(rhs, alpha));
        else if (alpha < 0) L = std::max(L, ceilDiv(rhs, alpha));
        else if (rhs < 0) H = L - 1;
      }
      if (L <= H) iv.push_back({L, H});
    }
    if (iv.empty()) return;
    std::sort(iv.begin(), iv.end());
    i64 cl = iv[0].first, ch = iv[0].second;
    for (std::size_t t = 1; t < iv.size(); ++t)
    {
      if (iv[t].first <= ch + 1)
        ch = std::max(ch, iv[t].second);
      else
      {
        total += ch - cl + 1;
        cl = iv[t].first;
        ch = iv[t].second;
      }
    }
    total += ch - cl + 1;
  };
  rec(0);
  out = total;
  return true;
}

bool allBoxes(const UnionCase& U)
{
  return std::all_of(U.polys.begin(), U.polys.end(), [](const IPoly& P) { return P.isBox; });
}

bool exactCount(const UnionCase& U, int b, double budget, i128& out)
{
  if (U.formula)
  {
    out = U.formula(b);
    return true;
  }
  if (allBoxes(U))
  {
    out = countBoxes(U.polys, U.n, b);
    return true;
  }
  return enumerate(U, b, budget, out);
}

void toDouble(const IPoly& P, int n, Eigen::MatrixXd& A, Eigen::VectorXd& b)
{
  A.resize(Eigen::Index(P.a.size()), n);
  b.resize(Eigen::Index(P.a.size()));
  for (std::size_t j = 0; j < P.a.size(); ++j)
  {
    for (int i = 0; i < n; ++i) A(Eigen::Index(j), i) = double(P.a[j][i]);
    b(Eigen::Index(j)) = double(P.B[j]) / double(kMilli);
  }
}

// Alg. 2 by hand, natural log (long double).  kAlg2RFactor is the manuscript's
// factor in sqrt(ln(8r/eta)).  It was 4, but the [proof] tests showed
// lem:nonconv needs 8; the paper and get_precision.cpp were changed to 8.
constexpr long double kAlg2RFactor = 8.0L;
long double alg2Gamma(long double n, long double r, long double m, long double eps)
{
  const long double eta = eps / 96.0L;
  return 64.0L * m * n / eta * std::sqrt(std::log(kAlg2RFactor * r / eta));
}
long double alg2LogRatio(long double n, long double r, long double m, long double gHat,
                         long double eps)
{
  return std::max(0.0L, std::log10(alg2Gamma(n, r, m, eps) / gHat));
}

std::string i128str(i128 v)
{
  if (v == 0) return "0";
  std::string s;
  bool neg = v < 0;
  if (neg) v = -v;
  while (v > 0)
  {
    s += char('0' + int(v % 10));
    v /= 10;
  }
  if (neg) s += '-';
  return std::string(s.rbegin(), s.rend());
}

}  // namespace

// ============================================================================
TEST_CASE("[T3][alg2] latticePrecision equals Alg. 2 (natural log) on a grid")
{
  int checked = 0, bad = 0, skipped = 0, diffLog10 = 0, diffLog2 = 0;
  for (double n : {1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 40.0})
    for (double r : {2.0, 4.0, 6.0, 12.0, 100.0, 1e4, 1e6})
      for (double m : {1.0, 2.0, 3.0, 10.0, 1000.0, 1e5})
        for (double g : {1e-12, 1e-6, 1e-3, 0.01, 0.5, 1.0, 10.0, 1e3, 1e6, 1e9})
          for (double eps : {0.05, 0.3, 0.8, 1.0})
          {
            const long double lr = alg2LogRatio(n, r, m, g, eps);
            if (lr > 0 && std::abs(lr - std::round(lr)) < 1e-9)
            {
              ++skipped;  // ceil is ill-conditioned at an integer
              continue;
            }
            const int want = int(std::ceil(lr));
            const int got = ttc::latticePrecision(std::size_t(n), std::size_t(r),
                                                    std::size_t(m), g, eps / 96.0);
            ++checked;
            if (got != want)
            {
              ++bad;
              if (bad < 10)
                MESSAGE("n=" << n << " r=" << r << " m=" << m << " gHat=" << g << " eps=" << eps
                             << ": got " << got << ", Alg. 2 gives " << want);
            }
            // The same formula with log10 / log2 inside the sqrt.
            const long double eta = eps / 96.0L;
            const long double base = 64.0L * m * n / eta;
            auto bWith = [&](long double lg) {
              return int(std::ceil(std::max(0.0L, std::log10(base * std::sqrt(lg) / g))));
            };
            diffLog10 += bWith(std::log10(kAlg2RFactor * r / eta)) != want;
            diffLog2 += bWith(std::log2(kAlg2RFactor * r / eta)) != want;
          }
  MESSAGE(checked << " grid points (" << skipped << " at integer log10 ratios skipped), " << bad
                  << " mismatches; with log10 in the sqrt b would differ at " << diffLog10
                  << ", with log2 at " << diffLog2);
  CHECK(bad == 0);
  // b = 0 when the inner radius already exceeds gamma (Alg. 2 starts at maxratio = 1).
  CHECK(ttc::latticePrecision(2, 4, 1, 1e12, 0.8 / 96.0) == 0);

  // Natural log is what the Azuma step of claim cl needs:
  // 2 exp(-kappa^2/2) <= eta/(4r), kappa = sqrt(2 ln(8/eta)) + sqrt(2 ln r).
  int azumaBad = 0, azumaBad10 = 0;
  for (double r : {1.0, 2.0, 4.0, 12.0, 100.0, 1e4, 1e6, 1e9})
    for (double eps : {0.01, 0.05, 0.3, 0.8, 1.0})
    {
      const long double eta = eps / 96.0L;
      const long double kLn = std::sqrt(2 * std::log(8 / eta)) + std::sqrt(2 * std::log((long double)r));
      const long double k10 = std::sqrt(2 * std::log10(8 / eta)) + std::sqrt(2 * std::log10((long double)r));
      azumaBad += !(2 * std::exp(-kLn * kLn / 2) <= eta / (4 * r) * (1 + 1e-12L));
      azumaBad10 += !(2 * std::exp(-k10 * k10 / 2) <= eta / (4 * r));
    }
  MESSAGE("Azuma tail bound: fails at " << azumaBad << " points with ln, at " << azumaBad10
                                        << "/40 with log10");
  CHECK(azumaBad == 0);
}

TEST_CASE("[T3][proof] lem:nonconv step kappa/gamma <= eta/(32 m n) with Alg. 2's gamma")
{
  // Alg. 2's gamma = 64 m n/eta sqrt(ln(4r/eta)), so the step is
  // kappa <= 2 sqrt(ln(4r/eta)).  The proof cites sqrt u + sqrt v <= sqrt(2(u+v)),
  // which only gives kappa <= 2 sqrt(ln(8r/eta)).
  int points = 0, failAlg2 = 0, fail8r = 0, bChanges = 0;
  long double worst = 0;
  std::string worstAt;
  for (double r : {2.0, 4.0, 6.0, 8.0, 12.0, 16.0, 100.0, 1e4, 1e6, 1e9})
    for (double eps : {0.01, 0.05, 0.3, 0.8, 1.0})
    {
      const long double eta = eps / 96.0L;
      const long double kappa = std::sqrt(2 * std::log(8 / eta)) + std::sqrt(2 * std::log((long double)r));
      const long double need = 2 * std::sqrt(std::log(kAlg2RFactor * r / eta));  // Alg. 2 (8r/eta since 2026-11-14)
      const long double need8 = 2 * std::sqrt(std::log(8 * r / eta));
      ++points;
      if (kappa > need)
      {
        ++failAlg2;
        if (kappa / need > worst)
        {
          worst = kappa / need;
          worstAt = "r=" + std::to_string(r) + " eps=" + std::to_string(eps);
        }
      }
      fail8r += kappa > need8 * (1 + 1e-15L);
      // Would b change with sqrt(ln(8r/eta))?  (n = 2, m = 2, gamma_min on a grid)
      for (double g : {1e-3, 0.0137, 0.5, 2.9})
      {
        const int b4 = int(std::ceil(std::max(0.0L, std::log10(64 * 2 * 2 / eta * need / 2 / g))));
        const int b8 = int(std::ceil(std::max(0.0L, std::log10(64 * 2 * 2 / eta * need8 / 2 / g))));
        bChanges += b4 != b8;
      }
    }
  MESSAGE("kappa > 2 sqrt(ln(F r/eta)) (F = Alg. 2 factor) at " << failAlg2 << "/" << points
          << " (r, eps) points, worst ratio = " << double(worst) << " at "
          << worstAt << "; with ln(8r/eta): " << fail8r << " failures; b would change at "
          << bChanges << "/" << 4 * points << " (r, eps, gamma_min) points");
  CHECK(fail8r == 0);
  CHECK_MESSAGE(failAlg2 == 0,
                "lem:nonconv: (1+kappa/gamma)^n <= (1+eta/(32mn))^n does not follow for "
                "Alg. 2's gamma; needs sqrt(ln(8r/eta))");
}

TEST_CASE("[T3][proof] ttc's b meets lem:nonconv's kappa/(10^b gamma_min) <= eta/(32 m n)")
{
  // What the proof of lem:lattice needs from the b ttc actually uses: every
  // 10^b P_i contains a ball of radius 10^b gamma_min with
  // kappa / (10^b gamma_min) <= eta / (32 m n).  ceil() in Alg. 2 usually
  // covers the <= 2.8% gap above, except when log10(gamma/gamma_min) sits just
  // below an integer.  gamma_min on a fine log grid (1000 points per decade).
  int points = 0, fail = 0;
  std::string first;
  for (double n : {1.0, 2.0, 3.0, 10.0, 40.0})
    for (double m : {1.0, 2.0, 7.0, 100.0})
      for (double rPer : {2.0, 3.0, 8.0})
        for (double eps : {0.05, 0.3, 0.8, 1.0})
          for (int k = 0; k < 3000; ++k)
          {
            const double gHat = std::pow(10.0, -3.0 + k / 1000.0);
            const long double r = std::max(2.0L, (long double)rPer * m);
            const long double eta = eps / 96.0L;
            const long double kappa =
                std::sqrt(2 * std::log(8 / eta)) + std::sqrt(2 * std::log(r));
            const int b = ttc::latticePrecision(std::size_t(n), std::size_t(r), std::size_t(m),
                                                  gHat, eps / 96.0);
            const long double radius = std::pow(10.0L, (long double)b) * gHat;
            ++points;
            if (kappa / radius > eta / (32 * m * n) * (1 + 1e-12L))
            {
              if (fail++ == 0)
                first = "n=" + std::to_string(n) + " m=" + std::to_string(m) + " r="
                        + std::to_string(double(r)) + " eps=" + std::to_string(eps)
                        + " gamma_min=" + std::to_string(gHat) + " b=" + std::to_string(b);
            }
          }
  MESSAGE("kappa/(10^b gamma_min) > eta/(32mn) at " << fail << "/" << points
                                                    << " points; first: " << first);
  CHECK(fail == 0);
}

TEST_CASE("[T3][proof] lem:nonconv conclusion from the proof's chain with Alg. 2's gamma")
{
  // vol Q' <= U vol Q, U = 1 + m((1+k/g)^n - 1); vol Q'' >= L vol Q, L = 1 - m(1-(1-k/g)^n);
  // (1-eta/4)|Q cap Z| ... gives vol Q in [(1-eta/4)/U, 1/((1-eta/4) L)] |Q cap Z|.
  int bad = 0, points = 0;
  for (double n : {1.0, 2.0, 3.0, 10.0, 40.0, 1000.0})
    for (double m : {1.0, 2.0, 10.0, 1000.0, 1e6})
      for (double rPer : {2.0, 4.0, 100.0})
        for (double eps : {0.01, 0.3, 0.8, 1.0})
        {
          const long double r = std::max(2.0L, (long double)rPer * m);
          const long double eta = eps / 96.0L;
          const long double kappa =
              std::sqrt(2 * std::log(8 / eta)) + std::sqrt(2 * std::log(r));
          const long double gamma = alg2Gamma(n, r, m, eps);
          const long double x = kappa / gamma;
          const long double U = 1 + m * std::expm1(n * std::log1p(x));
          const long double L = 1 - m * (-std::expm1(n * std::log1p(-x)));
          const long double lo = (1 - eta / 4) / U, hi = 1 / ((1 - eta / 4) * L);
          ++points;
          const bool ok = x <= 1 && L > 0 && lo >= 1 - eta && hi <= 1 + eta;
          if (!ok)
          {
            ++bad;
            MESSAGE("n=" << n << " m=" << m << " r=" << double(r) << " eps=" << eps
                         << ": vol/|Q cap Z| in [" << double(lo) << ", " << double(hi)
                         << "], eta=" << double(eta));
          }
        }
  MESSAGE(points << " points");
  CHECK(bad == 0);
}

TEST_CASE("[T3][counter] exact lattice counter agrees with closed forms")
{
  for (const UnionCase& U : cases())
  {
    const int maxB = U.n == 2 ? 3 : 2;
    for (int b = 0; b <= maxB; ++b)
    {
      i128 enumerated = 0;
      REQUIRE(enumerate(U, b, 1e9, enumerated));
      if (U.formula) CHECK_MESSAGE(i128str(enumerated) == i128str(U.formula(b)), U.name << " b=" << b);
      if (allBoxes(U)) CHECK_MESSAGE(i128str(enumerated) == i128str(countBoxes(U.polys, U.n, b)), U.name << " b=" << b);
    }
  }
  // 2D simplex: C(t+2, 2).
  UnionCase T{"tri", 2, {simplexM(2, kMilli)}, 0.5L, {}};
  for (int b = 0; b <= 3; ++b)
  {
    i128 c = 0;
    REQUIRE(enumerate(T, b, 1e9, c));
    CHECK(i128str(c) == i128str(binom(i128(pow10ll(b)) + 2, 2)));
  }
  // Closed boundaries: the two triangles of the square share the diagonal.
  i128 c = 0;
  REQUIRE(enumerate(cases()[5], 1, 1e9, c));
  CHECK(i128str(c) == "121");
}

TEST_CASE("[T3][lattice] lem:lattice at ttc's b on small 2D/3D unions; slack of GetPrecision")
{
  const double budget = 1.5e8;  // columns x polytopes per enumeration
  int checkedCases = 0, bad = 0, skipped = 0;
  for (double eps : {0.8, 0.3})
  {
    const double eta = eps / 96.0;
    for (const UnionCase& U : cases())
    {
      const int n = U.n;
      const std::size_t m = U.polys.size();
      std::size_t r = 0;
      double gMin = std::numeric_limits<double>::infinity();
      for (const auto& P : U.polys)
      {
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        toDouble(P, n, A, b);
        const ttc::InnerBall ball = ttc::certifiedInnerBall(A, b);
        REQUIRE(ball.radius > 0);
        gMin = std::min(gMin, ball.radius);
        r += std::size_t(A.rows());
      }
      const int bTtc = ttc::latticePrecision(n, r, m, gMin, eps / 96.0);
      CHECK(bTtc == int(std::ceil(alg2LogRatio(n, r, m, gMin, eps))));

      // Relative error of the lattice count for b = 0..bTtc.
      std::vector<long double> err(bTtc + 1, std::numeric_limits<long double>::quiet_NaN());
      std::string row;
      for (int b = 0; b <= bTtc; ++b)
      {
        i128 cnt = 0;
        if (!exactCount(U, b, budget, cnt)) continue;
        err[b] = (long double)cnt / std::pow(10.0L, (long double)b * n) / U.vol - 1.0L;
        char buf[64];
        std::snprintf(buf, sizeof buf, " b=%d:%+.2Le", b, err[b]);
        row += buf;
      }
      if (std::isnan(double(err[bTtc])))
      {
        ++skipped;
        MESSAGE(U.name << " eps=" << eps << ": b_ttc=" << bTtc << " too large to enumerate;"
                       << row);
        continue;
      }
      ++checkedCases;
      // lem:lattice at ttc's b.
      if (!(std::abs(err[bTtc]) <= eta))
      {
        ++bad;
        MESSAGE(U.name << " eps=" << eps << ": |err| = " << double(std::abs(err[bTtc]))
                       << " > eta = " << eta << " at b=" << bTtc);
      }
      // Smallest b* such that every known b' in [b*, bTtc] is eta-accurate.
      int bStar = bTtc;
      for (int b = bTtc; b >= 0; --b)
      {
        if (std::isnan(double(err[b]))) continue;
        if (std::abs(err[b]) <= eta) bStar = b;
        else break;
      }
      const double T = ttc::lraSampleThreshold(eps, 0.2, m);
      MESSAGE(U.name << " (n=" << n << ", m=" << m << ", r=" << r << ", gamma_min=" << gMin
                     << ") eps=" << eps << " eta=" << eta << ": b_ttc=" << bTtc
                     << " b*=" << bStar << " (slack " << bTtc - bStar << " digits; C x10^"
                     << (bTtc - bStar) * n << ", log10 C = " << std::log10(T) + bTtc * n
                     << " vs " << std::log10(T) + bStar * n << ", GMP digits " << bTtc + 32
                     << " vs " << bStar + 32 << ", first-polytope halvings ~"
                     << int(bTtc * n * std::log2(10.0)) << " vs "
                     << int(bStar * n * std::log2(10.0)) << ");" << row);
    }
  }
  MESSAGE(checkedCases << " (case, eps) pairs checked at b_ttc, " << skipped << " skipped");
  CHECK(bad == 0);
  CHECK(skipped == 0);
}

TEST_CASE("[T3][cost] GMP membership test time vs digits b + 32")
{
  // SampleMultiset<MpFloat>::removeInside over Thresh(0.8, 0.2) 2D points, as
  // phase 2 runs it once per polytope; digits for b = 5..20.
  const double T = ttc::lraSampleThreshold(0.8, 0.2, 2);
  Eigen::MatrixXd A(4, 2);
  A << 1, 0, -1, 0, 0, 1, 0, -1;
  Eigen::VectorXd bb(4);
  bb << 3, -2, 3, -2;  // [2,3]^2: misses every stored point
  std::string row;
  for (int b : {5, 7, 10, 15, 20, 40})
  {
    ttc::MpFloat::default_precision(unsigned(b + 32));
    ttc::SampleMultiset<ttc::MpFloat> store;
    std::mt19937 g(5);
    std::uniform_real_distribution<double> U(0, 1);
    const ttc::MpFloat scale = pow(ttc::MpFloat(10), b);
    for (int i = 0; i < int(T); ++i)
    {
      Eigen::Matrix<ttc::MpFloat, Eigen::Dynamic, 1> v(2);
      v(0) = floor(ttc::MpFloat(U(g)) * scale) / scale;
      v(1) = floor(ttc::MpFloat(U(g)) * scale) / scale;
      store.add(v);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t removed = store.removeInside(A, bb);
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(removed == 0);
    char buf[96];
    std::snprintf(buf, sizeof buf, " b=%d (%d digits): %.1f ms;", b, b + 32, 1e3 * secs);
    row += buf;
  }
  MESSAGE("removeInside over " << int(T) << " points:" << row);
}
