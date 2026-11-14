// T2 of text-for-agents/ttc-union-test-plan.md: error injection into the
// union estimator (Alg. 1 lines 6-22, ttc::runUnionAlgorithm).
//
// Exact components: every test union is a list of 2D parallelotopes
// P = {M u + t : lo <= u <= hi}.  The exact lattice sampler draws u uniform on
// the box inflated by 10^-b ||W_i||_1 (W = M^-1), maps x = M u + t, rounds each
// coordinate stochastically to 10^-b Z and rejects against P.  The rounding
// kernel of a lattice point k has support k + [-s, s]^n, which lies inside the
// inflated body for every k in P, so the accepted law is exactly uniform on
// P cap L (L = 10^-b Z^n, membership by the store's own double test A x - b <= 0).
// |P_i cap L|, |Q cap L| and |D_i cap L| are enumerated exactly.
//
// The bias model.  A sample drawn for P_i survives to the output iff it lies in
//   D_i = P_i \ (P_{i+1} u ... u P_m)          (removeInside, line 8-10),
// and every halving/thinning step is an independent Bernoulli(1/2) thinning
// with p halved (lines 11-19), so |X|/(p C) is a martingale and, exactly,
//   E[est] = sum_i Vhat_i * Q_i(D_i)                                    (*)
// where Vhat_i is the volume handed to the loop and Q_i the sampler's law on
// P_i.  With exact components (Vhat_i = |P_i cap L| 10^-bn, Q_i uniform) this is
// sum_i |D_i cap L| 10^-bn = |Q cap L| 10^-bn.  Consequences checked below:
//  (a) volume errors Vhat_i = V_i (1 + e_i): relative bias = sum e_i |D_i|/|Q|,
//      linear in e, = e for a fixed sign, independent of the overlap;
//  (b) a sampler with TV distance eta_i from uniform: |Q_i(D_i) - U_i(D_i)| <=
//      eta_i, so |bias| <= eta * rho, rho = sum V_i / |Q| (overlap ratio, <= m).
//      The manuscript's pointwise (1 +- eps') law gives |bias| <= eps' instead;
//      TV <= eps'/m keeps the TV bias below eps' for every overlap pattern.
// Every (*) prediction is checked with a z-test on the per-run residual
// est/pred - 1 (|z| <= 4: per-check alpha 6e-5).  Theorem-regime runs (errors
// within the manuscript's budgets) must fail |est/|Q| - 1| > eps/2 (Lemma
// numlattice) with rate <= delta: one-sided binomial test, alpha = 0.01.
//
// (c) ttc's real components (RoundedPolytope oracle + clipping,
// sampleLattice over sampleRestartWalk, b from precisionFromCubes) exactly as
// computeLraVolume wires them (src/volume/lra_volume.cpp), on the same unions.
//
// Slow sweeps:  bin/test_union_robust -ts=slow --no-skip
#include "doctest.h"
#include "common.hpp"
#include "volume/union_estimator.hpp"

#include <chrono>
#include <cstdio>
#include <cstdint>
#include <limits>

namespace
{
using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

constexpr double kZMean = 4.0;   // residual-mean z-tests, alpha ~ 6e-5 each
constexpr double kAlphaBinom = 0.01;

// --------------------------------------------------------------- geometry
struct Para
{
  Mat M, W;     // x = M u + t, u = W (x - t)
  Vec t, lo, hi;
  Mat A;        // H-representation, rows (W_i, -W_i)
  Vec b;
  double volume = 0;  // exact continuous volume
};

Para makePara(const Mat& M, const Vec& t, const Vec& lo, const Vec& hi)
{
  const int n = int(M.rows());
  Para P;
  P.M = M;
  P.W = M.inverse();
  P.t = t;
  P.lo = lo;
  P.hi = hi;
  P.A = Mat(2 * n, n);
  P.b = Vec(2 * n);
  const Vec Wt = P.W * t;
  P.volume = std::abs(M.determinant());
  for (int i = 0; i < n; ++i)
  {
    P.A.row(2 * i) = P.W.row(i);
    P.b(2 * i) = hi(i) + Wt(i);
    P.A.row(2 * i + 1) = -P.W.row(i);
    P.b(2 * i + 1) = -lo(i) - Wt(i);
    P.volume *= hi(i) - lo(i);
  }
  return P;
}

// Sub-parallelotope with u in lo + [a, c] (hi - lo), relative ranges.
Para subPara(const Para& P, const Vec& a, const Vec& c)
{
  const Vec w = P.hi - P.lo;
  return makePara(P.M, P.t, P.lo + a.cwiseProduct(w), P.lo + c.cwiseProduct(w));
}

Mat rot2(double theta)
{
  Mat R(2, 2);
  R << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
  return R;
}
Vec v2(double a, double b)
{
  Vec v(2);
  v << a, b;
  return v;
}
Para box2(double x0, double x1, double y0, double y1)
{
  return makePara(Mat::Identity(2, 2), Vec::Zero(2), v2(x0, y0), v2(x1, y1));
}

// Membership exactly as SampleStore::removeInside evaluates it (tolerance 0).
inline bool inside(const Mat& A, const Vec& b, const Vec& x)
{
  Vec diff = A * x - b;
  return (diff.array() <= 0.0).all();
}
// 2D fast path; identical arithmetic (two products, one sum, no FMA on x86-64).
inline bool inside2(const Mat& A, const Vec& b, double x0, double x1)
{
  for (Eigen::Index r = 0; r < A.rows(); ++r)
    if (A(r, 0) * x0 + A(r, 1) * x1 - b(r) > 0.0) return false;
  return true;
}

struct Union
{
  std::string name;
  std::vector<Para> P;
  int b = 2;  // lattice precision for the exact-component study
};

// The test unions (all 2D).  Axis boxes sit at half-lattice offsets so no
// lattice point lies on a facet; rotated bodies have irrational facets.
Union pair90()
{
  // two boxes overlapping in 90% of their area
  return {"pair90", {box2(0.005, 1.605, 0.005, 1.605), box2(0.165, 1.765, 0.005, 1.605)}, 2};
}
Union chain8()
{
  // 8 boxes shifted by 1/10 of the side: rho = 4.7, D_i = left strip (10%)
  Union U{"chain8", {}, 2};
  for (int i = 0; i < 8; ++i)
    U.P.push_back(box2(0.005 + 0.16 * i, 1.605 + 0.16 * i, 0.005, 1.605));
  return U;
}
Union nestedBigFirst()
{
  return {"nested", {box2(-1.995, 2.005, -1.995, 2.005), box2(-0.995, 1.005, -0.995, 1.005)}, 2};
}
Union slabs3()
{
  // three parallel skinny slabs 40 x 0.4 (aspect 100), shifted 0.1 across,
  // rotated by 30 degrees: rho = 2, D_1, D_2 = 0.1-wide side strips.
  Union U{"slabs3", {}, 2};
  const Mat R = rot2(M_PI / 6);
  const Vec t = v2(0.12345, 0.0567);
  for (int i = 0; i < 3; ++i)
    U.P.push_back(makePara(R, t, v2(0.0, 0.1 * i), v2(40.0, 0.4 + 0.1 * i)));
  return U;
}
Union crossX()
{
  // two skinny slabs 40 x 0.4 crossing at 60 degrees in the middle
  Union U{"crossX", {}, 2};
  const Vec t = v2(0.12345, 0.0567);
  U.P.push_back(makePara(rot2(0.0), t, v2(-20, -0.2), v2(20, 0.2)));
  U.P.push_back(makePara(rot2(M_PI / 3), t, v2(-20, -0.2), v2(20, 0.2)));
  return U;
}

// Real-component-only unions (no lattice enumeration):
Union slabsThin()
{
  // aspect 1000: three slabs 40 x 0.04 shifted 0.004 across (D_1, D_2 are
  // 10% boundary strips), rotated by 30 degrees.
  Union U{"slabsThin", {}, 3};
  const Mat R = rot2(M_PI / 6);
  const Vec t = v2(0.12345, 0.0567);
  for (int i = 0; i < 3; ++i)
    U.P.push_back(makePara(R, t, v2(0.0, 0.004 * i), v2(40.0, 0.04 + 0.004 * i)));
  return U;
}
Union chain5d()
{
  // 5 unit cubes in R^5 shifted 0.1 along x0: D_i = 10% end slices.
  Union U{"chain5d", {}, 0};
  for (int i = 0; i < 5; ++i)
  {
    Vec lo = Vec::Zero(5), hi = Vec::Ones(5);
    lo(0) += 0.1 * i, hi(0) += 0.1 * i;
    U.P.push_back(makePara(Mat::Identity(5, 5), Vec::Zero(5), lo, hi));
  }
  return U;
}
Union skinny5d()
{
  // the same chain, mapped by a random rotation times diag(1,1,1,1,0.01)
  // (aspect 100) -- rotated skinny overlapping bodies in R^5.
  Union U{"skinny5d", {}, 0};
  Vec d = Vec::Ones(5);
  d(4) = 0.01;
  const Mat M = vt::randomRotation(5, 77) * d.asDiagonal();
  const Vec t = Vec::LinSpaced(5, 0.1, 0.5);
  for (int i = 0; i < 5; ++i)
  {
    Vec lo = Vec::Zero(5), hi = Vec::Ones(5);
    lo(0) += 0.1 * i, hi(0) += 0.1 * i;
    U.P.push_back(makePara(M, t, lo, hi));
  }
  return U;
}

// ------------------------------------------------------- exact enumeration
// A probe body S (subset of P_owner) whose |S cap L| and |S cap D_owner cap L|
// are needed for the bias predictions.
struct Probe
{
  const Para* S;
  int owner;
  long long all = 0, inD = 0;
};
struct Counts
{
  std::vector<long long> P, D;  // |P_i cap L|, |D_i cap L|
  long long Q = 0;              // |Q cap L|
  double scaleN = 1;            // 10^{bn}
  double latticeVolume(long long c) const { return double(c) / scaleN; }
  double f(int i) const { return double(D[i]) / double(P[i]); }
};

Counts enumerate(const Union& U, std::vector<Probe>& probes)
{
  const int m = int(U.P.size());
  const double scale = std::pow(10.0, U.b);
  double lo0 = 1e300, lo1 = 1e300, hi0 = -1e300, hi1 = -1e300;
  for (const auto& P : U.P)
    for (int c = 0; c < 4; ++c)
    {
      Vec u(2);
      u << (c & 1 ? P.hi(0) : P.lo(0)), (c & 2 ? P.hi(1) : P.lo(1));
      Vec x = P.M * u + P.t;
      lo0 = std::min(lo0, x(0)), hi0 = std::max(hi0, x(0));
      lo1 = std::min(lo1, x(1)), hi1 = std::max(hi1, x(1));
    }
  Counts C;
  C.P.assign(m, 0);
  C.D.assign(m, 0);
  C.scaleN = scale * scale;
  const long long k0a = (long long)std::floor(lo0 * scale) - 2, k0b = (long long)std::ceil(hi0 * scale) + 2;
  const long long k1a = (long long)std::floor(lo1 * scale) - 2, k1b = (long long)std::ceil(hi1 * scale) + 2;
  std::vector<char> in(m);
  for (long long k0 = k0a; k0 <= k0b; ++k0)
    for (long long k1 = k1a; k1 <= k1b; ++k1)
    {
      const double x0 = double(k0) / scale, x1 = double(k1) / scale;
      int last = -1;
      for (int i = 0; i < m; ++i)
      {
        in[i] = inside2(U.P[i].A, U.P[i].b, x0, x1);
        if (in[i]) ++C.P[i], last = i;
      }
      if (last < 0) continue;
      ++C.Q;
      ++C.D[last];
      for (auto& pr : probes)
        if (in[pr.owner] && inside2(pr.S->A, pr.S->b, x0, x1))
        {
          ++pr.all;
          if (last == pr.owner) ++pr.inD;
        }
    }
  return C;
}

// --------------------------------------------------------- exact samplers
// One exactly uniform point of S cap P cap L (S subset of P), 2D.
Vec drawLattice(const Para& S, const Para& P, int b, std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  const double s = std::pow(10.0, -b), scale = std::pow(10.0, b);
  const double i0 = s * S.W.row(0).lpNorm<1>(), i1 = s * S.W.row(1).lpNorm<1>();
  const double ulo0 = S.lo(0) - i0, uw0 = S.hi(0) - S.lo(0) + 2 * i0;
  const double ulo1 = S.lo(1) - i1, uw1 = S.hi(1) - S.lo(1) + 2 * i1;
  const double m00 = S.M(0, 0), m01 = S.M(0, 1), m10 = S.M(1, 0), m11 = S.M(1, 1);
  for (;;)
  {
    const double u0 = ulo0 + uw0 * unit(rng), u1 = ulo1 + uw1 * unit(rng);
    const double x[2] = {m00 * u0 + m01 * u1 + S.t(0), m10 * u0 + m11 * u1 + S.t(1)};
    double k[2];
    for (int j = 0; j < 2; ++j)
    {
      const double raw = x[j] * scale;
      double base = std::floor(raw);
      if (unit(rng) < raw - base) base += 1;
      k[j] = base;
    }
    const double x0 = k[0] / scale, x1 = k[1] / scale;
    if (inside2(S.A, S.b, x0, x1) && inside2(P.A, P.b, x0, x1)) return v2(x0, x1);
  }
}

// Exact continuous uniform point of P.
Vec drawContinuous(const Para& P, std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  Vec u(P.lo.size());
  for (Eigen::Index i = 0; i < u.size(); ++i) u(i) = P.lo(i) + (P.hi(i) - P.lo(i)) * unit(rng);
  return P.M * u + P.t;
}

// The sampler law for polytope i: a mixture of uniform laws on sub-bodies,
// optionally reweighted pointwise: accept with probability wD on D_i, 1 off it.
struct Law
{
  std::vector<std::pair<double, Para>> comps;  // (weight, S); empty = exact uniform
  double acceptInD = 1.0;
};

bool inD(const Union& U, int i, const Vec& x)
{
  for (std::size_t j = i + 1; j < U.P.size(); ++j)
    if (inside(U.P[j].A, U.P[j].b, x)) return false;
  return true;
}

// --------------------------------------------------------- statistics
double binomUpper(int k, int n, double p)  // P(Bin(n,p) >= k)
{
  if (k <= 0) return 1.0;
  double s = 0;
  for (int j = k; j <= n; ++j)
    s += std::exp(std::lgamma(n + 1.0) - std::lgamma(j + 1.0) - std::lgamma(n - j + 1.0)
                  + j * std::log(p) + (n - j) * std::log1p(-p));
  return std::min(1.0, s);
}

struct Summary
{
  int R = 0;
  double meanRel = 0, sdRel = 0;   // est/exact - 1
  double meanRes = 0, seRes = 0;   // est/pred - 1
  int failEps = 0, failHalf = 0;   // |est/exact - 1| > eps, > eps/2
  double seconds = 0;
  double zRes() const { return seRes > 0 ? meanRes / seRes : 0.0; }
};

Summary summarize(const std::vector<double>& rel, const std::vector<double>& res, double eps)
{
  Summary S;
  S.R = int(rel.size());
  for (double r : rel) S.meanRel += r;
  S.meanRel /= S.R;
  for (double r : rel) S.sdRel += (r - S.meanRel) * (r - S.meanRel);
  S.sdRel = std::sqrt(S.sdRel / std::max(1, S.R - 1));
  double sdRes = 0;
  for (double r : res) S.meanRes += r;
  S.meanRes /= S.R;
  for (double r : res) sdRes += (r - S.meanRes) * (r - S.meanRes);
  S.seRes = std::sqrt(sdRes / std::max(1, S.R - 1)) / std::sqrt(double(S.R));
  for (double r : rel)
  {
    if (!(std::abs(r) <= eps)) ++S.failEps;
    if (!(std::abs(r) <= eps / 2)) ++S.failHalf;
  }
  return S;
}

// ----------------------------------------------- the exact-component driver
struct Study
{
  Union U;
  Counts C;
  std::vector<Probe> probes;  // filled per law
  explicit Study(Union u) : U(std::move(u))
  {
    std::vector<Probe> none;
    C = enumerate(U, none);
  }
  double exact() const { return C.latticeVolume(C.Q); }
  double V(int i) const { return C.latticeVolume(C.P[i]); }
  double rho() const
  {
    double s = 0;
    for (std::size_t i = 0; i < U.P.size(); ++i) s += V(int(i));
    return s / exact();
  }
};

// Q_i(D_i) of a law from enumerated probe counts.
struct LawSet
{
  std::vector<Law> laws;       // per polytope
  std::vector<double> qD;      // Q_i(D_i), exact
};

// Uniform law per polytope.
LawSet uniformLaws(const Study& St)
{
  LawSet L;
  for (std::size_t i = 0; i < St.U.P.size(); ++i)
  {
    L.laws.push_back(Law{});
    L.qD.push_back(St.C.f(int(i)));
  }
  return L;
}

// Mixture (1 - w) U + w Uniform(S_i), S_i = sub-box [a, c] of P_i (relative
// u-coordinates), with w set so that TV(law, U) = tv exactly; tv = -1 means
// w = 1 (the whole law is uniform on S_i, a "short walk").
LawSet mixtureLaws(const Study& St, const Vec& a, const Vec& c, double tv,
                   double* tvOut = nullptr)
{
  const int m = int(St.U.P.size());
  std::vector<Para> subs;
  for (int i = 0; i < m; ++i) subs.push_back(subPara(St.U.P[i], a, c));
  std::vector<Probe> probes;
  for (int i = 0; i < m; ++i) probes.push_back({&subs[i], i});
  enumerate(St.U, probes);
  LawSet L;
  double tvMax = 0;
  for (int i = 0; i < m; ++i)
  {
    const double fracS = double(probes[i].all) / double(St.C.P[i]);
    const double tvS = 1.0 - fracS;  // TV(Uniform(S), Uniform(P))
    const double w = tv < 0 ? 1.0 : std::min(1.0, tv / tvS);
    tvMax = std::max(tvMax, w * tvS);
    const double g = double(probes[i].inD) / double(probes[i].all);
    Law law;
    law.comps.push_back({1.0 - w, St.U.P[i]});
    law.comps.push_back({w, subs[i]});
    L.laws.push_back(law);
    L.qD.push_back((1.0 - w) * St.C.f(i) + w * g);
  }
  if (tvOut) *tvOut = tvMax;
  return L;
}

// The manuscript's pointwise law at the edge of eq. (almostuniform):
// Pr(x) = (1 - eps' a)/|P| on D_i, (1 + eps' c)/|P| off D_i, with f a = (1-f) c,
// a, c <= 1.  It under-weights every D_i as much as the lemma allows.
LawSet pointwiseLaws(const Study& St, double epsP)
{
  LawSet L;
  for (std::size_t i = 0; i < St.U.P.size(); ++i)
  {
    const double f = St.C.f(int(i));
    double a = 0, c = 0;
    if (f > 0 && f < 1)
    {
      a = f <= 0.5 ? 1.0 : (1.0 - f) / f;
      c = f * a / (1.0 - f);
    }
    Law law;
    law.acceptInD = (1.0 - epsP * a) / (1.0 + epsP * c);
    L.laws.push_back(law);
    L.qD.push_back(f * (1.0 - epsP * a));
  }
  return L;
}

// volume factors (1 + e_i) per run
using VolFn = std::function<std::vector<double>(unsigned seed)>;
VolFn exactVolumes(int m)
{
  return [m](unsigned) { return std::vector<double>(m, 1.0); };
}
VolFn fixedError(int m, double e)
{
  return [m, e](unsigned) { return std::vector<double>(m, 1.0 + e); };
}
VolFn randomSignError(int m, double e)
{
  return [m, e](unsigned seed) {
    std::mt19937 g(seed * 7919u + 17u);
    std::vector<double> f(m);
    for (auto& x : f) x = (g() & 1u) ? 1.0 + e : 1.0 - e;
    return f;
  };
}

Summary runStudy(const Study& St, const LawSet& L, const VolFn& volFn, double eps, double delta,
                 int R, unsigned seedBase)
{
  const int m = int(St.U.P.size());
  const double thresh = ttc::lraSampleThreshold(eps, delta, m);
  const double log2Scale = std::log2(thresh) + double(St.U.b) * 2.0 * std::log2(10.0);
  std::vector<double> rel, res;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < R; ++r)
  {
    const unsigned seed = seedBase + 1000003u * unsigned(r);
    const auto fac = volFn(seed);
    std::vector<ttc::PolytopeVolume> kept;
    double pred = 0;
    for (int i = 0; i < m; ++i)
    {
      kept.push_back({St.U.P[i].A, St.U.P[i].b, St.V(i) * fac[i]});
      pred += kept.back().volume * L.qD[i];
    }
    std::mt19937_64 rng(seed ^ 0x5bd1e995ULL);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto gen = [&](const Mat& A, const Vec&, long long n, double& t) {
      t = 0;
      int i = 0;
      while (&kept[i].A != &A) ++i;
      const Law& law = L.laws[i];
      std::vector<Vec> out;
      out.reserve(n);
      while ((long long)out.size() < n)
      {
        const Para* S = &St.U.P[i];
        if (!law.comps.empty())
        {
          double w = unit(rng);
          for (const auto& c : law.comps)
          {
            S = &c.second;
            if ((w -= c.first) < 0) break;
          }
        }
        Vec x = drawLattice(*S, St.U.P[i], St.U.b, rng);
        if (law.acceptInD < 1.0 && inD(St.U, i, x) && unit(rng) >= law.acceptInD) continue;
        out.push_back(std::move(x));
      }
      return out;
    };
    ttc::VolumeComputationResult result;
    ttc::runUnionAlgorithm<double>(kept, thresh, seed, log2Scale, gen, result, nullptr);
    rel.push_back(result.volumeEstimate / St.exact() - 1.0);
    res.push_back(result.volumeEstimate / pred - 1.0);
  }
  Summary S = summarize(rel, res, eps);
  S.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return S;
}

void printRow(const char* what, const Study& St, double eps, double inj,
              double predBias, const Summary& S)
{
  std::printf("  %-8s %-26s eps=%.1f inj=%+.4f predBias=%+.4f meanRel=%+.4f sd=%.4f "
              "res z=%+5.2f fail(eps/2)=%3d/%d fail(eps)=%3d/%d (%.1fs)\n",
              St.U.name.c_str(), what, eps, inj, predBias, S.meanRel, S.sdRel, S.zRes(),
              S.failHalf, S.R, S.failEps, S.R, S.seconds);
  std::fflush(stdout);
}

double predictedBias(const Study& St, const LawSet& L, const std::vector<double>& fac)
{
  double s = 0;
  for (std::size_t i = 0; i < St.U.P.size(); ++i) s += St.V(int(i)) * fac[i] * L.qD[i];
  return s / St.exact() - 1.0;
}

// The theorem-regime check: failure rate against eps/2 (Lemma numlattice) is
// <= delta (one-sided binomial, alpha 0.01).
void checkTheorem(const Summary& S, double delta)
{
  const double pv = binomUpper(S.failHalf, S.R, delta);
  CHECK_MESSAGE(pv >= kAlphaBinom, "failures(eps/2) " << S.failHalf << "/" << S.R
                                                       << " > delta=" << delta
                                                       << " (binomial p=" << pv << ")");
}
void checkPrediction(const Summary& S)
{
  CHECK_MESSAGE(std::abs(S.zRes()) <= kZMean,
                "mean residual vs E[est]=sum Vhat_i Q_i(D_i): " << S.meanRes << " (z=" << S.zRes()
                                                                << ")");
}

}  // namespace

// ===========================================================================
// (0) Baseline: exact components meet Lemma numlattice; E[est] = |Q cap L|.
// ===========================================================================
TEST_CASE("robust [baseline]: exact components on every test union")
{
  std::printf("== baseline (exact volumes, exact uniform lattice sampler)\n");
  for (auto U : {pair90(), chain8(), nestedBigFirst(), slabs3(), crossX()})
  {
    Study St(U);
    std::printf("  %s: m=%zu |Q cap L|=%lld rho=%.3f f_i=", St.U.name.c_str(), St.U.P.size(),
                St.C.Q, St.rho());
    for (std::size_t i = 0; i < St.U.P.size(); ++i) std::printf("%.3f ", St.C.f(int(i)));
    std::printf("\n");
    const auto L = uniformLaws(St);
    for (auto ed : {std::pair<double, double>{0.8, 0.2}, {0.3, 0.1}})
    {
      const int R = ed.first > 0.5 ? 60 : 12;
      Summary S = runStudy(St, L, exactVolumes(int(St.U.P.size())), ed.first, ed.second, R, 11);
      printRow("exact", St, ed.first, 0, 0, S);
      CAPTURE(St.U.name);
      checkPrediction(S);
      checkTheorem(S, ed.second);
    }
  }
}

// ===========================================================================
// (a) Volume errors: bias = sum e_i |D_i|/|Q|, linear; = e for a fixed sign.
// ===========================================================================
namespace
{
void volumeSweep(const Union& U, double eps, double delta, int R, bool full)
{
  Study St(U);
  const int m = int(St.U.P.size());
  const double epsP = eps / 12.0;
  const auto L = uniformLaws(St);
  std::vector<double> mults = full ? std::vector<double>{0.5, 1, 2, 5} : std::vector<double>{1, 5};
  for (double k : mults)
  {
    const double e = k * epsP;
    for (int sign : {+1, -1, 0})
    {
      VolFn fn = sign == 0 ? randomSignError(m, e) : fixedError(m, sign * e);
      Summary S = runStudy(St, L, fn, eps, delta, R, 101 + unsigned(k * 10) + 7u * unsigned(sign + 1));
      // average predicted bias over the runs' volume factors
      double pb = 0;
      for (int r = 0; r < R; ++r) pb += predictedBias(St, L, fn(101 + unsigned(k * 10) + 7u * unsigned(sign + 1) + 1000003u * unsigned(r)));
      pb /= R;
      char what[64];
      std::snprintf(what, sizeof what, "vol x(1%s%.1feps')", sign > 0 ? "+" : sign < 0 ? "-" : "+-", k);
      printRow(what, St, eps, sign == 0 ? e : sign * e, pb, S);
      CAPTURE(St.U.name);
      CAPTURE(e);
      CAPTURE(sign);
      checkPrediction(S);
      // In the C4 budget (|e| <= eps') the theorem must hold.
      if (k <= 1.0) checkTheorem(S, delta);
      // Fixed sign: the bias is exactly e, whatever the overlap.
      if (sign != 0) CHECK(std::abs(pb - sign * e) < 1e-9);
    }
  }
}
}  // namespace

TEST_CASE("robust [volume]: volume errors give bias = sum e_i |D_i|/|Q| (eps=0.8)")
{
  std::printf("== (a) volume-error injection\n");
  volumeSweep(chain8(), 0.8, 0.2, 60, false);
  volumeSweep(slabs3(), 0.8, 0.2, 60, false);
}

TEST_CASE("robust [volume]: full sweep e in {eps'/2, eps', 2eps', 5eps'}" * doctest::test_suite("slow") * doctest::skip())
{
  std::printf("== (a) volume-error injection, full sweep\n");
  for (auto U : {pair90(), chain8(), nestedBigFirst(), slabs3()})
  {
    volumeSweep(U, 0.8, 0.2, 200, true);
    volumeSweep(U, 0.3, 0.1, 100, true);
  }
}

// ===========================================================================
// (b) Sampler non-uniformity.
// ===========================================================================
namespace
{
struct LawSpec
{
  const char* name;
  Vec a, c;
};
std::vector<LawSpec> badLaws()
{
  return {{"corner0.1", v2(0, 0), v2(0.1, 0.1)},
          {"corner0.3", v2(0, 0), v2(0.3, 0.3)},
          {"centre0.3", v2(0.35, 0.35), v2(0.65, 0.65)},
          {"lowside0.1", v2(0, 0), v2(1, 0.1)}};
}

// Sensitivity kappa = d bias / d TV for a law, and the worst over the specs.
void tvSweep(const Union& U, double eps, double delta, int R, bool full)
{
  Study St(U);
  const int m = int(St.U.P.size());
  const double epsP = eps / 12.0;
  const double rho = St.rho();
  // Pointwise (1 +- eps') law at the edge of eq. (almostuniform): lemma regime.
  {
    const auto L = pointwiseLaws(St, epsP);
    const double pb = predictedBias(St, L, std::vector<double>(m, 1.0));
    Summary S = runStudy(St, L, exactVolumes(m), eps, delta, R, 501);
    printRow("pointwise(1-eps')onD", St, eps, epsP, pb, S);
    CAPTURE(St.U.name);
    checkPrediction(S);
    checkTheorem(S, delta);
    CHECK(pb >= -epsP - 1e-12);  // the lemma's bound
  }
  // Everything at the manuscript's budget at once: Vhat = V(1-eps') and the
  // pointwise law: bias ~ -2 eps' (lemma regime, must hold).
  {
    const auto L = pointwiseLaws(St, epsP);
    std::vector<double> fac(m, 1.0 - epsP);
    const double pb = predictedBias(St, L, fac);
    Summary S = runStudy(St, L, fixedError(m, -epsP), eps, delta, R, 601);
    printRow("vol(1-eps')+pointwise", St, eps, -epsP, pb, S);
    checkPrediction(S);
    checkTheorem(S, delta);
  }
  for (const auto& spec : badLaws())
  {
    std::vector<double> tvs = full ? std::vector<double>{epsP / m, epsP / rho, epsP, 2 * epsP, 5 * epsP}
                                   : std::vector<double>{epsP / m, epsP};
    for (double tv : tvs)
    {
      double tvReal = 0;
      const auto L = mixtureLaws(St, spec.a, spec.c, tv, &tvReal);
      const double pb = predictedBias(St, L, std::vector<double>(m, 1.0));
      Summary S = runStudy(St, L, exactVolumes(m), eps, delta, R, 701 + unsigned(tv * 1e5));
      char what[64];
      std::snprintf(what, sizeof what, "%s TV=%.4f", spec.name, tvReal);
      printRow(what, St, eps, tvReal, pb, S);
      CAPTURE(St.U.name);
      CAPTURE(spec.name);
      CAPTURE(tv);
      checkPrediction(S);
      // TV bound: |bias| <= TV * rho
      CHECK(std::abs(pb) <= tvReal * rho + 1e-12);
      // The plan's sufficient condition TV <= eps'/m: theorem must hold.
      if (tv <= epsP / m + 1e-15) checkTheorem(S, delta);
    }
  }
}
}  // namespace

TEST_CASE("robust [sampler]: TV mixtures, pointwise law; bias = sum V_i (Q_i-U_i)(D_i)/|Q|")
{
  std::printf("== (b) sampler-law injection\n");
  tvSweep(chain8(), 0.8, 0.2, 40, false);
  tvSweep(slabs3(), 0.8, 0.2, 40, false);
}

TEST_CASE("robust [sampler]: full TV sweep" * doctest::test_suite("slow") * doctest::skip())
{
  std::printf("== (b) sampler-law injection, full sweep\n");
  for (auto U : {pair90(), chain8(), nestedBigFirst(), slabs3()})
  {
    tvSweep(U, 0.8, 0.2, 200, true);
    tvSweep(U, 0.3, 0.1, 60, true);
  }
}

// Short-walk laws: every sample from the central sub-box of relative size r
// (mimics the old restart-at-centre bias).  Prediction must hold; the table
// gives the bias and failure rate of each r.
TEST_CASE("robust [sampler]: short-walk (central sub-box) laws on large overlaps")
{
  std::printf("== (b) short-walk laws\n");
  for (auto U : {chain8(), nestedBigFirst(), slabs3(), crossX()})
  {
    Study St(U);
    const int m = int(St.U.P.size());
    for (double r : {0.95, 0.8, 0.5})
    {
      const Vec a = v2(0.5 - r / 2, 0.5 - r / 2), c = v2(0.5 + r / 2, 0.5 + r / 2);
      double tvReal = 0;
      const auto L = mixtureLaws(St, a, c, -1.0, &tvReal);
      const double pb = predictedBias(St, L, std::vector<double>(m, 1.0));
      Summary S = runStudy(St, L, exactVolumes(m), 0.8, 0.2, 40, 801 + unsigned(r * 100));
      char what[64];
      std::snprintf(what, sizeof what, "short r=%.2f TV=%.3f", r, tvReal);
      printRow(what, St, 0.8, tvReal, pb, S);
      CAPTURE(St.U.name);
      CAPTURE(r);
      checkPrediction(S);
    }
  }
}

// ===========================================================================
// (c) ttc's real components on the same unions.
// ===========================================================================
namespace
{
// ttc's streamSeed (lra_volume.cpp computeLraVolume).
struct StreamSeed
{
  unsigned base, index = 0;
  unsigned operator()()
  {
    std::uint64_t z = base + 0x9E3779B97F4A7C15ULL * (++index);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return static_cast<unsigned>(z ^ (z >> 31));
  }
};

struct RealSetup
{
  std::vector<ttc::PolytopeVolume> kept;  // clipped oracle volumes, as ttc
  std::vector<double> volErr;             // oracle / exact - 1 (before clipping)
  int precision = 0;
  double exact = 0;                       // continuous vol(Q)
  std::vector<double> uD;                 // U_i(D_i), continuous, MC with 4e6 pts
  std::vector<double> uDse;
};

double continuousUnionVolume(const Union& U, std::vector<double>& uD, std::vector<double>& uDse)
{
  // vol(Q) = sum_i vol(P_i) U_i(D_i); U_i(D_i) by exact uniform sampling.
  std::mt19937_64 rng(4242);
  const int N = 2000000 / int(U.P.size());
  double vol = 0;
  uD.clear();
  uDse.clear();
  for (std::size_t i = 0; i < U.P.size(); ++i)
  {
    long long hit = 0;
    for (int k = 0; k < N; ++k)
      if (inD(U, int(i), drawContinuous(U.P[i], rng))) ++hit;
    const double f = double(hit) / N;
    uD.push_back(f);
    uDse.push_back(std::sqrt(std::max(f * (1 - f), 1.0 / N) / N));
    vol += U.P[i].volume * f;
  }
  return vol;
}

RealSetup realSetup(const Union& U, double eps)
{
  RealSetup R;
  const int n = int(U.P[0].A.cols());
  const double epsP = eps / 12.0;
  std::size_t facets = 0;
  double minRadius = std::numeric_limits<double>::infinity();
  for (const auto& P : U.P)
  {
    const auto ball = ttc::certifiedInnerBall(P.A, P.b);
    const double outer = ttc::outerRadius(P.A, P.b, ball.center);
    minRadius = std::min(minRadius, ball.radius);
    facets += P.A.rows();
    const ttc::RoundedPolytope body(P.A, P.b, ball, true);
    double vol = body.volume(ttc::kDefaultVolumeOracle, epsP,
                             ttc::defaultOracleWalkLength(ttc::kDefaultVolumeOracle, n));
    R.volErr.push_back(vol / P.volume - 1.0);
    const double lower = (1.0 - epsP) * std::pow(2.0 * ball.radius / std::sqrt(double(n)), double(n));
    const double upper = (1.0 + epsP) * std::pow(2.0 * outer, double(n));
    vol = std::clamp(vol, lower, upper);
    R.kept.push_back({P.A, P.b, vol});
  }
  R.precision = ttc::precisionFromCubes(n, facets, U.P.size(), minRadius, eps);
  R.exact = continuousUnionVolume(U, R.uD, R.uDse);
  return R;
}

// ttc's GenerateSamples exactly as computeLraVolume (--nogmp) calls it.
auto ttcGenerator(int precision, double eps, double delta, std::size_t m, StreamSeed& ss, int n,
                  bool round = true)
{
  const double epsP = eps / 12.0;
  const double oracleDelta = delta / (4.0 * std::max<double>(1.0, double(m)));
  const unsigned wl = ttc::defaultSampleWalkLength(ttc::SamplerWalk::Billiard, n);
  return [=, &ss](const Mat& A, const Vec& b, long long n, double& time) {
    auto inner = [&](const Mat& Ai, const Vec& bi, long long k, double& t) {
      double ballSeconds = 0;
      t = 0;
      return ttc::sampleRestartWalk(Ai, bi, k, ttc::SamplerWalk::Billiard, wl, ss(), ballSeconds,
                                    round);
    };
    return ttc::sampleLattice<double>(A, b, n, precision, epsP, oracleDelta, ss(), inner, time);
  };
}
}  // namespace

TEST_CASE("robust [ttc]: real oracle volumes are within eps' (C4 budget)")
{
  std::printf("== (c) ttc volume oracle (RoundedPolytope, fixed seed 3) on the test unions\n");
  for (double eps : {0.8, 0.3})
    for (auto U : {pair90(), chain8(), nestedBigFirst(), slabs3(), crossX(), slabsThin(), chain5d(), skinny5d()})
    {
      const double epsP = eps / 12.0;
      double worst = 0;
      int over = 0;
      const int n = int(U.P[0].A.cols());
      for (const auto& P : U.P)
      {
        const auto ball = ttc::certifiedInnerBall(P.A, P.b);
        const ttc::RoundedPolytope body(P.A, P.b, ball, true);
        const double e = body.volume(ttc::kDefaultVolumeOracle, epsP,
                                     ttc::defaultOracleWalkLength(ttc::kDefaultVolumeOracle, n))
                             / P.volume - 1.0;
        worst = std::max(worst, std::abs(e));
        if (std::abs(e) > epsP) ++over;
      }
      std::printf("  %-9s eps=%.1f eps'=%.4f worst |e_i|=%.2e over eps': %d/%zu\n",
                  U.name.c_str(), eps, epsP, worst, over, U.P.size());
      CAPTURE(U.name);
      CAPTURE(eps);
      // Seeded spread: P(|e| > eps') <= delta' = delta/(4m) per call.
      const double delta = eps > 0.5 ? 0.2 : 0.1;
      const double deltaP = delta / (4.0 * U.P.size());
      const int K = n > 2 ? 8 : 16;
      int bad = 0;
      double sum = 0, sum2 = 0;
      for (int s = 0; s < K; ++s)
      {
        const auto& P = U.P.back();
        const auto ball = ttc::certifiedInnerBall(P.A, P.b);
        const ttc::RoundedPolytope body(P.A, P.b, ball, true);
        const double e = body.volumeWithSeed(ttc::kDefaultVolumeOracle, epsP,
                                             ttc::defaultOracleWalkLength(ttc::kDefaultVolumeOracle, n),
                                             1000u + 7919u * s)
                             / P.volume - 1.0;
        sum += e, sum2 += e * e;
        if (std::abs(e) > epsP) ++bad;
      }
      const double mu = sum / K, sd = std::sqrt(std::max(0.0, sum2 / K - mu * mu));
      std::printf("           seeded (last P, K=%d): mean e=%+.4f sd=%.4f |e|>eps': %d/%d (delta'=%.4f)\n",
                  K, mu, sd, bad, K, deltaP);
      CHECK_MESSAGE(binomUpper(bad, K, deltaP) >= kAlphaBinom,
                    "oracle misses eps' in " << bad << "/" << K << " seeded runs, delta'=" << deltaP);
    }
}

// Q_i(D_i) of ttc's sampler vs U_i(D_i): the quantity (*) says matters.
// eq. (almostuniform) implies |Q_i(D)/U_i(D) - 1| <= eps' for every region D.
namespace
{
void samplerMass(const std::vector<Union>& unions, bool round)
{
  const double eps = 0.8, delta = 0.2, epsP = eps / 12.0;
  for (const auto& U : unions)
  {
    RealSetup Rs = realSetup(U, eps);
    const int m = int(U.P.size());
    StreamSeed ss{12345u};
    auto gen = ttcGenerator(Rs.precision, eps, delta, U.P.size(), ss, int(U.P[0].A.cols()), round);
    double unionBias = 0, unionVar = 0;
    for (int i = 0; i < m; ++i)
    {
      if (Rs.uD[i] == 0.0) continue;  // D_i empty: the sampler is irrelevant
      const int N = 20000;
      double t = 0;
      auto pts = gen(U.P[i].A, U.P[i].b, N, t);
      long long hit = 0;
      for (const auto& x : pts)
        if (inD(U, i, x)) ++hit;
      const double q = double(hit) / N;
      const double se = std::sqrt(q * (1 - q) / N + Rs.uDse[i] * Rs.uDse[i]);
      const double ratio = q / Rs.uD[i] - 1.0;
      std::printf("  %-8s P%d b=%d U(D)=%.4f Q(D)=%.4f ratio-1=%+.4f (z vs eps' band %+.2f)\n",
                  U.name.c_str(), i + 1, Rs.precision, Rs.uD[i], q, ratio,
                  (std::abs(q - Rs.uD[i]) - epsP * Rs.uD[i]) / se);
      unionBias += U.P[i].volume * (q - Rs.uD[i]);
      unionVar += U.P[i].volume * U.P[i].volume * se * se;
      CAPTURE(U.name);
      CAPTURE(i);
      CHECK_MESSAGE((std::abs(q - Rs.uD[i]) - epsP * Rs.uD[i]) / se <= kZMean,
                    "sampler mass on D_i off by " << ratio << " (> eps'=" << epsP << ")");
    }
    std::printf("  %-8s sampler-induced union bias %+.4f +- %.4f (eps'=%.4f)\n", U.name.c_str(),
                unionBias / Rs.exact, std::sqrt(unionVar) / Rs.exact, epsP);
  }
}
}  // namespace

TEST_CASE("robust [ttc]: sampleLattice o sampleRestartWalk mass on D_i within eps'")
{
  std::printf("== (c) ttc sampler mass on D_i (N=20000 per polytope)\n");
  samplerMass({pair90(), chain8(), nestedBigFirst(), slabs3(), crossX(), slabsThin(), chain5d(),
               skinny5d()},
              true);
}

// ttc --no-round (and every ttc before 08026a2): the restart walk runs on the
// unrounded body.  Same check; expected to FAIL on the skinny unions.
TEST_CASE("robust [ttc]: --no-round sampler mass on D_i within eps'" * doctest::test_suite("slow")
          * doctest::skip())
{
  std::printf("== (c) ttc --no-round sampler mass on D_i (N=20000 per polytope)\n");
  samplerMass({chain8(), slabs3(), crossX(), slabsThin(), skinny5d()}, false);
}

namespace
{
void realUnion(const Union& U, double eps, double delta, int R)
{
  RealSetup Rs = realSetup(U, eps);
  const int m = int(U.P.size());
  const double thresh = ttc::lraSampleThreshold(eps, delta, m);
  const double log2Scale = std::log2(thresh) + double(Rs.precision) * double(U.P[0].A.cols()) * std::log2(10.0);
  // prediction (*) with the real volumes and the exact uniform sampler
  double pred = 0;
  for (int i = 0; i < m; ++i) pred += Rs.kept[i].volume * Rs.uD[i];
  std::vector<double> rel, res;
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < R; ++r)
  {
    const unsigned seed = 1u + 7919u * unsigned(r);
    StreamSeed ss{seed};
    auto gen = ttcGenerator(Rs.precision, eps, delta, U.P.size(), ss, int(U.P[0].A.cols()));
    ttc::VolumeComputationResult result;
    ttc::runUnionAlgorithm<double>(Rs.kept, thresh, seed, log2Scale, gen, result, nullptr);
    rel.push_back(result.volumeEstimate / Rs.exact - 1.0);
    res.push_back(result.volumeEstimate / pred - 1.0);
  }
  Summary S = summarize(rel, res, eps);
  S.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  double worstE = 0;
  for (double e : Rs.volErr) worstE = std::max(worstE, std::abs(e));
  std::printf("  %-8s eps=%.1f b=%d max|e_vol|=%.4f volBias=%+.4f meanRel=%+.4f sd=%.4f "
              "samplerBias(res)=%+.4f z=%+5.2f fail(eps/2)=%d/%d fail(eps)=%d/%d (%.1fs)\n",
              U.name.c_str(), eps, Rs.precision, worstE, pred / Rs.exact - 1.0, S.meanRel, S.sdRel,
              S.meanRes, S.zRes(), S.failHalf, S.R, S.failEps, S.R, S.seconds);
  std::fflush(stdout);
  CAPTURE(U.name);
  CAPTURE(eps);
  // Theorem: P(|est/vol(Q) - 1| > eps) <= delta.
  CHECK_MESSAGE(binomUpper(S.failEps, S.R, delta) >= kAlphaBinom,
                "failures(eps) " << S.failEps << "/" << S.R << " > delta");
  // Sampler-induced bias (residual vs the exact-sampler prediction) within the
  // eps' budget of eq. (almostuniform).
  CHECK_MESSAGE(std::abs(S.meanRes) - eps / 12.0 <= kZMean * S.seRes,
                "sampler-induced union bias " << S.meanRes << " exceeds eps'");
}
}  // namespace

TEST_CASE("robust [ttc]: full union with ttc's real components (eps=0.8, delta=0.2)")
{
  std::printf("== (c) full union, ttc components\n");
  for (auto U : {pair90(), chain8(), nestedBigFirst(), slabs3(), crossX(), slabsThin()})
    realUnion(U, 0.8, 0.2, 16);
}

TEST_CASE("robust [ttc]: full union with ttc's real components, more seeds" * doctest::test_suite("slow") * doctest::skip())
{
  std::printf("== (c) full union, ttc components (slow)\n");
  for (auto U : {pair90(), chain8(), nestedBigFirst(), slabs3(), crossX(), slabsThin()})
    realUnion(U, 0.3, 0.1, 20);
  for (auto U : {pair90(), chain8(), nestedBigFirst(), slabs3(), crossX(), slabsThin(), chain5d(), skinny5d()})
    realUnion(U, 0.8, 0.2, 100);
  for (auto U : {chain5d(), skinny5d()})
    realUnion(U, 0.3, 0.1, 20);
}
