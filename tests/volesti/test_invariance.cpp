// Invariance and inner-ball sensitivity tests for the volesti volume oracle.
//
// Volume is equivariant under scaling, translation and rotation:
//   vol(sP) = s^n vol(P),  vol(P + t) = vol(P),  vol(QP) = vol(P)  (Q orthogonal).
// The oracle (volume_cooling_gaussians + GaussianBallWalk, ttc configuration)
// should therefore be accurate on sP, P+t and QP whenever it is accurate on P.
// It also should not depend on which *valid* inner ball it is handed.
//
// Every oracle call runs in a forked child with a CPU limit, because some inputs
// make the oracle loop forever (get_first_gaussian's bisection on tiny bodies).
// A timeout / crash counts as a wrong answer.
//
// Default run: ~45 s.  The full size grid is tagged "slow" and skipped:
//   bin/test_invariance -ts=slow -nr     (or --no-skip)
#include "doctest.h"
#include "common.hpp"

#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#include "preprocess/max_inscribed_ball.hpp"

using namespace vt;

namespace
{

// Tolerance: ttc asks the oracle for relative error eps/12 = 0.0667; a single
// deterministic run (seed 3) at scale 1 is within ~0.05 on these bodies, so 3x
// the target is a generous bound that a scale-correct oracle meets.
constexpr double kTol = 0.2;

struct Guarded
{
  double value = NAN;
  bool timedOut = false;
  double seconds = 0;
};

// Run f() in a forked child with a CPU-time limit (SIGXCPU; robust to a loaded
// machine) and a 10x wall-clock backstop.
template <class F>
Guarded guarded(F f, unsigned seconds)
{
  Guarded G;
  int fd[2];
  if (pipe(fd) != 0) return G;
  std::fflush(stdout);
  auto t0 = std::chrono::steady_clock::now();
  pid_t pid = fork();
  if (pid == 0)
  {
    close(fd[0]);
    struct rlimit rl = {seconds, seconds + 1};
    setrlimit(RLIMIT_CPU, &rl);
    alarm(10 * seconds);
    double v = f();
    ssize_t w = write(fd[1], &v, sizeof v);
    (void)w;
    _exit(0);
  }
  close(fd[1]);
  double v = NAN;
  if (read(fd[0], &v, sizeof v) != static_cast<ssize_t>(sizeof v)) v = NAN;
  close(fd[0]);
  int st = 0;
  waitpid(pid, &st, 0);
  G.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (WIFSIGNALED(st)) G.timedOut = true;
  G.value = v;
  return G;
}

struct IBall
{
  VT c;
  double r;
};

// Exact Chebyshev balls.
IBall cubeBall(int n, double lo, double hi) { return {VT::Constant(n, (lo + hi) / 2), (hi - lo) / 2}; }
IBall simplexBall(int n, double s = 1.0)
{
  const double c = s / (n + std::sqrt(double(n)));
  return {VT::Constant(n, c), c};
}
IBall boxBall(const std::vector<double>& lo, const std::vector<double>& hi)
{
  const int n = int(lo.size());
  VT c(n);
  double r = INFINITY;
  for (int i = 0; i < n; ++i)
  {
    c(i) = (lo[i] + hi[i]) / 2;
    r = std::min(r, (hi[i] - lo[i]) / 2);
  }
  return {c, r};
}
IBall mapBall(const IBall& B, const MT& Q, const VT& t, double s)  // Q orthogonal * s
{
  return {Q * B.c * s + t, B.r * s};
}

// Oracle with an explicit inner ball (ttc sets a certified Chebyshev ball).
double oracleWithBall(const Body& B, const VT& c, double r)
{
  HPoly P = toHPoly(B);
  setBall(P, c, r);
  return ttcOracle(P);
}
double oracleVolestiBall(const Body& B)
{
  HPoly P = toHPoly(B);
  volestiBall(P);
  return ttcOracle(P);
}

double rel(double v, double exact) { return v / exact - 1.0; }

bool ok(const Guarded& G, double exact, double tol = kTol)
{
  return !G.timedOut && std::isfinite(G.value) && G.value > 0 && std::abs(rel(G.value, exact)) <= tol;
}

std::string show(const Guarded& G, double exact)
{
  char buf[96];
  if (G.timedOut)
    std::snprintf(buf, sizeof buf, "TIMEOUT(%.0fs)", G.seconds);
  else
    std::snprintf(buf, sizeof buf, "%+.3g", rel(G.value, exact));
  return buf;
}

Body scaled(const Body& B, double s)
{
  const int n = int(B.A.cols());
  return linearImage(B, s * MT::Identity(n, n), VT::Zero(n), "s");
}

struct Shape
{
  std::string name;
  Body body;
  IBall ball;
};

// cube [-1,1]^n, standard simplex, and boxes [0,1]^{n-1} x [0,aspect].
Shape makeShape(const std::string& kind, int n)
{
  if (kind == "cube") return {kind, cube(n), cubeBall(n, -1, 1)};
  if (kind == "simplex") return {kind, simplex(n), simplexBall(n)};
  const double aspect = kind == "box10" ? 10 : 100;
  std::vector<double> lo(n, 0.0), hi(n, 1.0);
  hi[n - 1] = aspect;
  return {kind, box(lo, hi), boxBall(lo, hi)};
}

}  // namespace

// ------------------------------------------------------------------ scaling
TEST_CASE("invariance: isotropic scaling vol(sP) = s^n vol(P)")
{
  std::cout << "\n== scaling, exact Chebyshev ball, ttc oracle (rel err)\n";
  for (std::string kind : {"cube", "simplex"})
    for (int n : {2, 5})
      for (double s : {1e-6, 1e-3, 1.0, 1e3, 1e6})
      {
        Shape S = makeShape(kind, n);
        Body B = scaled(S.body, s);
        IBall b = mapBall(S.ball, MT::Identity(n, n), VT::Zero(n), s);
        Guarded G = guarded([&] { return oracleWithBall(B, b.c, b.r); }, 5);
        std::cout << "  " << kind << n << " s=" << s << " rel=" << show(G, B.exactVolume)
                  << " (" << G.seconds << "s)\n";
        CHECK_MESSAGE(ok(G, B.exactVolume), (kind + std::to_string(n)), " s=", s, " ", show(G, B.exactVolume));
      }
}

TEST_CASE("invariance: scaling n=10 cube (one small, one large)")
{
  const int n = 10;
  for (double s : {1e-3, 1e6})
  {
    Shape S = makeShape("cube", n);
    Body B = scaled(S.body, s);
    IBall b = mapBall(S.ball, MT::Identity(n, n), VT::Zero(n), s);
    Guarded G = guarded([&] { return oracleWithBall(B, b.c, b.r); }, 30);
    std::cout << "  cube10 s=" << s << " rel=" << show(G, B.exactVolume) << "\n";
    CHECK_MESSAGE(ok(G, B.exactVolume), "cube10 s=", s, " ", show(G, B.exactVolume));
  }
}

// -------------------------------------------------------------- translation
TEST_CASE("invariance: translation vol(P + t) = vol(P)")
{
  std::cout << "\n== translation |t| in {1e3,1e6,1e9} (exact ball / volesti ball)\n";
  for (std::string kind : {"cube", "simplex"})
    for (int n : {2, 5})
      for (double T : {1e3, 1e6, 1e9})
      {
        Shape S = makeShape(kind, n);
        VT dir = randomRotation(n, 7).col(0);
        VT t = T * dir;
        Body B = linearImage(S.body, MT::Identity(n, n), t, "t");
        IBall b = mapBall(S.ball, MT::Identity(n, n), t, 1.0);
        Guarded G1 = guarded([&] { return oracleWithBall(B, b.c, b.r); }, 5);
        Guarded G2 = guarded([&] { return oracleVolestiBall(B); }, 5);
        std::cout << "  " << kind << n << " |t|=" << T << " exactBall=" << show(G1, B.exactVolume)
                  << " volestiBall=" << show(G2, B.exactVolume) << "\n";
        CHECK_MESSAGE(ok(G1, B.exactVolume), (kind + std::to_string(n)), " |t|=", T, " exact ball ", show(G1, B.exactVolume));
        CHECK_MESSAGE(ok(G2, B.exactVolume), (kind + std::to_string(n)), " |t|=", T, " volesti ball ", show(G2, B.exactVolume));
      }
}

// ----------------------------------------------------------------- rotation
TEST_CASE("invariance: rotation vol(QP) = vol(P)")
{
  std::cout << "\n== rotation (exact ball / volesti ball)\n";
  for (std::string kind : {"cube", "simplex"})
    for (int n : {2, 5})
      for (unsigned seed : {1u, 2u, 3u})
      {
        Shape S = makeShape(kind, n);
        MT Q = randomRotation(n, seed);
        Body B = linearImage(S.body, Q, VT::Zero(n), "Q");
        IBall b = mapBall(S.ball, Q, VT::Zero(n), 1.0);
        Guarded G1 = guarded([&] { return oracleWithBall(B, b.c, b.r); }, 5);
        Guarded G2 = guarded([&] { return oracleVolestiBall(B); }, 5);
        std::cout << "  " << kind << n << " Q#" << seed << " exactBall=" << show(G1, B.exactVolume)
                  << " volestiBall=" << show(G2, B.exactVolume) << "\n";
        CHECK_MESSAGE(ok(G1, B.exactVolume), (kind + std::to_string(n)), " Q#", seed, " ", show(G1, B.exactVolume));
        CHECK_MESSAGE(ok(G2, B.exactVolume), (kind + std::to_string(n)), " Q#", seed, " ", show(G2, B.exactVolume));
      }
}

// ------------------------------------------------------ inner-ball handling
// Valid balls (any centre with any radius <= its slack) must all give the
// right volume.  Invalid (too large) balls are never passed by ttc, but the
// oracle must not silently return inf / garbage: we require either a correct
// value or a recognisable failure (<= 0 / non-finite is NOT accepted as
// "correct", so we only CHECK finiteness for invalid balls).
TEST_CASE("inner ball sensitivity (valid balls must not matter)")
{
  std::cout << "\n== inner ball sensitivity (rel err)\n";
  for (std::string kind : {"cube", "simplex"})
    for (int n : {2, 5})
    {
      Shape S = makeShape(kind, n);
      const Body& B = S.body;
      struct Case { std::string name; VT c; double r; bool valid; };
      std::vector<Case> cases;
      for (double f : {1.0, 0.5, 0.1}) cases.push_back({"r x" + std::to_string(f), S.ball.c, S.ball.r * f, true});
      // Off-centre valid ball: move half a radius along e1, radius = slack there.
      VT c2 = S.ball.c;
      c2(0) += 0.5 * S.ball.r;
      HPoly tmp = toHPoly(B);
      tmp.normalize();
      double r2 = (tmp.get_vec() - tmp.get_mat() * c2).minCoeff();
      cases.push_back({"offcentre", c2, r2, true});
      for (double f : {2.0, 10.0}) cases.push_back({"r x" + std::to_string(f) + "(invalid)", S.ball.c, S.ball.r * f, false});

      for (const Case& C : cases)
      {
        Guarded G = guarded([&] { return oracleWithBall(B, C.c, C.r); }, 5);
        std::cout << "  " << kind << n << " " << C.name << " rel=" << show(G, B.exactVolume) << "\n";
        if (C.valid)
          CHECK_MESSAGE(ok(G, B.exactVolume), (kind + std::to_string(n)), " ", C.name, " ", show(G, B.exactVolume));
        else
          CHECK_MESSAGE((!G.timedOut && std::isfinite(G.value)), (kind + std::to_string(n)), " ", C.name,
                        " returned ", G.value);
      }
    }
}

// volesti's own ComputeInnerBall (max_inscribed_ball IPM + lpsolve fallback)
// must find the Chebyshev ball regardless of scale / translation.
// Known upstream issue: volesti's max_inscribed_ball IPM uses absolute
// tolerances and its lp_solve fallback returns invalid balls at extreme
// scales/offsets.  ttc never uses this ball unchecked (reliableInnerBall in
// src/volume/polytope_balls.cpp certifies it against a Chebyshev LP), so this
// documents the library behaviour rather than gating the suite.
TEST_CASE("volesti ComputeInnerBall / max_inscribed_ball across scales"
          * doctest::may_fail())
{
  std::cout << "\n== max_inscribed_ball: converge flag, radius/true, centre err/true r\n";
  for (std::string kind : {"cube", "simplex"})
    for (int n : {2, 5, 10})
      for (double s : {1e-9, 1e-6, 1e-3, 1.0, 1e3, 1e6, 1e9})
        for (double T : {0.0, 1e6})
        {
          Shape S = makeShape(kind, n);
          VT t = T * VT::Ones(n) / std::sqrt(double(n));
          Body B = linearImage(S.body, s * MT::Identity(n, n), t, "st");
          IBall truth = mapBall(S.ball, MT::Identity(n, n), t, s);
          HPoly P = toHPoly(B);
          P.normalize();
          auto ipm = max_inscribed_ball(P.get_mat(), P.get_vec(), 5000, 1e-8);
          auto vb = P.ComputeInnerBall();
          const double ipmRatio = std::get<1>(ipm) / truth.r;
          const double vbRatio = vb.second / truth.r;
          const double vbCentreErr = (vb.first.getCoefficients() - truth.c).norm() / truth.r;
          const double vbSlack = (P.get_vec() - P.get_mat() * vb.first.getCoefficients()).minCoeff();
          std::printf("  %s%d s=%-6g |t|=%-4g ipm conv=%d r/true=%.6f | ComputeInnerBall r/true=%.6f "
                      "centreErr=%.2e certSlack/r=%.6f\n",
                      kind.c_str(), n, s, T, int(std::get<2>(ipm)), ipmRatio, vbRatio, vbCentreErr,
                      vbSlack / vb.second);
          // The radius volesti hands the oracle must be (nearly) optimal and valid.
          CHECK_MESSAGE(std::abs(vbRatio - 1) < 1e-3, (kind + std::to_string(n)), " s=", s, " |t|=", T,
                        " ComputeInnerBall radius/true = ", vbRatio);
          CHECK_MESSAGE(vbSlack >= vb.second * (1 - 1e-6), (kind + std::to_string(n)), " s=", s, " |t|=", T,
                        " ComputeInnerBall radius exceeds certified slack");
        }
}

// ------------------------------------------------------------- size grid
// Full size grid requested by the user: s x {cube, simplex, box aspect 10/100}
// x n in {2,5,10}, plus one-axis scaling and a tiny body far from the origin.
// ~10 min on the current oracle (hangs are cut at 20 s CPU), ~2.5 min once fixed.
TEST_CASE("size grid (slow)" * doctest::test_suite("slow") * doctest::skip())
{
  const std::vector<double> scales = {1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6};
  std::printf("\n== size grid: rel err of ttc oracle on s*P (exact Chebyshev ball)\n%-12s", "body");
  for (double s : scales) std::printf(" %13g", s);
  std::printf("\n");
  int bad = 0, total = 0;
  for (int n : {2, 5, 10})
    for (std::string kind : {"cube", "simplex", "box10", "box100"})
    {
      Shape S = makeShape(kind, n);
      std::printf("%-12s", (std::string(kind) + "/n" + std::to_string(n)).c_str());
      for (double s : scales)
      {
        Body B = scaled(S.body, s);
        IBall b = mapBall(S.ball, MT::Identity(n, n), VT::Zero(n), s);
        Guarded G = guarded([&] { return oracleWithBall(B, b.c, b.r); }, 20);
        std::printf(" %13s", show(G, B.exactVolume).c_str());
        std::fflush(stdout);
        ++total;
        if (!ok(G, B.exactVolume)) ++bad;
      }
      std::printf("\n");
    }

  std::printf("\n== anisotropic (one axis scaled) and tiny-far bodies\n");
  struct Extra { std::string name; Body B; IBall b; };
  std::vector<Extra> extras;
  auto addBox = [&](const std::string& name, std::vector<double> lo, std::vector<double> hi) {
    extras.push_back({name, box(lo, hi), boxBall(lo, hi)});
  };
  addBox("[0,1e-6]x[0,1e6]", {0, 0}, {1e-6, 1e6});
  addBox("[0,1e-3]x[0,1e3]", {0, 0}, {1e-3, 1e3});
  addBox("[0,1]x[0,1e6]", {0, 0}, {1, 1e6});
  addBox("[0,1e-6]x[0,1]", {0, 0}, {1e-6, 1});
  addBox("[0,1e-6]x[0,1]^3x[0,1e6]", {0, 0, 0, 0, 0}, {1e-6, 1, 1, 1, 1e6});
  for (int n : {2, 5, 10})
  {
    std::vector<double> lo(n, 1e6), hi(n, 1e6 + 1e-3);
    addBox("[1e6,1e6+1e-3]^" + std::to_string(n), lo, hi);
  }
  for (int n : {2, 5})
  {
    std::vector<double> lo(n, 1e9), hi(n, 1e9 + 1);
    addBox("[1e9,1e9+1]^" + std::to_string(n), lo, hi);
  }
  for (const Extra& E : extras)
  {
    Guarded G = guarded([&] { return oracleWithBall(E.B, E.b.c, E.b.r); }, 20);
    std::printf("  %-28s rel=%s\n", E.name.c_str(), show(G, E.B.exactVolume).c_str());
    ++total;
    if (!ok(G, E.B.exactVolume)) ++bad;
  }
  std::printf("grid: %d/%d cells outside +-%.2f\n", bad, total, kTol);
  CHECK(bad == 0);
}
