#pragma once

// Single-polytope geometry and volume for the LRA engine (lra_volume.cpp) and
// the standalone tests in tests/volesti/.  Plain functions on Eigen matrices;
// volesti and lp_solve stay inside polytope_volume.cpp.
//
// Cousins & Vempala's practical volume algorithm (Math. Prog. Comp. 2016,
// Fig. 2) starts with `T = Round(K)` and returns |det T| * vol(T^-1 K): every
// volesti walk steps in units of the inner radius, so on a skinny body
// (R/r >> sqrt n) it never crosses the body within its step budget and the
// estimate collapses (box [0,2]x[0,1000]: 7.78 instead of 2000).  None of
// volesti's volume routines rounds internally, so RoundedPolytope applies the
// max-inscribed-ellipsoid rounding of Cousins-Vempala's preprocess.m
// (volesti's inscribed_ellipsoid_rounding) first.  See
// text-for-agents/ttc-volume-vs-theory.md.

#include <memory>
#include <string>

#include <Eigen/Dense>

#include "volume/volume_oracle_kind.hpp"

namespace ttc
{

struct InnerBall
{
  Eigen::VectorXd center;
  double radius = -1.0;  // <= 0: no full-dimensional ball (empty / flat polytope)
};

// Certified Chebyshev ball of {x : A x <= b}.  volesti's ComputeInnerBall (an
// IPM whose convergence flag it ignores, with an lp_solve fallback) is
// re-checked against the slacks at its centre and against a fresh Chebyshev
// LP, keeping the larger certified ball; a valid volesti ball within 0.1% of
// the LP optimum is returned untouched.
InnerBall certifiedInnerBall(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

// Radius of a ball about `center` containing {x : A x <= b}, from coordinate
// extent LPs (each on a fresh lp_solve model, retried scaled and via the
// dual).  Throws if the polytope is unbounded or the LPs fail.
double outerRadius(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                   const Eigen::VectorXd& center);

// A polytope prepared for the volume oracle: row-normalised, holding its
// certified inner ball, and (unless rounding is off or fails) rounded so that
// its maximum inscribed ellipsoid is the unit ball.  vol(P) = det() *
// vol(rounded P); volume() already includes the det() factor.
class RoundedPolytope
{
 public:
  // `ball` must be certifiedInnerBall(A, b).
  RoundedPolytope(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const InnerBall& ball,
                  bool round);
  ~RoundedPolytope();
  RoundedPolytope(RoundedPolytope&&) noexcept;
  RoundedPolytope& operator=(RoundedPolytope&&) noexcept;

  bool rounded() const;
  double det() const;                          // 1 when not rounded
  const std::string& roundingFailure() const;  // why rounding was skipped
  double innerRadius() const;                  // of the (rounded) body

  // The oracle on volesti's deterministic stream (fixed seed 3), as ttc runs it.
  double volume(VolumeOracle oracle, double error, unsigned walkLength) const;
  // The oracle on an explicitly seeded stream (retries, statistics).
  double volumeWithSeed(VolumeOracle oracle, double error, unsigned walkLength,
                        unsigned seed) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace ttc
