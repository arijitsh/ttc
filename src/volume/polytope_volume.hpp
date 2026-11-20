#pragma once

// The engine of VolEst (Algorithm 1, line 8; vol_est.hpp): VolEsti's
// practical Cousins-Vempala volume algorithms, run on a rounded polytope.
//
// Cousins & Vempala (Math. Prog. Comp. 2016, Fig. 2) start with T = Round(K)
// and return |det T| vol(T^-1 K).  Every volesti walk steps in units of the
// inner radius, so on a skinny body (R/r >> sqrt n) it never crosses the body
// within its step budget and the estimate collapses (box [0,2]x[0,1000]: 7.78
// instead of 2000).  None of volesti's volume routines rounds internally, so
// RoundedPolytope first applies the maximum-inscribed-ellipsoid rounding of
// Cousins-Vempala's preprocess.m (volesti's inscribed_ellipsoid_rounding).

#include <memory>
#include <string>

#include <Eigen/Dense>

#include "volume/polytope_balls.hpp"
#include "volume/volume_oracle_kind.hpp"

namespace ttc
{

// Affine rounding of P = {x : A x <= b}: the maximum inscribed ellipsoid of P
// (iterated, 5 iterations, axis-ratio stop 6) is mapped to the unit ball.  The
// rounded body {y : A y <= b} (fields A, b) satisfies x = shift + T y for x in
// P, and vol(P) = det * vol(rounded).  `ball` must be certifiedInnerBall(A, b).
// On failure applied == false and `failure` says why.
struct PolytopeRounding
{
  bool applied = false;
  std::string failure;
  Eigen::MatrixXd T;
  Eigen::VectorXd shift;
  double det = 1.0;
  Eigen::MatrixXd A;  // rounded body, rows normalized
  Eigen::VectorXd b;
  InnerBall ball;     // certified ball of the rounded body
};
PolytopeRounding roundPolytope(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                               const InnerBall& ball);

// A polytope prepared for the volume oracle: rows normalized, its certified
// inner ball stored, and (unless rounding is off or fails) rounded.
// volume() and volumeWithSeed() already include the det() factor.
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

  // The oracle on volesti's deterministic stream (fixed seed 3).
  double volume(VolumeOracle oracle, double error, unsigned walkLength) const;
  // The oracle on an explicitly seeded stream, as ttc runs it.  Non-finite or
  // non-positive on failure.
  double volumeWithSeed(VolumeOracle oracle, double error, unsigned walkLength,
                        unsigned seed) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace ttc
