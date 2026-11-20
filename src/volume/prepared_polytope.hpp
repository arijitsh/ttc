#pragma once

// The polytopes P_i = P(c_i) of Algorithm 1, ready for the main loop:
// half-spaces (cdd-canonical) together with the two balls
//   B(z_i, gamma_i) ⊆ P_i ⊆ B(z_i, R_i)
// that GetPrecision (gamma_i), the VolEst clipping (gamma_i, R_i) and the
// sampling walks (z_i) use.  Cubes whose polytope has zero n-dimensional
// volume (empty or flat) contribute nothing to the union and are dropped.

#include <cstddef>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "volume/polytope_balls.hpp"

namespace ttc
{

struct Cube;           // cube.hpp
class RealVariables;   // linear_constraints.hpp

struct PreparedPolytope
{
  Eigen::MatrixXd A;  // P_i = {x : A x <= b}
  Eigen::VectorXd b;
  InnerBall ball;            // certified Chebyshev ball: centre z_i, radius gamma_i
  double outerRadius = 0.0;  // R_i
  long rawFacets = 0;        // rows before cdd (reported only)
};

struct PreparedPolytopes
{
  std::vector<PreparedPolytope> kept;
  std::size_t zeroVolume = 0;  // cubes dropped as empty or flat
  double ballSeconds = 0.0;    // inner-ball and outer-radius LPs
};

struct PrepareOptions
{
  bool cddSimplify = true;     // canonicalize with cddlib (--no-cdd-simp clears it)
  std::string dumpInePrefix;   // if set, write each P_i to <prefix>_cubeN.ine
};

// P(c_i) for every cube, in cube order.  Throws on non-finite coefficients and
// on an unbounded polytope.
PreparedPolytopes preparePolytopes(const std::vector<Cube>& cubes, const RealVariables& variables,
                                   const PrepareOptions& options);

// A bound on |x_j| over every point of every kept polytope.
double maxAbsCoordinate(const std::vector<PreparedPolytope>& P);

}  // namespace ttc
