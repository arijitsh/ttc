#pragma once

// The two balls of a polytope P = {x : A x <= b} that Algorithm 1 relies on,
//   B(z, gamma) ⊆ P ⊆ B(z, R):
// gamma is MaxInscribedBall of GetPrecision (Alg. 2, line 7) and, with R,
// bounds the VolEst clipping of line 8; z is where the sampling walks start.
// Both come from LPs in double arithmetic (lp_solve).

#include <Eigen/Dense>

namespace ttc
{

struct InnerBall
{
  Eigen::VectorXd center;
  double radius = -1.0;  // <= 0: no full-dimensional ball (empty or flat polytope)
};

// Certified Chebyshev ball of P.  volesti's ComputeInnerBall (an IPM whose
// convergence flag it ignores, with an lp_solve fallback) is re-checked
// against the slacks at its centre and against a fresh Chebyshev LP, and the
// larger certified ball is kept; a valid volesti ball within 0.1% of the LP
// optimum is returned untouched.
InnerBall certifiedInnerBall(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

// Scales each nonzero row of A x <= b (and its b) to unit norm, in place.
void normalizeRows(Eigen::MatrixXd& A, Eigen::VectorXd& b);

// Radius of a ball about `center` that contains P, from the 2n coordinate
// extent LPs (each on a fresh lp_solve model, retried scaled and via the
// dual).  Throws if P is unbounded or the LPs fail.
double outerRadius(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                   const Eigen::VectorXd& center);

}  // namespace ttc
