#pragma once

// cddlib preprocessing of P(c) (manuscript, Implementation, "Preprocessing"):
// removes redundant inequalities and detects hidden equalities, which make
// the polytope flat (zero n-dimensional volume).

#include <Eigen/Dense>

namespace ttc
{

// Canonicalizes {x : A x <= b} in place (same set, fewer rows).  Returns
// false when the polytope is empty or flat; A and b are then left as they
// were.  On a cddlib error no row is removed and true is returned.
bool canonicalizePolytope(Eigen::MatrixXd& A, Eigen::VectorXd& b);

}  // namespace ttc
