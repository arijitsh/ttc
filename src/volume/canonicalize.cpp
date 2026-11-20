#include "volume/canonicalize.hpp"

#include <cmath>
#include <cstdlib>
#include <utility>

#include "volume/polytope_balls.hpp"

// cddlib (double precision, not the GMP build).  Its set_int clashes with
// lp_solve's, hence the rename.
#define set_int cdd_set_int
extern "C" {
#include <cddlib/setoper.h>
#include <cddlib/cdd.h>
}
#undef set_int

namespace ttc
{
namespace
{

// dd_MatrixCanonicalize on {x : A x <= b}, as the Python prototype's
// pycddlib canonicalize.  Same contract as canonicalizePolytope.
bool cddCanonicalize(Eigen::MatrixXd& A, Eigen::VectorXd& b)
{
  const long m = A.rows();
  const long n = A.cols();
  if (m == 0) return true;

  dd_set_global_constants();
  // A cdd row (g_0, g_1, ..., g_n) means g_0 + g.x >= 0, so A_i x <= b_i is
  // the row (b_i, -A_i1, ..., -A_in).
  dd_MatrixPtr M = dd_CreateMatrix(m, n + 1);
  M->representation = dd_Inequality;
  for (long i = 0; i < m; ++i)
  {
    dd_set_d(M->matrix[i][0], b(i));
    for (long j = 0; j < n; ++j) dd_set_d(M->matrix[i][j + 1], -A(i, j));
  }

  dd_rowset implicitLinearity = nullptr;
  dd_rowset redundant = nullptr;
  dd_rowindex newPosition = nullptr;
  dd_ErrorType err = dd_NoError;
  dd_MatrixCanonicalize(&M, &implicitLinearity, &redundant, &newPosition, &err);

  bool keep = true;
  if (err == dd_NoError && M != nullptr)
  {
    if (set_card(M->linset) > 0)
      keep = false;  // an implicit equality: flat
    else
    {
      const long rows = M->rowsize;
      Eigen::MatrixXd newA(rows, n);
      Eigen::VectorXd newB(rows);
      for (long i = 0; i < rows; ++i)
      {
        newB(i) = dd_get_d(M->matrix[i][0]);
        for (long j = 0; j < n; ++j) newA(i, j) = -dd_get_d(M->matrix[i][j + 1]);
      }
      A = std::move(newA);
      b = std::move(newB);
    }
  }

  if (implicitLinearity != nullptr) set_free(implicitLinearity);
  if (redundant != nullptr) set_free(redundant);
  if (newPosition != nullptr) free(newPosition);
  if (M != nullptr) dd_FreeMatrix(M);
  dd_free_global_constants();
  return keep;
}

}  // namespace

bool canonicalizePolytope(Eigen::MatrixXd& A, Eigen::VectorXd& b)
{
  if (A.rows() == 0) return true;

  // cddlib decides redundancy and implicit equalities with absolute
  // tolerances (~1e-7), so a full-dimensional polytope of width ~1e-7 used to
  // be reported flat.  Canonicalize in the frame x = c + r y of the certified
  // Chebyshev ball (inner radius 1, unit-norm rows) and map back, so that the
  // verdicts do not depend on the formula's units.
  const InnerBall frame = certifiedInnerBall(A, b);
  if (!(frame.radius > 0.0) || !std::isfinite(frame.radius))
    return frame.radius > 0.0;  // flat or empty: zero volume; unbounded: keep raw
  Eigen::MatrixXd Ay = A * frame.radius;
  Eigen::VectorXd by = b - A * frame.center;
  normalizeRows(Ay, by);
  if (!cddCanonicalize(Ay, by)) return false;
  // alpha.y <= beta with y = (x - c)/r  <=>  (alpha/r).x <= beta + (alpha/r).c
  A = Ay / frame.radius;
  b = by + A * frame.center;
  return true;
}

}  // namespace ttc
