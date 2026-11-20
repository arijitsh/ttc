#include "volume/polytope_balls.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <lp_lib.h>

#include "logger.hpp"
#include "volume/volesti_body.hpp"

namespace ttc
{
namespace
{

using LpPtr = std::unique_ptr<lprec, decltype(&delete_lp)>;

LpPtr makeLp(int columns)
{
  LpPtr lp(make_lp(0, columns), delete_lp);
  if (!lp) throw std::runtime_error("outer-radius LP allocation failed");
  set_verbose(lp.get(), CRITICAL);
  return lp;
}

// ------------------------------------------------------------ outer radius

// max sign*y_j  s.t.  A y <= s, y free, on a fresh lp_solve model.
int solvePrimalExtent(const Eigen::MatrixXd& A, const Eigen::VectorXd& s, int j, double sign,
                      bool scaled, double& value)
{
  const int dim = static_cast<int>(A.cols());
  LpPtr lp = makeLp(dim);
  if (!scaled) set_scaling(lp.get(), SCALE_NONE);
  for (int k = 1; k <= dim; ++k) set_unbounded(lp.get(), k);
  std::vector<double> row(dim + 1, 0.0);
  set_add_rowmode(lp.get(), TRUE);
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    for (int k = 0; k < dim; ++k) row[k + 1] = A(i, k);
    if (!add_constraint(lp.get(), row.data(), LE, s(i)))
      throw std::runtime_error("outer-radius LP constraint failed");
  }
  set_add_rowmode(lp.get(), FALSE);
  std::fill(row.begin(), row.end(), 0.0);
  row[j + 1] = sign;
  if (!set_obj_fn(lp.get(), row.data())) throw std::runtime_error("outer-radius LP objective failed");
  set_maxim(lp.get());
  const int status = solve(lp.get());
  value = get_objective(lp.get());
  return status;
}

// The same extent via the LP dual  min s.lambda  s.t.  A^T lambda = sign*e_j,
// lambda >= 0.  It has no free variables, and any dual-feasible lambda
// certifies sign*y_j <= s.lambda on the whole polytope.
int solveDualExtent(const Eigen::MatrixXd& A, const Eigen::VectorXd& s, int j, double sign,
                    double& value)
{
  const int dim = static_cast<int>(A.cols());
  const int rows = static_cast<int>(A.rows());
  LpPtr lp = makeLp(rows);
  std::vector<double> row(rows + 1, 0.0);
  set_add_rowmode(lp.get(), TRUE);
  for (int k = 0; k < dim; ++k)
  {
    for (int i = 0; i < rows; ++i) row[i + 1] = A(i, k);
    if (!add_constraint(lp.get(), row.data(), EQ, k == j ? sign : 0.0))
      throw std::runtime_error("outer-radius LP constraint failed");
  }
  set_add_rowmode(lp.get(), FALSE);
  for (int i = 0; i < rows; ++i) row[i + 1] = s(i);
  if (!set_obj_fn(lp.get(), row.data())) throw std::runtime_error("outer-radius LP objective failed");
  set_minim(lp.get());
  const int status = solve(lp.get());
  if (status != OPTIMAL) return status;
  std::vector<double> raw(rows);
  get_variables(lp.get(), raw.data());
  const Eigen::VectorXd lambda = Eigen::Map<Eigen::VectorXd>(raw.data(), rows).cwiseMax(0.0);
  Eigen::VectorXd target = Eigen::VectorXd::Zero(dim);
  target(j) = sign;
  if (!((A.transpose() * lambda - target).norm() <= 1e-7)) return NUMFAILURE;
  value = s.dot(lambda);
  return status;
}

// -------------------------------------------------------------- inner ball

// Chebyshev ball  max r  s.t.  a_i.x + r <= b_i  (unit-norm rows, x free,
// r >= 0) on a fresh lp_solve model; unscaled first, then default scaling.
bool solveChebyshevBall(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                        Eigen::VectorXd& center, double& radius)
{
  const int dim = static_cast<int>(A.cols());
  for (const bool scaled : {false, true})
  {
    LpPtr lp = makeLp(dim + 1);
    if (!scaled) set_scaling(lp.get(), SCALE_NONE);
    for (int k = 1; k <= dim; ++k) set_unbounded(lp.get(), k);
    std::vector<double> row(dim + 2, 0.0);
    row[dim + 1] = 1.0;
    set_add_rowmode(lp.get(), TRUE);
    for (Eigen::Index i = 0; i < A.rows(); ++i)
    {
      for (int k = 0; k < dim; ++k) row[k + 1] = A(i, k);
      if (!add_constraint(lp.get(), row.data(), LE, b(i)))
        throw std::runtime_error("Chebyshev-ball LP constraint failed");
    }
    set_add_rowmode(lp.get(), FALSE);
    std::fill(row.begin(), row.end(), 0.0);
    row[dim + 1] = 1.0;
    if (!set_obj_fn(lp.get(), row.data()))
      throw std::runtime_error("Chebyshev-ball LP objective failed");
    set_maxim(lp.get());
    if (solve(lp.get()) != OPTIMAL) continue;
    std::vector<double> values(dim + 1);
    get_variables(lp.get(), values.data());
    center = Eigen::Map<Eigen::VectorXd>(values.data(), dim);
    radius = values[dim];
    return true;
  }
  return false;
}

}  // namespace

void normalizeRows(Eigen::MatrixXd& A, Eigen::VectorXd& b)
{
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    const double norm = A.row(i).norm();
    if (norm > 0.0)
    {
      A.row(i) /= norm;
      b(i) /= norm;
    }
  }
}

// Coordinate extrema bound P by a box, hence by a ball about `center`.  The
// LPs are posed in y = x - center with unit-norm rows, so y = 0 is strictly
// feasible, and every extent gets a fresh model: re-optimising one lp_solve
// model under a new objective reported spurious UNBOUNDED on ill-conditioned
// zonotopes (max y_j "unbounded" while min y_j = -229 on a centrally
// symmetric parallelotope).  Rows are already equilibrated, so the first
// solve is unscaled; a non-optimal status is retried with lp_solve's default
// scaling and then through the dual before the polytope is rejected.
double outerRadius(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                   const Eigen::VectorXd& center)
{
  const int dim = static_cast<int>(A.cols());
  Eigen::MatrixXd unitA = A;
  Eigen::VectorXd slack = b - A * center;
  normalizeRows(unitA, slack);
  double radius = 0.0;
  for (int j = 0; j < dim; ++j)
  {
    double extent = 0.0;
    for (const double sign : {1.0, -1.0})
    {
      double value = 0.0;
      int primal = solvePrimalExtent(unitA, slack, j, sign, false, value);
      if (primal != OPTIMAL) primal = solvePrimalExtent(unitA, slack, j, sign, true, value);
      if (primal != OPTIMAL)
      {
        const int dual = solveDualExtent(unitA, slack, j, sign, value);
        if (dual == INFEASIBLE && primal == UNBOUNDED) throw std::runtime_error("polytope is unbounded");
        if (dual != OPTIMAL)
          throw std::runtime_error("outer-radius LP failed (lp_solve status " + std::to_string(primal)
                                   + "/" + std::to_string(dual) + ")");
      }
      extent = std::max(extent, std::abs(value));
    }
    radius = std::hypot(radius, extent);
  }
  if (!std::isfinite(radius) || radius <= 0.0) throw std::runtime_error("invalid polytope outer radius");
  return radius;
}

namespace detail
{

// volesti's ComputeInnerBall returns the max_inscribed_ball IPM iterate
// without checking that the IPM converged (it did not on 116 of 393
// syntheticLRA1131 polytopes, and stalled at a ball 11x too small on a
// dim-27 simplex), and its lp_solve fallback trusts any OPTIMAL status.  Yet
// the oracle's walk step 4r/sqrt(a n), the clipping bounds and GetPrecision
// all scale with the radius.  An overestimated radius stalls the Gaussian walk
// at the centre, the annealing schedule collapses to a = 0, and the oracle
// returns (pi/0)^(n/2) = inf; an underestimated one slows the walk and biases
// the estimate low.  So re-derive the radius from the slacks at each
// candidate centre, compare with a fresh Chebyshev LP, and keep the larger
// certified ball.  A valid volesti ball within 0.1% of the LP optimum is
// returned untouched, so well-behaved inputs give bit-identical results.
std::pair<Point, double> reliableInnerBall(HPoly& P)
{
  constexpr double kRelTol = 1e-6;      // claimed radius vs. certified slack
  constexpr double kOptimalTol = 1e-3;  // valid but clearly suboptimal ball
  const auto ball = P.ComputeInnerBall();  // also row-normalizes P
  const Eigen::MatrixXd A = P.get_mat();
  const Eigen::VectorXd b = P.get_vec();
  const double none = -std::numeric_limits<double>::infinity();
  auto slackRadius = [&](const Eigen::VectorXd& c) {
    return c.size() == A.cols() && c.allFinite() ? (b - A * c).minCoeff() : none;
  };

  const double volestiCert = slackRadius(ball.first.getCoefficients());
  Eigen::VectorXd lpCenter;
  double lpRadius = 0.0;
  double lpCert = none;
  if (solveChebyshevBall(A, b, lpCenter, lpRadius)) lpCert = std::min(lpRadius, slackRadius(lpCenter));

  const bool volestiValid = std::isfinite(ball.second) && ball.second > 0.0
                            && ball.second <= volestiCert * (1.0 + kRelTol);
  if (volestiValid && !(ball.second < lpCert * (1.0 - kOptimalTol))) return ball;

  std::pair<Point, double> fixed{Point(static_cast<unsigned>(A.cols())), -1.0};
  if (lpCert >= volestiCert && lpCert > 0.0)
    fixed = {Point(lpCenter), lpCert};
  else if (volestiCert > 0.0)
    fixed = {ball.first, volestiCert};
  Log(2) << "inner ball: volesti radius " << ball.second << " (certified " << volestiCert
         << "), Chebyshev LP radius " << lpCert << "; using " << fixed.second << std::endl;
  if (fixed.second > 0.0) P.set_InnerBall(fixed);
  return fixed;
}

HPoly makeBody(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const InnerBall& ball)
{
  HPoly P(static_cast<unsigned>(A.cols()), A, b);
  P.normalize();
  if (ball.radius > 0.0) P.set_InnerBall({Point(ball.center), ball.radius});
  return P;
}

}  // namespace detail

InnerBall certifiedInnerBall(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
{
  detail::HPoly P(static_cast<unsigned>(A.cols()), A, b);
  const auto ball = detail::reliableInnerBall(P);
  return {ball.first.getCoefficients(), ball.second};
}

}  // namespace ttc
