#pragma once

// The polytope P(c) of a cube c (manuscript, Implementation, "Constructing the
// Polytopes"): the conjunction of its literals as {x : A x <= b}.
//
//   l1 <= l2, l1 < l2     one row  (l1 - l2 <= 0; strictness changes no volume)
//   l1 >= l2, l1 > l2     one row  (l2 - l1 <= 0)
//   l1 = l2               two rows
//   not (l1 = l2)         no row   (it removes a hyperplane, a null set)
//   not (l1 <= l2) etc.   the flipped comparison
//
// l1 and l2 are linear terms over the real variables; constants may be
// rational, products need all but one factor constant, and divisions a
// constant divisor.

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <cvc5/cvc5.h>

#include "volume/cube.hpp"

namespace ttc
{

// The coordinates of R^n: the distinct real variables of F in first-occurrence
// order (non-real terms are ignored).  Column j of A is variable j.
class RealVariables
{
 public:
  explicit RealVariables(const std::vector<cvc5::Term>& variables);

  std::size_t size() const { return d_terms.size(); }
  const std::vector<std::string>& names() const { return d_names; }
  // Column of `term`, or -1 if it is not one of the variables.
  long column(const cvc5::Term& term) const;

 private:
  std::vector<cvc5::Term> d_terms;
  std::vector<std::string> d_names;
  std::unordered_map<cvc5::Term, std::size_t> d_column;
};

struct LinearConstraints
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  bool empty = false;  // a literal is a false comparison of constants: P(c) = {}
};

// P(c) as {x : A x <= b}, one row per comparison in the order of c's literals.
// Throws std::runtime_error on a literal that is not a linear comparison and on
// a cube without constraints (P(c) = R^n has no finite volume).
LinearConstraints cubeConstraints(const Cube& cube, const RealVariables& variables);

}  // namespace ttc
