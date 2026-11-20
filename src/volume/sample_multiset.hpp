#pragma once

// The multiset X of Algorithm 1: the lattice points sampled so far that no
// later polytope has claimed.  Templated on the coordinate type, double
// (--nogmp) or MpFloat (GMP modes), so the membership test A x <= b runs at
// the stored precision without double cancellation near a facet.

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include <Eigen/Dense>

namespace ttc
{

template <typename Coord>
class SampleMultiset
{
 public:
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;

  // Lines 9-11: remove every s in X with s in {x : A x <= b}.  Returns the
  // number removed.  A and b are cast to Coord, so the test runs at X's precision.
  std::size_t removeInside(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
  {
    if (A.rows() == 0) return 0;
    using Mat = Eigen::Matrix<Coord, Eigen::Dynamic, Eigen::Dynamic>;
    const Mat Ac = A.template cast<Coord>();
    const Vec bc = b.template cast<Coord>();
    const Coord zero(0);
    return eraseIf([&](const Vec& s) { return ((Ac * s - bc).array() <= zero).all(); });
  }

  // Lines 13 and 17: remove every element independently with probability 1/2.
  // Draws one uniform from rng per element, in order.
  template <typename URNG>
  std::size_t halve(URNG& rng)
  {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    return eraseIf([&](const Vec&) { return unit(rng) < 0.5; });
  }

  // Line 20.
  template <typename Points>
  void append(const Points& S)
  {
    for (const auto& s : S) add(s);
  }
  void add(const Vec& s) { d_points.push_back(s); }

  std::size_t size() const { return d_points.size(); }
  const std::vector<Vec>& points() const { return d_points; }

 private:
  template <typename Predicate>
  std::size_t eraseIf(Predicate&& remove)
  {
    const std::size_t before = d_points.size();
    d_points.erase(std::remove_if(d_points.begin(), d_points.end(), remove), d_points.end());
    return before - d_points.size();
  }

  std::vector<Vec> d_points;
};

}  // namespace ttc
