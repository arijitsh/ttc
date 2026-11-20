#pragma once

// Internal to src/volume/*.cpp: the volesti H-polytope type that the ball,
// volume and sampling code hand to volesti.  Including this header pulls in
// volesti; the public headers stay free of it.

#include <utility>

#include <Eigen/Dense>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "volume/polytope_balls.hpp"

namespace ttc::detail
{

using Kernel = Cartesian<double>;
using Point = Kernel::Point;
using HPoly = HPolytope<Point>;

// volesti's inner ball of P, certified as described at certifiedInnerBall,
// and stored in P where volesti's walks read it.  Normalizes P's rows.
std::pair<Point, double> reliableInnerBall(HPoly& P);

// {x : A x <= b} with rows normalized and `ball` stored (if it has positive
// radius): the state reliableInnerBall leaves a polytope in.
HPoly makeBody(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const InnerBall& ball);

}  // namespace ttc::detail
