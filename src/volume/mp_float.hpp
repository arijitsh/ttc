#pragma once

// GMP floats for the sampling side (GmpMode::PointRepr and Full).  Their
// precision is global: MpFloat::default_precision(digits) is set once before
// sampling starts (unionvolume.cpp).

#include <Eigen/Dense>
#include <boost/multiprecision/eigen.hpp>
#include <boost/multiprecision/gmp.hpp>

namespace ttc
{

using MpFloat = boost::multiprecision::mpf_float;
using MpVector = Eigen::Matrix<MpFloat, Eigen::Dynamic, 1>;

}  // namespace ttc
