#pragma once

// Line 3 of Algorithm 1, b <- GetPrecision(F_bd, M, eta) with eta = eps/96:
// the lattice 10^-b Z^n is fine enough that every P_i holds a ball of radius
// gamma in lattice units, which makes lattice-point counts track volume
// (lem:nonconv).
//
//   Algorithm 2  GetPrecision(F_bd, M, eta)
//    1  maxratio <- 1, r <- 0
//    2  for c_i in F_bd: r <- r + facets(P(c_i))
//    6  for i = 1 to m
//    7    gamma_hat <- MaxInscribedBall(P^i)
//    8    gamma <- 64 m n / eta * sqrt(ln(8 r / eta))
//    9    maxratio <- max(maxratio, gamma / gamma_hat)
//   13  return ceil(log10 maxratio)
//
// gamma does not depend on i, so maxratio = max(1, gamma / min_i gamma_hat_i).
// r and min gamma_hat are taken over the kept polytopes (zero-volume ones
// have no inscribed ball); m counts every cube.

#include <cstddef>
#include <vector>

#include "volume/prepared_polytope.hpp"

namespace ttc
{

// b for the kept polytopes P (inner radii = MaxInscribedBall) out of m cubes
// in R^n; 0 when P is empty.
int getPrecision(const std::vector<PreparedPolytope>& P, std::size_t n, std::size_t m, double eta);

// Lines 8-13 in closed form, ceil(log10 max(1, gamma / minInnerRadius)) for r
// facets in total, computed in logs so gamma need not fit in a double.
// Throws std::runtime_error when b exceeds the range of int.
int latticePrecision(std::size_t n, std::size_t r, std::size_t m, double minInnerRadius,
                     double eta);

}  // namespace ttc
