#pragma once

#include <cvc5/cvc5.h>

#include <cstddef>
#include <optional>
#include <random>
#include <vector>

#include "pact/hash_constraint.hpp"

// Generates random parity constraints over bit-vector projection variables.
class BvHash
{
 public:
  explicit BvHash(cvc5::Solver& solver);

  void setVariables(std::vector<cvc5::Term> vars);

  std::optional<XorClause> randomClause(std::mt19937& rng) const;

 private:
  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_vars;
  std::size_t d_totalBits;
};

