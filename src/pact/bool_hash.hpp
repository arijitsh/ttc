#pragma once

#include <cvc5/cvc5.h>

#include <optional>
#include <random>
#include <vector>

#include "pact/hash_constraint.hpp"

// Generates random parity constraints over Boolean projection variables.
class BoolHash
{
 public:
  explicit BoolHash(cvc5::Solver& solver);

  void setVariables(std::vector<cvc5::Term> boolVars,
                    std::vector<cvc5::Term> bvVars);

  std::optional<XorClause> randomClause(std::mt19937& rng) const;

 private:
  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_boolVars;
  std::vector<cvc5::Term> d_bvVars;
};

