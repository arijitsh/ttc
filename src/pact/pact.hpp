#pragma once

#include <cvc5/cvc5.h>

#include <cstdint>
#include <optional>
#include <random>
#include <vector>

#include "pact/bool_hash.hpp"
#include "pact/bv_hash.hpp"
#include "pact/hash_constraint.hpp"
#include "pact/next_index.hpp"
#include "pact/saturatingcount.hpp"

// Implements the core PACT approximate model counting loop based on random
// parity constraints over Boolean and bit-vector projection variables.
class Pact
{
 public:
  Pact(cvc5::Solver& solver,
       const std::vector<cvc5::Term>& projectionVars,
       std::uint64_t seed = 42,
       bool useNativeXor = false);

  std::uint64_t count();

  std::uint64_t getSmtCallCount() const { return d_counter.getSmtCallCount(); }

 private:
  struct Parameters
  {
    std::size_t threshold;
    std::size_t iterations;
  };

  Parameters getParameters() const;
  std::uint64_t median(std::vector<std::uint64_t>& values) const;
  HashConstraint generateHashConstraint();

  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_projectionVars;
  std::vector<cvc5::Term> d_booleanVars;
  std::vector<cvc5::Term> d_booleanBvVars;
  std::vector<cvc5::Term> d_bitvectorVars;
  std::mt19937 d_rng;
  SaturatingCounter d_counter;
  BoolHash d_boolHash;
  BvHash d_bvHash;
  NextIndex d_nextIndex;
};

