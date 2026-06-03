#pragma once

#include <cvc5/cvc5.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <random>
#include <vector>

#include "pact/bool_hash.hpp"
#include "pact/bv_hash.hpp"
#include "pact/hash_constraint.hpp"
#include "pact/saturatingcount.hpp"

// Implements the core PACT approximate model counting loop based on random
// parity constraints over Boolean and bit-vector projection variables.
class Pact
{
 public:
  Pact(cvc5::Solver& solver,
       const std::vector<cvc5::Term>& projectionVars,
       std::uint64_t seed = 42,
       bool useNativeXor = false,
       bool bvPact = false);

  std::uint64_t count();

  std::uint64_t getSmtCallCount() const { return d_counter.getSmtCallCount(); }

 private:
  struct Parameters
  {
    std::size_t threshold;
    std::size_t iterations;
  };

  // Outcome of a single galloping measurement: the model count found within the
  // pivot at the winning number of hashes, plus that number of hashes. The
  // estimate for the measurement is cellCount << hashCount.
  struct MeasurementResult
  {
    std::size_t hashCount;
    std::size_t cellCount;
  };

  // Called for each saturating-count evaluation inside a measurement so the
  // caller can emit progress. Arguments: number of hashes evaluated, the count
  // result (nullopt when the pivot was saturated), and the next number of
  // hashes the galloping search will jump to.
  using ReportFn = std::function<void(
      std::size_t, const std::optional<std::size_t>&, std::int64_t)>;

  Parameters getParameters() const;
  std::uint64_t median(std::vector<std::uint64_t>& values) const;
  HashConstraint generateHashConstraint();

  // One galloping (exponential-then-binary) search over the number of hashes,
  // mirroring ApproxMC's one_measurement_count. prevMeasure carries the winning
  // hash count between measurements so later rounds start near the answer.
  MeasurementResult oneMeasurement(std::size_t iter,
                                   std::int64_t& prevMeasure,
                                   std::vector<HashConstraint>& hashPool,
                                   std::size_t threshold,
                                   const ReportFn& report);

  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_projectionVars;
  std::vector<cvc5::Term> d_booleanVars;
  std::vector<cvc5::Term> d_booleanBvVars;
  std::vector<cvc5::Term> d_bitvectorVars;
  std::mt19937 d_rng;
  SaturatingCounter d_counter;
  BoolHash d_boolHash;
  BvHash d_bvHash;
  // Upper bound on the number of useful hashes (total projection bits); the
  // galloping search's hiIndex sentinel.
  std::size_t d_maxHashes = 1;

  // Whether parity hashes are handed to the SAT solver as native XOR clauses
  // (CaDiCaL Gauss-Jordan) instead of CNF. Requires the core SAT solver to be
  // CaDiCaL, which rebuildCountSolver enforces on the counting solver.
  bool d_useNativeXor = false;

  // The base formula (after any BV-PACT substitution), captured once so we can
  // rebuild a fresh counting solver between measurements.
  std::vector<cvc5::Term> d_baseAssertions;

  // A dedicated solver (sharing d_solver's term manager) on which all hashed
  // counting happens. It is rebuilt from scratch each round: cvc5's incremental
  // SAT backend accumulates the bit-blasted hash clauses of every round and
  // never reclaims them across push/pop or resetAssertions, so reusing one
  // solver makes each successive round dramatically slower. Discarding the
  // solver per round keeps the per-round cost flat.
  std::optional<cvc5::Solver> d_countSolver;
  void rebuildCountSolver();
};

