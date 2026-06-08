#pragma once

#include <cvc5/cvc5.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <random>
#include <string>
#include <utility>
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
       bool bvPact = false,
       double epsilon = 0.8,
       double delta = 0.2);

  std::uint64_t count();

  struct ProbeResult
  {
    std::size_t hashCount = 0;
    bool saturated = false;
    std::size_t cellCount = 0;
  };

  struct ProgressEvent
  {
    std::size_t hashCount = 0;
    std::optional<std::size_t> cellCount;
    std::string next;
    std::size_t reuseSat = 0;
    std::size_t reuseChecked = 0;
    std::size_t threshold = 0;
  };

  ProbeResult probeHashCount(std::size_t hashCount);
  double estimateFromMeasurement(std::size_t hashCount,
                                 std::size_t cellCount) const;

  std::uint64_t getSmtCallCount() const { return d_counter.getSmtCallCount(); }
  std::size_t plannedIterations() const { return getParameters().iterations; }
  std::size_t lastWinningHash() const { return d_lastWinningHash; }
  std::size_t lastWinningCellCount() const { return d_lastWinningCellCount; }

  void setIterationOverride(std::optional<std::size_t> iterations)
  {
    d_iterationOverride = iterations;
  }

  void setInitialHashCount(std::size_t hashCount)
  {
    d_initialPrevMeasure = static_cast<std::int64_t>(hashCount);
  }

  void setBaseSaturatedKnown(bool saturated) { d_assumeBaseSaturated = saturated; }

  void setProgressCallback(std::function<void(const ProgressEvent&)> callback)
  {
    d_progressCallback = std::move(callback);
  }

  void setQuiet(bool quiet) { d_quiet = quiet; }

 private:
  struct Parameters
  {
    std::size_t threshold;
    std::size_t iterations;
    bool appmc7;
    double alpha;
    double beta;
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
  // result (nullopt when the pivot was saturated), the next number of hashes
  // the galloping search will jump to, and the reuse_models tally for this
  // evaluation (reuseSat stored solutions still satisfied the active hashes out
  // of reuseChecked examined, so they were reused without an SMT call).
  using ReportFn = std::function<void(std::size_t,
                                      const std::optional<std::size_t>&,
                                      std::int64_t,
                                      std::size_t,
                                      std::size_t)>;

  Parameters getParameters() const;

  // Tail of the binomial used by ApproxMC to bound the error probability of a
  // single measurement: sum_{k=ceil(t/2)}^{t} C(t,k) p^k (1-p)^(t-k). Mirrors
  // Counter::calc_error_bound in approxmc. Used to pick the (odd) number of
  // measurements so the combined under/over-estimation probability <= delta.
  double calcErrorBound(std::size_t t, double p) const;

  // ApproxMC6 rounding (CAV23 "Rounding Meets Approximate Model Counting",
  // Algorithm 5): each measurement's cell count is rounded up to an
  // epsilon-dependent multiple of the pivot before the median is taken. This is
  // what makes the p_L/p_U probabilities -- and hence getParameters()'s
  // iteration count -- give the (epsilon, delta) guarantee.
  double roundCount(double cellCount) const;

  HashConstraint generateHashConstraint();

  // A satisfying projection assignment found during the current measurement,
  // tagged with the number of hashes that were active when it was found.
  // ApproxMC's SavedModel / glob_model: a model found with `hashCount` hashes
  // active also satisfies any prefix of those hashes, so it can be reused (as an
  // already-counted, pre-blocked solution) at every lower hash level without a
  // fresh SMT call -- and at higher levels too whenever it still satisfies the
  // extra parities (checked cheaply by modelSatisfiesHashes).
  struct SavedModel
  {
    std::vector<cvc5::Term> values;  // values in d_projectionVars order
    std::size_t hashCount;
  };

  // True iff `model` (projection-variable values) satisfies every parity in
  // `hashes`. Pure term evaluation (substitute + simplify), no SMT solving --
  // the cheap check that lets ApproxMC reuse stored solutions.
  bool modelSatisfiesHashes(const std::vector<cvc5::Term>& model,
                            const std::vector<HashConstraint>& hashes);

  // One galloping (exponential-then-binary) search over the number of hashes,
  // mirroring ApproxMC's one_measurement_count. prevMeasure carries the winning
  // hash count between measurements so later rounds start near the answer.
  MeasurementResult oneMeasurement(std::size_t iter,
                                   std::int64_t& prevMeasure,
                                   std::vector<HashConstraint>& hashPool,
                                   std::size_t threshold,
                                   const ReportFn& report);

  // ApproxMC7 large-epsilon measurement (AAAI25 Algorithm 4): instead of
  // counting up to a pivot, each level only tests whether the hashed cell is
  // non-empty, and the final cell count is adjusted by alpha/beta.
  MeasurementResult oneMeasurementAppmc7(std::size_t iter,
                                         std::int64_t& prevMeasure,
                                         std::vector<HashConstraint>& hashPool,
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

  // Approximation guarantee. epsilon is the multiplicative tolerance and delta
  // the failure probability; getParameters() derives the pivot/threshold and the
  // number of measurements from them exactly as ApproxMC does.
  double d_epsilon = 0.8;
  double d_delta = 0.2;
  std::optional<std::size_t> d_iterationOverride;
  std::int64_t d_initialPrevMeasure = 0;
  std::size_t d_lastWinningHash = 0;
  std::size_t d_lastWinningCellCount = 0;
  bool d_assumeBaseSaturated = false;
  bool d_quiet = false;
  std::function<void(const ProgressEvent&)> d_progressCallback;

  // ApproxMC's reuse_models optimisation: keep the satisfying assignments found
  // so far in the current measurement and reuse the still-valid ones at the next
  // hash count instead of re-deriving them with the SMT solver. On by default;
  // disable with TTC_NO_REUSE for A/B comparison.
  bool d_reuseModels = true;
  std::vector<SavedModel> d_savedModels;
  // Upper bound on the number of useful hashes (total projection bits); the
  // galloping search's hiIndex sentinel.
  std::size_t d_maxHashes = 1;

  // Whether parity hashes are handed to the SAT solver as native XOR clauses
  // (CaDiCaL Gauss-Jordan) instead of CNF. Requires the core SAT solver to be
  // CaDiCaL, which rebuildCountSolver enforces on the counting solver.
  bool d_useNativeXor = false;

  // --xor-activation literal: hold the whole hash pool on a single solver and
  // toggle hashes with indicator-variable assumptions instead of rebuilding the
  // count solver per galloping level (the rebuild mode). Only meaningful with
  // native XOR (CaDiCaL Gauss-Jordan); ignored otherwise. d_assertedHashes
  // tracks how many pool hashes are already asserted on the current solver.
  bool d_xorActivationLiteral = false;
  std::size_t d_assertedHashes = 0;
  std::size_t d_activationCounter = 0;

  // True once the base (no-hash) count has been shown to saturate the pivot.
  // Level 0 of every measurement re-counts that identical hash-free formula, so
  // each measurement (after the first) would otherwise pay the full base
  // enumeration again -- expensive on large formulas with a small model count
  // (e.g. netrel), where the search hovers at 0/1 hashes. Knowing the result is
  // invariant lets oneMeasurement short-circuit level 0.
  bool d_baseSaturated = false;

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
