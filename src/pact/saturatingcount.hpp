#pragma once

#include <cvc5/cvc5.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "pact/hash_constraint.hpp"

// SaturatingCounter evaluates the number of satisfying assignments for the
// current solver assertions with respect to a set of projection variables. The
// counter stops once the provided threshold is reached and reports
// std::nullopt in that case.
class SaturatingCounter
{
 public:
  explicit SaturatingCounter(
      cvc5::Solver& solver,
      const std::vector<cvc5::Term>* projectionVars = nullptr);

  void setProjectionVars(const std::vector<cvc5::Term>* projectionVars);

  // Rebind the counter to a different solver instance. Used to swap in a fresh
  // solver between rounds so accumulated incremental SAT state does not pile up.
  void setSolver(cvc5::Solver& solver) { d_solver = &solver; }

  void setUseNativeXor(bool useNativeXor) { d_useNativeXor = useNativeXor; }

  void resetStatistics();

  // Drops the incremental model/constraint cache. Call this after the solver's
  // assertion stack has been reset out from under the counter.
  void resetCache();

  std::uint64_t getSmtCallCount() const { return d_smtCalls; }

  std::optional<std::size_t> count(
      const std::vector<HashConstraint>& additionalConstraints,
      std::size_t threshold);

 private:
  cvc5::Solver* d_solver;
  const std::vector<cvc5::Term>* d_projectionVars;
  std::vector<HashConstraint> d_cachedConstraints;
  std::vector<std::vector<cvc5::Term>> d_cachedModels;
  std::uint64_t d_smtCalls = 0;
  bool d_useNativeXor = false;
};

