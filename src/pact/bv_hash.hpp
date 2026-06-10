#pragma once

#include <cvc5/cvc5.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>
#include <vector>

#include "pact/hash_constraint.hpp"

// Hash family used over the bit-vector projection variables. The default XOR
// family is ApproxMC-style parity hashing; the Prime and Lemire families are the
// word-level 2-universal hashes ported from cvc5's SmtApproxMc (a*x mod p == c
// and a multiply-shift hash). The word-level families are only meaningful when
// the projection variables are bit-vectors, so they are gated on that in main.
enum class BvHashMode
{
  Xor,     // parity over a random subset of projection bits (default)
  Prime,   // sum_i a_i * x_i  (mod p) == c, p the smallest prime > 2^slice
  Lemire,  // high `slice` bits of sum_i a_i * x_i == c  (multiply-shift)
  PrimeGj  // same constraint as Prime, but emitted as a mod-p row for cvc5's
           // Z_p Gauss-Jordan propagator (assertModpClause) instead of as a
           // bit-vector arithmetic term, avoiding the multiply/urem bit-blast.
};

// Generates random hash constraints over bit-vector projection variables.
class BvHash
{
 public:
  explicit BvHash(cvc5::Solver& solver);

  void setVariables(std::vector<cvc5::Term> vars);

  std::optional<XorClause> randomClause(std::mt19937& rng) const;

  // Number of cells a single hash of the active family splits the space into:
  // 2 for XOR, the prime modulus for Prime, 2^slice for Lemire. Pact scales each
  // per-measurement estimate by this value raised to the number of active hashes
  // (the XOR machinery's 2^hashCount generalised to an arbitrary base).
  double perHashMultiplier() const;

 private:
  // The original ApproxMC-style parity over projection bits.
  std::optional<XorClause> randomXorClause(std::mt19937& rng) const;

  // Word-level hashes ported from cvc5's SmtApproxMc. They return a full
  // bit-vector equality term; the caller wraps it in an XorClause with no native
  // XOR literals so HashConstraint always asserts the (arithmetic) fallback.
  cvc5::Term generatePrimeHash(std::mt19937& rng) const;
  cvc5::Term generateLemireHash(std::mt19937& rng) const;

  // PrimeGj: build the same prime hash as generatePrimeHash but return it as a
  // mod-p row -- per-projection-bit literals + Z_p weights, the modulus p and
  // the rhs c -- together with the equivalent bit-vector arithmetic term in
  // `fallback` (used only for the model-reuse check, never asserted).
  XorClause generatePrimeRow(std::mt19937& rng) const;

  // Minimum bit-width (mirroring SmtApproxMc::getMinBW) needed for the Prime
  // hash's accumulator so the modular sum does not overflow.
  std::uint32_t getMinBitwidth() const;

  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_vars;
  std::size_t d_totalBits;
  BvHashMode d_mode;
  std::uint32_t d_maxBitwidth;  // widest projection bit-vector
  std::uint32_t d_sliceSize;    // word size used by the Prime/Lemire hashes
};
