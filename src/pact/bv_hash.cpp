#include "pact/bv_hash.hpp"

#include "features.hpp"

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

BvHash::BvHash(cvc5::Solver& solver)
    : d_solver(solver), d_totalBits(0)
{}

void BvHash::setVariables(std::vector<cvc5::Term> vars)
{
  d_vars = std::move(vars);
  d_totalBits = 0;
  for (const cvc5::Term& v : d_vars)
  {
    if (v.getSort().isBitVector())
    {
      d_totalBits += v.getSort().getBitVectorSize();
    }
  }
}

std::optional<XorClause> BvHash::randomClause(std::mt19937& rng) const
{
  if (d_totalBits == 0)
  {
    return std::nullopt;
  }

  auto& tm = ttc::getTermBuilder(d_solver);

  // A proper pairwise-independent hash includes each projection bit in the
  // parity independently with probability 1/2 (ApproxMC-style XOR hashing).
  // Sampling only a handful of bits (as an earlier version did) yields
  // degenerate hashes for formulas with few wide projection variables: once a
  // few XORs are added the surviving models become correlated on those few
  // bits, so the next low-bit parity either keeps every model or kills every
  // model instead of roughly halving the count.
  auto makeBit = [&](const cvc5::Term& var, std::size_t bitIndex) {
    std::vector<uint32_t> params = {static_cast<uint32_t>(bitIndex),
                                    static_cast<uint32_t>(bitIndex)};
    cvc5::Op extract = tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, params);
    return tm.mkTerm(extract, {var});
  };

  std::vector<cvc5::Term> bits;
  bits.reserve(d_totalBits);
  std::bernoulli_distribution coin(0.5);
  for (const cvc5::Term& var : d_vars)
  {
    if (!var.getSort().isBitVector())
    {
      continue;
    }
    std::size_t width = var.getSort().getBitVectorSize();
    for (std::size_t b = 0; b < width; ++b)
    {
      if (coin(rng))
      {
        bits.push_back(makeBit(var, b));
      }
    }
  }

  // The coin flips can leave the parity empty; fall back to a single random bit
  // so the hash still constrains the space.
  if (bits.empty())
  {
    std::uniform_int_distribution<std::size_t> varDist(0, d_vars.size() - 1);
    while (bits.empty())
    {
      const cvc5::Term& var = d_vars[varDist(rng)];
      if (!var.getSort().isBitVector())
      {
        continue;
      }
      std::size_t width = var.getSort().getBitVectorSize();
      if (width == 0)
      {
        continue;
      }
      std::uniform_int_distribution<std::size_t> bitDist(0, width - 1);
      bits.push_back(makeBit(var, bitDist(rng)));
    }
  }

  cvc5::Term parity = bits.size() == 1 ? bits.front() : ttc::mkBvXor(tm, bits);

  unsigned rhsBit = static_cast<unsigned>(
      std::uniform_int_distribution<int>(0, 1)(rng));
  cvc5::Term rhs = tm.mkBitVector(1, rhsBit);
  std::vector<cvc5::Term> literals;
  literals.reserve(bits.size());
  cvc5::Term one = tm.mkBitVector(1, 1);
  for (const cvc5::Term& bit : bits)
  {
    literals.push_back(tm.mkTerm(
        cvc5::Kind::EQUAL, {bit, one}));
  }

  XorClause clause;
  clause.terms = std::move(literals);
  clause.rhs = rhsBit != 0;
  clause.fallback = tm.mkTerm(cvc5::Kind::EQUAL, {parity, rhs});
  return clause;
}

