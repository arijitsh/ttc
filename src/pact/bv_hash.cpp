#include "pact/bv_hash.hpp"

#include "features.hpp"

#include <algorithm>
#include <cstdint>
#include <random>
#include <unordered_set>
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

  std::size_t maxBits = std::min<std::size_t>(d_totalBits, static_cast<std::size_t>(4));
  std::uniform_int_distribution<std::size_t> countDist(1, maxBits);
  std::size_t targetCount = countDist(rng);

  std::vector<cvc5::Term> bits;
  bits.reserve(targetCount);
  std::unordered_set<std::uint64_t> seen;

  std::uniform_int_distribution<std::size_t> varDist(0, d_vars.size() - 1);

  auto encode = [](std::size_t varIndex, std::size_t bitIndex) {
    return (static_cast<std::uint64_t>(varIndex) << 32)
           | static_cast<std::uint64_t>(bitIndex);
  };

  std::size_t attempts = 0;
  while (bits.size() < targetCount && attempts < targetCount * 16)
  {
    ++attempts;
    std::size_t varIdx = varDist(rng);
    const cvc5::Term& var = d_vars[varIdx];
    std::size_t width = var.getSort().getBitVectorSize();
    if (width == 0)
    {
      continue;
    }
    std::uniform_int_distribution<std::size_t> bitDist(0, width - 1);
    std::size_t bitIndex = bitDist(rng);
    std::uint64_t key = encode(varIdx, bitIndex);
    if (!seen.insert(key).second)
    {
      continue;
    }
    std::vector<uint32_t> params = {static_cast<uint32_t>(bitIndex),
                                    static_cast<uint32_t>(bitIndex)};
    cvc5::Op extract = tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, params);
    bits.push_back(tm.mkTerm(extract, {var}));
  }

  if (bits.empty())
  {
    const cvc5::Term& var = d_vars.front();
    std::size_t width = var.getSort().getBitVectorSize();
    if (width == 0)
    {
      return std::nullopt;
    }
    std::vector<uint32_t> params = {static_cast<uint32_t>(width - 1),
                                    static_cast<uint32_t>(width - 1)};
    cvc5::Op extract = tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, params);
    bits.push_back(tm.mkTerm(extract, {var}));
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

