#include "pact/bv_hash.hpp"

#include "features.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace
{
// Smallest prime greater than 2^i, indexed by slice width i. Copied verbatim
// from cvc5's SmtApproxMc::populatePrimes so the Prime hash modulus matches.
const std::uint64_t kPrimes[] = {
    2,          2,          5,          11,        17,        37,
    67,         131,        257,        521,       1031,      2053,
    4099,       8209,       16411,      32771,     65537,     131101,
    262147,     524309,     1048583,    2097169,   4194319,   8388617,
    16777259,   33554467,   67108879,   134217757, 268435459, 536870923,
    1073741827, 2147483659, 4294967311, 8589934609};
constexpr std::size_t kPrimesCount = sizeof(kPrimes) / sizeof(kPrimes[0]);

// (1 << n) - 1 without invoking undefined behaviour for n >= 64.
std::uint64_t lowMask(std::uint32_t n)
{
  if (n >= 64)
  {
    return ~std::uint64_t(0);
  }
  return (std::uint64_t(1) << n) - 1;
}
}  // namespace

BvHash::BvHash(cvc5::Solver& solver)
    : d_solver(solver),
      d_totalBits(0),
      d_mode(BvHashMode::Xor),
      d_maxBitwidth(0),
      d_sliceSize(1)
{
  if (const char* mode = std::getenv("TTC_HASH_MODE"))
  {
    std::string m(mode);
    if (m == "prime")
    {
      d_mode = BvHashMode::Prime;
    }
    else if (m == "lemire")
    {
      d_mode = BvHashMode::Lemire;
    }
    else if (m == "prime-gj")
    {
      d_mode = BvHashMode::PrimeGj;
    }
  }
}

void BvHash::setVariables(std::vector<cvc5::Term> vars)
{
  d_vars = std::move(vars);
  d_totalBits = 0;
  d_maxBitwidth = 0;
  for (const cvc5::Term& v : d_vars)
  {
    if (v.getSort().isBitVector())
    {
      std::uint32_t width = v.getSort().getBitVectorSize();
      d_totalBits += width;
      d_maxBitwidth = std::max(d_maxBitwidth, width);
    }
  }

  // Word size used by the Prime/Lemire hashes, mirroring SmtApproxMc's
  // slice_size selection: half the widest variable, never wider than that
  // variable, and clamped so the prime table and (for Lemire) the 2*slice
  // accumulator stay in range. An explicit TTC_SLICE_SIZE overrides it.
  std::uint32_t slice = 0;
  if (const char* env = std::getenv("TTC_SLICE_SIZE"))
  {
    slice = static_cast<std::uint32_t>(std::atoi(env));
  }
  if (slice == 0)
  {
    slice = d_maxBitwidth / 2;
  }
  if (d_maxBitwidth > 0 && slice > d_maxBitwidth)
  {
    slice = d_maxBitwidth;
  }
  if (slice > 16)
  {
    // Cap so Lemire's 2*slice accumulator and the prime index stay safe in 64
    // bits. (cvc5 caps at 16 once the slice would exceed 32.)
    slice = 16;
  }
  if (slice == 0)
  {
    slice = 1;
  }
  if (slice >= kPrimesCount)
  {
    slice = static_cast<std::uint32_t>(kPrimesCount - 1);
  }
  d_sliceSize = slice;
}

double BvHash::perHashMultiplier() const
{
  switch (d_mode)
  {
    case BvHashMode::Prime:
    case BvHashMode::PrimeGj:
      return static_cast<double>(kPrimes[d_sliceSize]);
    case BvHashMode::Lemire:
      return std::pow(2.0, static_cast<double>(d_sliceSize));
    case BvHashMode::Xor:
    default:
      return 2.0;
  }
}

std::optional<XorClause> BvHash::randomClause(std::mt19937& rng) const
{
  if (d_totalBits == 0)
  {
    return std::nullopt;
  }

  if (d_mode == BvHashMode::Xor)
  {
    return randomXorClause(rng);
  }

  // PrimeGj: emit a mod-p row (literals + weights + modulus) for cvc5's Z_p
  // Gauss-Jordan propagator instead of a bit-vector arithmetic term.
  if (d_mode == BvHashMode::PrimeGj)
  {
    return generatePrimeRow(rng);
  }

  // Prime / Lemire are word-level arithmetic hashes: they have no native-XOR
  // representation, so they are carried as an XorClause with no parity literals
  // and an arithmetic fallback equality. assertToSolver sees the empty `terms`
  // and always asserts the fallback, regardless of the native-XOR setting.
  cvc5::Term constraint = (d_mode == BvHashMode::Prime)
                              ? generatePrimeHash(rng)
                              : generateLemireHash(rng);
  XorClause clause;
  clause.rhs = false;
  clause.fallback = std::move(constraint);
  return clause;
}

std::optional<XorClause> BvHash::randomXorClause(std::mt19937& rng) const
{
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
  double density = std::getenv("TTC_XOR_DENSITY")
                       ? std::atof(std::getenv("TTC_XOR_DENSITY"))
                       : 0.5;
  std::bernoulli_distribution coin(density);
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

std::uint32_t BvHash::getMinBitwidth() const
{
  // Minimum width of the modular accumulator: 2*slice + 1 plus enough bits to
  // hold the sum of all the slice products without overflow before the urem.
  std::uint32_t minBw = 2 * d_sliceSize + 1;
  std::uint32_t numSlicedVar = 0;
  for (const cvc5::Term& x : d_vars)
  {
    if (!x.getSort().isBitVector())
    {
      continue;
    }
    std::uint32_t width = x.getSort().getBitVectorSize();
    numSlicedVar += (width + d_sliceSize - 1) / d_sliceSize;
  }
  if (numSlicedVar == 0)
  {
    numSlicedVar = 1;
  }
  std::uint32_t extensionForSum = static_cast<std::uint32_t>(
      std::ceil(std::log2(static_cast<double>(numSlicedVar))));
  return minBw + extensionForSum;
}

cvc5::Term BvHash::generatePrimeHash(std::mt19937& rng) const
{
  auto& tm = ttc::getTermBuilder(d_solver);
  const std::uint32_t bitwidth = d_sliceSize;
  const std::uint32_t newBvWidth = getMinBitwidth();
  const std::uint64_t prime = kPrimes[bitwidth];

  cvc5::Term p = tm.mkBitVector(newBvWidth, prime);
  std::uniform_int_distribution<std::uint64_t> coefDist(0, prime - 1);
  const std::uint64_t cVal = coefDist(rng);

  cvc5::Term axpb = tm.mkBitVector(newBvWidth, std::uint64_t(0));
  cvc5::Term c = tm.mkBitVector(newBvWidth, cVal);

  for (const cvc5::Term& x : d_vars)
  {
    if (!x.getSort().isBitVector())
    {
      continue;
    }
    const std::uint32_t thisBvWidth = x.getSort().getBitVectorSize();
    std::uint32_t numSlices = thisBvWidth / bitwidth;
    if (bitwidth > thisBvWidth || numSlices == 0)
    {
      numSlices = 1;
    }
    for (std::uint32_t slice = 0; slice < numSlices; ++slice)
    {
      std::uint32_t sliceStart = slice * bitwidth;
      std::uint32_t sliceEnd = (slice + 1) * bitwidth - 1;
      std::uint32_t extendBy = newBvWidth - bitwidth;
      if (sliceEnd >= thisBvWidth)
      {
        extendBy = sliceEnd - thisBvWidth + 1 + newBvWidth - bitwidth;
        sliceEnd = thisBvWidth - 1;
      }
      const std::uint64_t aVal = coefDist(rng);

      cvc5::Op extractOp =
          tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, {sliceEnd, sliceStart});
      cvc5::Term xSliced = tm.mkTerm(extractOp, {x});
      cvc5::Op zeroExtOp =
          tm.mkOp(cvc5::Kind::BITVECTOR_ZERO_EXTEND, {extendBy});
      xSliced = tm.mkTerm(zeroExtOp, {xSliced});
      cvc5::Term a = tm.mkBitVector(newBvWidth, aVal);
      cvc5::Term ax = tm.mkTerm(cvc5::Kind::BITVECTOR_MULT, {a, xSliced});
      axpb = tm.mkTerm(cvc5::Kind::BITVECTOR_ADD, {ax, axpb});
    }
  }

  axpb = tm.mkTerm(cvc5::Kind::BITVECTOR_UREM, {axpb, p});
  return tm.mkTerm(cvc5::Kind::EQUAL, {axpb, c});
}

XorClause BvHash::generatePrimeRow(std::mt19937& rng) const
{
  // Same hash as generatePrimeHash:  sum_i a_i * x_i  ==  c  (mod p).  We emit
  // it two ways from one set of random draws: (1) as a mod-p row over the
  // individual projection bits -- bit (slice coeff a, position pos = start+k)
  // contributes (a * 2^k mod p) -- handed to the Z_p Gauss-Jordan propagator;
  // and (2) as the equivalent bit-vector arithmetic term in `fallback`, used
  // only by the cheap model-reuse check (never asserted, so never blasted).
  auto& tm = ttc::getTermBuilder(d_solver);
  const std::uint32_t bitwidth = d_sliceSize;
  const std::uint32_t newBvWidth = getMinBitwidth();
  const std::uint64_t prime = kPrimes[bitwidth];

  cvc5::Term p = tm.mkBitVector(newBvWidth, prime);
  std::uniform_int_distribution<std::uint64_t> coefDist(0, prime - 1);
  const std::uint64_t cVal = coefDist(rng);

  cvc5::Term axpb = tm.mkBitVector(newBvWidth, std::uint64_t(0));
  cvc5::Term c = tm.mkBitVector(newBvWidth, cVal);
  cvc5::Term one = tm.mkBitVector(1, std::uint64_t(1));

  std::vector<cvc5::Term> literals;
  std::vector<std::uint64_t> weights;

  for (const cvc5::Term& x : d_vars)
  {
    if (!x.getSort().isBitVector())
    {
      continue;
    }
    const std::uint32_t thisBvWidth = x.getSort().getBitVectorSize();
    std::uint32_t numSlices = thisBvWidth / bitwidth;
    if (bitwidth > thisBvWidth || numSlices == 0)
    {
      numSlices = 1;
    }
    for (std::uint32_t slice = 0; slice < numSlices; ++slice)
    {
      std::uint32_t sliceStart = slice * bitwidth;
      std::uint32_t sliceEnd = (slice + 1) * bitwidth - 1;
      std::uint32_t extendBy = newBvWidth - bitwidth;
      if (sliceEnd >= thisBvWidth)
      {
        extendBy = sliceEnd - thisBvWidth + 1 + newBvWidth - bitwidth;
        sliceEnd = thisBvWidth - 1;
      }
      const std::uint64_t aVal = coefDist(rng);

      // (1) Fallback bit-vector arithmetic term (identical to generatePrimeHash).
      cvc5::Op extractOp =
          tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, {sliceEnd, sliceStart});
      cvc5::Term xSliced = tm.mkTerm(extractOp, {x});
      cvc5::Op zeroExtOp =
          tm.mkOp(cvc5::Kind::BITVECTOR_ZERO_EXTEND, {extendBy});
      xSliced = tm.mkTerm(zeroExtOp, {xSliced});
      cvc5::Term a = tm.mkBitVector(newBvWidth, aVal);
      cvc5::Term ax = tm.mkTerm(cvc5::Kind::BITVECTOR_MULT, {a, xSliced});
      axpb = tm.mkTerm(cvc5::Kind::BITVECTOR_ADD, {ax, axpb});

      // (2) Mod-p row: one literal per bit of this slice, weight a*2^k mod p.
      std::uint64_t pow2 = 1 % prime;  // 2^k mod p, starting at k = 0
      for (std::uint32_t pos = sliceStart; pos <= sliceEnd; ++pos)
      {
        std::uint64_t w = (aVal % prime) * pow2 % prime;
        if (w != 0)
        {
          cvc5::Op bitOp = tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, {pos, pos});
          cvc5::Term bit = tm.mkTerm(bitOp, {x});
          literals.push_back(tm.mkTerm(cvc5::Kind::EQUAL, {bit, one}));
          weights.push_back(w);
        }
        pow2 = pow2 * 2 % prime;
      }
    }
  }

  axpb = tm.mkTerm(cvc5::Kind::BITVECTOR_UREM, {axpb, p});

  XorClause clause;
  clause.terms = std::move(literals);
  clause.weights = std::move(weights);
  clause.modulus = prime;
  clause.rhsValue = cVal;
  clause.rhs = false;
  clause.fallback = tm.mkTerm(cvc5::Kind::EQUAL, {axpb, c});
  return clause;
}

cvc5::Term BvHash::generateLemireHash(std::mt19937& rng) const
{
  auto& tm = ttc::getTermBuilder(d_solver);
  const std::uint32_t bitwidth = d_sliceSize;
  const std::uint32_t newBvWidth = bitwidth * 2;

  std::uniform_int_distribution<std::uint64_t> cDist(0, lowMask(bitwidth));
  const std::uint64_t cVal = cDist(rng);
  std::uniform_int_distribution<std::uint64_t> aDist(0, lowMask(2 * bitwidth));

  cvc5::Term axpb = tm.mkBitVector(newBvWidth, std::uint64_t(0));
  cvc5::Term c = tm.mkBitVector(bitwidth, cVal);

  for (const cvc5::Term& x : d_vars)
  {
    if (!x.getSort().isBitVector())
    {
      continue;
    }
    const std::uint32_t thisBvWidth = x.getSort().getBitVectorSize();
    std::uint32_t numSlices = thisBvWidth / bitwidth;
    if (bitwidth > thisBvWidth || numSlices == 0)
    {
      numSlices = 1;
    }
    for (std::uint32_t slice = 0; slice < numSlices; ++slice)
    {
      std::uint32_t sliceStart = slice * bitwidth;
      std::uint32_t sliceEnd = (slice + 1) * bitwidth - 1;
      std::uint32_t extendBy = bitwidth;
      if (sliceEnd >= thisBvWidth)
      {
        extendBy = sliceEnd - thisBvWidth + bitwidth - 1;
        sliceEnd = thisBvWidth - 1;
      }
      const std::uint64_t aVal = aDist(rng);

      cvc5::Op extractOp =
          tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, {sliceEnd, sliceStart});
      cvc5::Term xSliced = tm.mkTerm(extractOp, {x});
      cvc5::Op zeroExtOp =
          tm.mkOp(cvc5::Kind::BITVECTOR_ZERO_EXTEND, {extendBy});
      xSliced = tm.mkTerm(zeroExtOp, {xSliced});
      cvc5::Term a = tm.mkBitVector(newBvWidth, aVal);
      cvc5::Term ax = tm.mkTerm(cvc5::Kind::BITVECTOR_MULT, {a, xSliced});
      axpb = tm.mkTerm(cvc5::Kind::BITVECTOR_ADD, {ax, axpb});
    }
  }

  // Keep the high `bitwidth` bits of the accumulator (multiply-shift): the
  // product fills 2*bitwidth bits and the hash value is its top half.
  cvc5::Op highOp =
      tm.mkOp(cvc5::Kind::BITVECTOR_EXTRACT, {bitwidth * 2 - 1, bitwidth});
  cvc5::Term axpbHigh = tm.mkTerm(highOp, {axpb});
  return tm.mkTerm(cvc5::Kind::EQUAL, {axpbHigh, c});
}
