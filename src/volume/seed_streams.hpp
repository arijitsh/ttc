#pragma once

// Independent random streams derived from --seed.  Each polytope's sampling
// walks and volume-oracle runs get their own seed: one base seed for every
// polytope would replay the same direction sequence in each, correlating the
// samples the union estimator combines (and, for the oracle, giving congruent
// polytopes identical errors that never average out).  Seeds are mixed with
// the splitmix64 finaliser so consecutive streams are not near-adjacent
// mt19937 states.

#include <cstddef>
#include <cstdint>

namespace ttc
{

inline unsigned splitmix64Finalize(std::uint64_t z)
{
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return static_cast<unsigned>(z ^ (z >> 31));
}

// The seeds of the sampling side: one per GenerateSamples call (randomized
// rounding) and one per batch of walk proposals, in call order.
class SeedStream
{
 public:
  explicit SeedStream(unsigned base) : d_base(base) {}
  unsigned next() { return splitmix64Finalize(d_base + 0x9E3779B97F4A7C15ULL * (++d_index)); }

 private:
  unsigned d_base;
  unsigned d_index = 0;
};

// Seed of the volume oracle's run `attempt` (0 first) on the polytope at
// 1-based position `polytope` among the kept ones.
inline unsigned oracleSeed(unsigned base, std::size_t polytope, unsigned attempt)
{
  return splitmix64Finalize((std::uint64_t(base) << 32) ^ 0xD1B54A32D192ED03ULL
                            ^ (0x9E3779B97F4A7C15ULL * (polytope + 1))
                            ^ (0xBF58476D1CE4E5B9ULL * attempt));
}

}  // namespace ttc
