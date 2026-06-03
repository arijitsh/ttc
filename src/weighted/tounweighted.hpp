#pragma once

#include <cstdint>
#include <vector>

namespace ttc::weighted {

struct LiteralWeight
{
  int var = 0;
  long double positive = 0.5L;
  long double negative = 0.5L;
};

struct WeightedCnf
{
  int varCount = 0;
  std::vector<std::vector<int>> clauses;
  std::vector<int> samplingVars;
  std::vector<LiteralWeight> weights;
};

struct UnweightedCnf
{
  int originalVarCount = 0;
  int originalClauseCount = 0;
  int varCount = 0;
  int divisorPower = 0;
  long double multiplier = 1.0L;
  std::vector<std::vector<int>> clauses;
  std::vector<int> samplingVars;
};

class ToUnweightedConverter
{
 public:
  explicit ToUnweightedConverter(int precision = 7);

  UnweightedCnf convert(const WeightedCnf& input) const;

 private:
  int d_precision;

  std::pair<std::uint64_t, int> quantizeWeight(long double weight) const;
};

}  // namespace ttc::weighted
