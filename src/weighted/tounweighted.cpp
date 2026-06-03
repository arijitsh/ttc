#include "weighted/tounweighted.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace ttc::weighted {
namespace {

void appendUnique(std::vector<int>& values, int value)
{
  if (std::find(values.begin(), values.end(), value) == values.end())
  {
    values.push_back(value);
  }
}

void pushVar(int var, std::vector<std::vector<int>>& clauses)
{
  for (auto& clause : clauses)
  {
    clause.push_back(var);
  }
}

std::string binaryWithoutLeastSignificantBit(std::uint64_t value)
{
  if (value == 0)
  {
    return {};
  }
  std::string bits;
  while (value > 0)
  {
    bits.push_back((value & 1U) ? '1' : '0');
    value >>= 1U;
  }
  std::reverse(bits.begin(), bits.end());
  bits.pop_back();
  return bits;
}

std::vector<std::vector<int>> getCnf(int var,
                                     const std::string& binStr,
                                     bool sign,
                                     int origVars)
{
  std::vector<std::vector<int>> clauses;
  const int binLen = static_cast<int>(binStr.size());
  clauses.push_back({binLen + 1 + origVars});
  for (int i = 0; i < binLen; ++i)
  {
    int newVar = binLen - i + origVars;
    if (!sign)
    {
      newVar = -newVar;
    }
    if (binStr[binLen - i - 1] == '0')
    {
      clauses.push_back({newVar});
    }
    else
    {
      pushVar(newVar, clauses);
    }
  }
  pushVar(var, clauses);
  return clauses;
}

bool containsClause(const std::vector<std::vector<int>>& clauses,
                    const std::vector<int>& clause)
{
  return std::find(clauses.begin(), clauses.end(), clause) != clauses.end();
}

}  // namespace

ToUnweightedConverter::ToUnweightedConverter(int precision)
    : d_precision(precision)
{
  if (d_precision < 2)
  {
    throw std::runtime_error("weighted-to-unweighted precision must be at least 2");
  }
}

std::pair<std::uint64_t, int> ToUnweightedConverter::quantizeWeight(
    long double weight) const
{
  if (weight < 0.0L || weight > 1.0L)
  {
    throw std::runtime_error("literal weight is outside [0, 1]");
  }

  const long double scaled = weight * std::ldexp(1.0L, d_precision);
  std::uint64_t numerator =
      static_cast<std::uint64_t>(std::llround(scaled));
  int precision = d_precision;
  while (numerator > 0 && (numerator % 2U) == 0U && precision > 0)
  {
    numerator /= 2U;
    --precision;
  }
  if (numerator == 0)
  {
    throw std::runtime_error(
        "literal weight rounded to zero; increase conversion precision");
  }
  return {numerator, precision};
}

UnweightedCnf ToUnweightedConverter::convert(const WeightedCnf& input) const
{
  if (input.varCount < 0)
  {
    throw std::runtime_error("CNF variable count cannot be negative");
  }

  UnweightedCnf output;
  output.originalVarCount = input.varCount;
  output.originalClauseCount = static_cast<int>(input.clauses.size());
  output.varCount = input.varCount;
  output.clauses = input.clauses;

  for (int var : input.samplingVars)
  {
    if (var <= 0 || var > input.varCount)
    {
      throw std::runtime_error("sampling variable is outside the CNF variable range");
    }
    appendUnique(output.samplingVars, var);
  }

  std::unordered_set<int> samplingSet(output.samplingVars.begin(),
                                      output.samplingVars.end());
  std::unordered_map<int, LiteralWeight> weightsByVar;
  for (const LiteralWeight& weight : input.weights)
  {
    if (weight.var <= 0 || weight.var > input.varCount)
    {
      throw std::runtime_error("weighted variable is outside the CNF variable range");
    }
    if (samplingSet.find(weight.var) == samplingSet.end())
    {
      continue;
    }
    weightsByVar[weight.var] = weight;
  }

  for (const auto& [var, rawWeight] : weightsByVar)
  {
    long double positive = rawWeight.positive;
    long double negative = rawWeight.negative;
    if (positive < 0.0L || negative < 0.0L)
    {
      throw std::runtime_error("negative literal weights are not supported");
    }
    if (positive == 1.0L && negative == 1.0L)
    {
      continue;
    }
    if (positive == 0.0L || negative == 0.0L)
    {
      throw std::runtime_error(
          "zero literal weights require preprocessing and are not supported in --PB");
    }

    const long double total = positive + negative;
    if (total <= 0.0L)
    {
      throw std::runtime_error("literal weights for a variable sum to zero");
    }
    output.multiplier *= total;
    positive /= total;
    negative /= total;

    if (positive == 1.0L && negative == 1.0L)
    {
      continue;
    }

    auto [bitMult, bitPrec] = quantizeWeight(positive);
    if (bitPrec == 1 && bitMult == 1)
    {
      ++output.divisorPower;
      continue;
    }
    if (bitPrec == 0)
    {
      throw std::runtime_error(
          "weight cannot be encoded without zero-weight preprocessing");
    }

    const int baseVars = output.varCount;
    for (int i = 1; i <= bitPrec; ++i)
    {
      appendUnique(output.samplingVars, baseVars + i);
    }

    std::string binStr = binaryWithoutLeastSignificantBit(bitMult);
    if (static_cast<int>(binStr.size()) > bitPrec - 1)
    {
      throw std::runtime_error("quantized weight exceeds selected precision");
    }
    while (static_cast<int>(binStr.size()) < bitPrec - 1)
    {
      binStr.insert(binStr.begin(), '0');
    }

    std::string complementStr;
    complementStr.reserve(binStr.size());
    for (char bit : binStr)
    {
      complementStr.push_back(bit == '0' ? '1' : '0');
    }

    std::vector<std::vector<int>> positiveClauses =
        getCnf(-var, binStr, true, baseVars);
    output.clauses.insert(output.clauses.end(),
                          positiveClauses.begin(),
                          positiveClauses.end());

    std::vector<std::vector<int>> negativeClauses =
        getCnf(var, complementStr, false, baseVars);
    for (const auto& clause : negativeClauses)
    {
      if (!containsClause(positiveClauses, clause))
      {
        output.clauses.push_back(clause);
      }
    }

    output.varCount += bitPrec;
    output.divisorPower += bitPrec;
  }

  return output;
}

}  // namespace ttc::weighted
