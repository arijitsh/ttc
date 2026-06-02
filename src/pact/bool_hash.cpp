#include "pact/bool_hash.hpp"

#include "features.hpp"

#include <random>
#include <utility>
#include <vector>

BoolHash::BoolHash(cvc5::Solver& solver) : d_solver(solver) {}

void BoolHash::setVariables(std::vector<cvc5::Term> boolVars,
                            std::vector<cvc5::Term> bvVars)
{
  d_boolVars = std::move(boolVars);
  d_bvVars = std::move(bvVars);
}

std::optional<XorClause> BoolHash::randomClause(std::mt19937& rng) const
{
  if (d_boolVars.empty())
  {
    return std::nullopt;
  }

  if (d_bvVars.size() != d_boolVars.size())
  {
    return std::nullopt;
  }

  auto& tm = ttc::getTermBuilder(d_solver);

  std::vector<std::size_t> selectedIndices;
  selectedIndices.reserve(d_boolVars.size());
  std::bernoulli_distribution include(0.5);
  for (std::size_t i = 0; i < d_boolVars.size(); ++i)
  {
    if (include(rng))
    {
      selectedIndices.push_back(i);
    }
  }
  if (selectedIndices.empty())
  {
    std::uniform_int_distribution<std::size_t> pick(0, d_boolVars.size() - 1);
    selectedIndices.push_back(pick(rng));
  }

  std::vector<cvc5::Term> boolTerms;
  boolTerms.reserve(selectedIndices.size());
  std::vector<cvc5::Term> bvTerms;
  bvTerms.reserve(selectedIndices.size());
  for (std::size_t index : selectedIndices)
  {
    boolTerms.push_back(d_boolVars[index]);
    if (index < d_bvVars.size())
    {
      bvTerms.push_back(d_bvVars[index]);
    }
  }

  cvc5::Term parity = bvTerms.size() == 1 ? bvTerms.front()
                                           : ttc::mkBvXor(tm, bvTerms);
  bool rhsBit = std::uniform_int_distribution<int>(0, 1)(rng) != 0;
  cvc5::Term rhs = tm.mkBitVector(1, rhsBit ? 1u : 0u);
  XorClause clause;
  clause.terms = std::move(boolTerms);
  clause.rhs = rhsBit;
  clause.fallback = tm.mkTerm(cvc5::Kind::EQUAL, {parity, rhs});
  return clause;
}

