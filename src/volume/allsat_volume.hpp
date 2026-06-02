#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <cvc5/cvc5.h>

#include "features.hpp"

namespace ttc
{

struct BooleanLiteral
{
  std::uint32_t variable = 0;
  std::optional<bool> value;
};

struct TermLiteral
{
  cvc5::Term term;
  std::optional<bool> value;
};

struct Polytope
{
  std::vector<BooleanLiteral> booleanLiterals;
  std::vector<TermLiteral> termLiterals;
  cvc5::Term formula;
  bool islChecked = false;
  bool islNonEmpty = true;
};

struct AllSatVolumeResult
{
  std::vector<Polytope> polytopes;
};

AllSatVolumeResult enumeratePolytopes(const BooleanAbstractionAigData& data,
                                      cvc5::Solver& solver);

BooleanAbstractionAigData buildBooleanAigFromTerm(const cvc5::Term& term);

}  // namespace ttc

