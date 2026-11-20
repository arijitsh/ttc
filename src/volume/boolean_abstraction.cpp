#include "volume/boolean_abstraction.hpp"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <utility>

namespace ttc
{
namespace
{

// AIGER literals: 2 v is variable v, 2 v + 1 its negation; 0 is false, 1 true.
constexpr std::uint32_t kFalse = 0U;
constexpr std::uint32_t kTrue = 1U;

std::uint32_t negate(std::uint32_t literal) { return literal ^ 1U; }

struct LiteralPairHash
{
  std::size_t operator()(const std::pair<std::uint32_t, std::uint32_t>& p) const noexcept
  {
    return (static_cast<std::size_t>(p.first) << 32) ^ static_cast<std::size_t>(p.second);
  }
};

class AigBuilder
{
 public:
  BooleanAbstractionAigData build(const cvc5::Term& formula)
  {
    d_data.maxVariable = 1;
    d_data.outputs.push_back(visit(formula));
    std::sort(d_data.mapping.begin(), d_data.mapping.end(),
              [](const auto& x, const auto& y) { return x.first < y.first; });
    std::sort(d_data.inputs.begin(), d_data.inputs.end());
    return std::move(d_data);
  }

 private:
  std::uint32_t newLiteral()
  {
    const std::uint32_t literal = d_nextLiteral;
    d_nextLiteral += 2;
    d_data.maxVariable = std::max(d_data.maxVariable, literal >> 1U);
    return literal;
  }

  // An input of the AIG, recorded in M.
  std::uint32_t atom(const cvc5::Term& term)
  {
    if (auto it = d_atoms.find(term); it != d_atoms.end()) return it->second;
    const std::uint32_t literal = newLiteral();
    d_atoms.emplace(term, literal);
    d_data.inputs.push_back(literal);
    d_data.mapping.emplace_back(literal >> 1U, term);
    return literal;
  }

  std::uint32_t mkAnd(std::uint32_t x, std::uint32_t y)
  {
    if (x == kFalse || y == kFalse || x == negate(y)) return kFalse;
    if (x == kTrue || x == y) return y;
    if (y == kTrue) return x;
    if (x > y) std::swap(x, y);
    const auto key = std::make_pair(x, y);
    if (auto it = d_gates.find(key); it != d_gates.end()) return it->second;
    const std::uint32_t literal = newLiteral();
    d_data.andGates.push_back({literal, x, y});
    d_gates.emplace(key, literal);
    return literal;
  }

  std::uint32_t mkOr(std::uint32_t x, std::uint32_t y)
  {
    if (x == kTrue || y == kTrue || x == negate(y)) return kTrue;
    if (x == kFalse || x == y) return y;
    if (y == kFalse) return x;
    return negate(mkAnd(negate(x), negate(y)));
  }

  // (a1 & b1) | (a2 & b2), creating the first gate first.  The order fixes
  // the literal numbering, hence AllSAT's cube order and the estimate on a
  // given seed, so every visit and gate below is sequenced explicitly (C++
  // leaves the evaluation order of function arguments unspecified).
  std::uint32_t mkOrOfAnds(std::uint32_t a1, std::uint32_t b1, std::uint32_t a2, std::uint32_t b2)
  {
    const std::uint32_t first = mkAnd(a1, b1);
    const std::uint32_t second = mkAnd(a2, b2);
    return mkOr(first, second);
  }

  std::uint32_t visit(const cvc5::Term& term)
  {
    if (auto it = d_visited.find(term); it != d_visited.end()) return it->second;

    std::uint32_t result = kFalse;
    switch (term.getKind())
    {
      case cvc5::Kind::CONST_BOOLEAN:
        return term.getBooleanValue() ? kTrue : kFalse;
      case cvc5::Kind::NOT: result = negate(visit(term[0])); break;
      case cvc5::Kind::AND:
        result = term.getNumChildren() == 0 ? kTrue : visit(term[0]);
        for (std::size_t i = 1; i < term.getNumChildren(); ++i) result = mkAnd(result, visit(term[i]));
        break;
      case cvc5::Kind::OR:
        result = term.getNumChildren() == 0 ? kFalse : visit(term[0]);
        for (std::size_t i = 1; i < term.getNumChildren(); ++i) result = mkOr(result, visit(term[i]));
        break;
      case cvc5::Kind::IMPLIES:
      {
        const std::uint32_t lhs = visit(term[0]);
        const std::uint32_t rhs = visit(term[1]);
        result = mkOr(negate(lhs), rhs);
        break;
      }
      case cvc5::Kind::ITE:
      {
        const std::uint32_t cond = visit(term[0]);
        const std::uint32_t thenLit = visit(term[1]);
        const std::uint32_t elseLit = visit(term[2]);
        result = mkOrOfAnds(cond, thenLit, negate(cond), elseLit);
        break;
      }
      case cvc5::Kind::XOR:
      {
        const std::uint32_t a = visit(term[0]);
        const std::uint32_t b = visit(term[1]);
        result = mkOrOfAnds(negate(a), b, a, negate(b));
        break;
      }
      case cvc5::Kind::EQUAL:
        // Boolean equality is (a & b) | (!a & !b); an equality of numbers is an atom.
        if (term[0].getSort().isBoolean())
        {
          const std::uint32_t a = visit(term[0]);
          const std::uint32_t b = visit(term[1]);
          result = mkOrOfAnds(a, b, negate(a), negate(b));
        }
        else
          result = atom(term);
        break;
      default: result = atom(term); break;
    }
    if (term.getSort().isBoolean()) d_visited.emplace(term, result);
    return result;
  }

  BooleanAbstractionAigData d_data;
  std::uint32_t d_nextLiteral = 2;
  std::unordered_map<cvc5::Term, std::uint32_t> d_visited;
  std::unordered_map<cvc5::Term, std::uint32_t> d_atoms;
  std::unordered_map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t, LiteralPairHash>
      d_gates;
};

}  // namespace

BooleanAbstractionAigData buildBooleanAigFromTerm(const cvc5::Term& formula)
{
  return AigBuilder().build(formula);
}

}  // namespace ttc
