#pragma once

#include <cvc5/cvc5.h>

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <iostream>
namespace ttc
{
#if defined(TTC_ENABLE_DDNNF)
constexpr bool kDdnnfEnabled = true;
#else
constexpr bool kDdnnfEnabled = false;
#endif

inline bool ddnnfEnabled()
{
  return kDdnnfEnabled;
}

template <typename SolverT, typename = void>
struct TermBuilderHelper
{
  using storage_type = std::monostate;

  static SolverT createSolver(storage_type&)
  {
    return SolverT();
  }

  static SolverT duplicateSolver(SolverT&)
  {
    return SolverT();
  }

  static SolverT& builder(SolverT& solver)
  {
    return solver;
  }
};

template <typename SolverT>
struct TermBuilderHelper<SolverT,
                         std::void_t<decltype(std::declval<SolverT&>().getTermManager())>>
{
  using builder_type = decltype(std::declval<SolverT&>().getTermManager());
  using storage_type = std::remove_reference_t<builder_type>;

  static SolverT createSolver(storage_type& storage)
  {
    return SolverT(storage);
  }

  static SolverT duplicateSolver(SolverT& solver)
  {
    return SolverT(solver.getTermManager());
  }

  static builder_type builder(SolverT& solver)
  {
    return solver.getTermManager();
  }
};

template <typename SolverT>
auto getTermBuilder(SolverT& solver)
    -> decltype(TermBuilderHelper<SolverT>::builder(solver))
{
  return TermBuilderHelper<SolverT>::builder(solver);
}

template <typename SolverT>
SolverT makeSolverWithBuilder(SolverT& solver)
{
  return TermBuilderHelper<SolverT>::duplicateSolver(solver);
}

template <typename SolverT>
SolverT createSolverWithStorage(typename TermBuilderHelper<SolverT>::storage_type& storage)
{
  return TermBuilderHelper<SolverT>::createSolver(storage);
}

struct BooleanAbstractionData
{
  std::vector<uint32_t> clauses;
  std::vector<std::pair<uint32_t, cvc5::Term>> mapping;
};

struct BooleanAbstractionAigData
{
  struct AndGate
  {
    uint32_t lhs = 0;
    uint32_t rhs0 = 0;
    uint32_t rhs1 = 0;
  };

  uint32_t maxVariable = 0;
  std::vector<uint32_t> inputs;
  std::vector<uint32_t> outputs;
  std::vector<AndGate> andGates;
  std::vector<std::pair<uint32_t, cvc5::Term>> mapping;
};

template <typename SolverT, typename = void>
struct BooleanAbstractionHelper
{
  static BooleanAbstractionData get(SolverT&, const cvc5::Term&)
  {
    throw std::runtime_error("cvc5 Boolean abstraction is not supported in this build");
  }
};

template <typename SolverT>
struct BooleanAbstractionHelper<
    SolverT,
    std::void_t<decltype(std::declval<SolverT&>().getBooleanAbstraction(std::declval<cvc5::Term>()))>>
{
  static BooleanAbstractionData get(SolverT& solver, const cvc5::Term& term)
  {
    auto abstraction = solver.getBooleanAbstraction(term);
    BooleanAbstractionData data;
    data.clauses = abstraction.first;
    for (const auto& entry : abstraction.second)
    {
      data.mapping.emplace_back(static_cast<uint32_t>(entry.first), entry.second);
    }
    return data;
  }
};

template <typename SolverT>
BooleanAbstractionData getBooleanAbstraction(SolverT& solver, const cvc5::Term& term)
{
  return BooleanAbstractionHelper<SolverT>::get(solver, term);
}

template <typename SolverT, typename = void>
struct BooleanAbstractionAigHelper
{
  static BooleanAbstractionAigData get(SolverT&, const cvc5::Term&)
  {
    throw std::runtime_error(
        "cvc5 Boolean AIG abstraction is not supported in this build");
  }
};

template <typename SolverT>
struct BooleanAbstractionAigHelper<
    SolverT,
    std::void_t<decltype(
        std::declval<SolverT&>().getBooleanAbstractionAig(
            std::declval<cvc5::Term>()))>>
{
  static BooleanAbstractionAigData get(SolverT& solver, const cvc5::Term& term)
  {
    auto abstraction = solver.getBooleanAbstractionAig(term);
    BooleanAbstractionAigData data;
    data.maxVariable = abstraction.maxVariable;
    data.inputs = abstraction.inputs;
    data.outputs = abstraction.outputs;

    data.andGates.reserve(abstraction.andGates.size());
    for (const auto& gate : abstraction.andGates)
    {
      data.andGates.push_back({gate.lhs, gate.rhs0, gate.rhs1});
    }
    data.mapping.reserve(abstraction.variableMap.size());
    for (const auto& entry : abstraction.variableMap)
    {
      data.mapping.emplace_back(entry.first, entry.second);
    }
    std::sort(data.mapping.begin(), data.mapping.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::cout << "c AIG variables: " << data.maxVariable
                << " and gates: " << data.andGates.size()
                << std::endl;
    return data;
  }
};

template <typename SolverT>
BooleanAbstractionAigData getBooleanAbstractionAig(SolverT& solver,
                                                   const cvc5::Term& term)
{
  return BooleanAbstractionAigHelper<SolverT>::get(solver, term);
}

template <typename Builder>
auto mkBvXorImpl(Builder& builder, const std::vector<cvc5::Term>& terms, int)
    -> decltype(builder.mkBvXor(terms))
{
  return builder.mkBvXor(terms);
}

template <typename Builder>
cvc5::Term mkBvXorImpl(Builder& builder, const std::vector<cvc5::Term>& terms, long)
{
  return builder.mkTerm(cvc5::Kind::BITVECTOR_XOR, terms);
}

template <typename Builder>
cvc5::Term mkBvXor(Builder& builder, const std::vector<cvc5::Term>& terms)
{
  return mkBvXorImpl(builder, terms, 0);
}

}  // namespace ttc
