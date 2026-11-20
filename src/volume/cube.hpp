#pragma once

// A cube c_i of the DNF F_bd (Algorithm 1, line 2), read through the map M
// from abstraction variables to LRA atoms.  Its literals are the linear
// constraints whose conjunction is the polytope P(c_i) (linear_constraints.hpp).

#include <vector>

#include <cvc5/cvc5.h>

namespace ttc
{

// The atom `term` (a linear comparison) or, if !positive, its negation.
struct TermLiteral
{
  cvc5::Term term;
  bool positive = true;
};

struct Cube
{
  // The atoms the cube fixes, each once; atoms it leaves open are absent.
  std::vector<TermLiteral> literals;
};

}  // namespace ttc
