#include "cache.hpp"

#include "features.hpp"

bool termsEquivalent(cvc5::Solver& solver, const cvc5::Term& a, const cvc5::Term& b)
{
    if (a == b)
    {
        return true;
    }
    auto& tm = ttc::getTermBuilder(solver);
    cvc5::Term xorTerm = tm.mkTerm(cvc5::Kind::XOR, {a, b});
    cvc5::Solver check = ttc::makeSolverWithBuilder(solver);
    check.assertFormula(xorTerm);
    cvc5::Result r = check.checkSat();
    return r.isUnsat();
}
