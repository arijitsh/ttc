#include "arjun.hpp"

#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "features.hpp"

namespace arjun {

void minimizeProjectionSet(cvc5::Solver& solver,
                           const cvc5::Term& formula,
                           std::vector<cvc5::Term>& projVars)
{
    if (projVars.empty())
    {
        return;
    }

    auto& tm = ttc::getTermBuilder(solver);

    // Preserve solver state
    solver.push();

    // Collect all free variables in the formula
    std::vector<cvc5::Term> allVars;
    std::vector<cvc5::Term> allCopies;
    std::unordered_map<std::string, cvc5::Term> copyMap;
    std::unordered_set<std::string> seen;
    std::function<void(const cvc5::Term&)> collect = [&](const cvc5::Term& t) {
        if (t.getNumChildren() == 0 && t.hasSymbol())
        {
            std::string name = t.toString();
            if (seen.insert(name).second)
            {
                allVars.push_back(t);
                cvc5::Term copy = tm.mkConst(t.getSort(), "z_" + name);
                allCopies.push_back(copy);
                copyMap.emplace(name, copy);
            }
        }
        for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
        {
            collect(t[i]);
        }
    };
    collect(formula);

    // Build phi with all variables replaced by their fresh copies
    cvc5::Term phiZ = formula.substitute(allVars, allCopies);

    // Assert base formula
    solver.assertFormula(tm.mkTerm(cvc5::Kind::AND, {formula, phiZ}));

    std::vector<cvc5::Term> support;
    support.reserve(projVars.size());

    // Padoa check for each projection variable
    for (const cvc5::Term& v : projVars)
    {
        solver.push();
        // Equate variables already in the support
        for (const cvc5::Term& s : support)
        {
            cvc5::Term zs = copyMap[s.toString()];
            solver.assertFormula(tm.mkTerm(cvc5::Kind::EQUAL, {s, zs}));
        }
        // Assert inequality for current variable
        cvc5::Term zv = copyMap[v.toString()];
        cvc5::Term neq = tm.mkTerm(cvc5::Kind::NOT,
                                   {tm.mkTerm(cvc5::Kind::EQUAL, {v, zv})});
        solver.assertFormula(neq);

        cvc5::Result res = solver.checkSat();
        solver.pop();

        if (!res.isUnsat())
        {
            support.push_back(v);
        }
    }

    projVars.swap(support);
    solver.pop();
}

} // namespace arjun

