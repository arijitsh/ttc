#include "arjun.hpp"

#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "features.hpp"

namespace arjun {

Reduction analyzeProjectionSet(cvc5::Solver& solver,
                               const cvc5::Term& formula,
                               const std::vector<cvc5::Term>& projVars)
{
    Reduction result;
    if (projVars.empty())
    {
        return result;
    }

    auto& tm = ttc::getTermBuilder(solver);

    // ---- Pass 0: unconstrained (free) variable detection -----------------
    // Gather the names of every variable that actually occurs in the formula.
    // A projection variable absent from this set is not constrained by any
    // assertion: both of its values extend any model, so it is removed here.
    // (This also keeps such variables out of the definability pass below,
    // whose copy map only covers variables that occur in the formula.)
    std::unordered_set<std::string> occurring;
    std::function<void(const cvc5::Term&)> gather = [&](const cvc5::Term& t) {
        if (t.getNumChildren() == 0 && t.hasSymbol())
        {
            occurring.insert(t.toString());
        }
        for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
        {
            gather(t[i]);
        }
    };
    gather(formula);

    std::vector<cvc5::Term> present;
    present.reserve(projVars.size());
    for (const cvc5::Term& v : projVars)
    {
        if (occurring.count(v.toString()) != 0)
        {
            present.push_back(v);
        }
        else
        {
            result.free.push_back(v);
        }
    }
    if (present.empty())
    {
        return result;
    }

    // ---- Pass 1: backbone (forced-constant) detection --------------------
    // A Boolean projection variable is backbone if the formula has no model in
    // which it takes the value opposite to the one seen in some model. We grab
    // a single model, then for each variable test only the opposite value: one
    // extra SAT call per variable, instead of two. Bit-vector projection vars
    // are skipped here (only definability applies to them).
    std::vector<cvc5::Term> candidates;
    candidates.reserve(present.size());
    {
        solver.push();
        solver.assertFormula(formula);
        cvc5::Result base = solver.checkSat();
        if (!base.isSat())
        {
            // Unsatisfiable (or unknown): no safe reduction; keep everything.
            solver.pop();
            result.free.clear();
            result.support = projVars;
            return result;
        }
        // Snapshot the model values of every Boolean projection variable now,
        // while this SAT result is current. getValue() is only legal right
        // after a SAT/UNKNOWN checkSat, and the per-variable push/pop below
        // ends with an (often UNSAT) result, so the values cannot be queried
        // lazily inside the loop.
        std::unordered_map<std::string, bool> modelValues;
        for (const cvc5::Term& v : present)
        {
            if (v.getSort().isBoolean())
            {
                modelValues.emplace(v.toString(),
                                    solver.getValue(v).getBooleanValue());
            }
        }
        for (const cvc5::Term& v : present)
        {
            if (!v.getSort().isBoolean())
            {
                candidates.push_back(v);
                continue;
            }
            bool modelValue = modelValues.at(v.toString());
            cvc5::Term opposite =
                modelValue ? tm.mkTerm(cvc5::Kind::NOT, {v}) : v;
            solver.push();
            solver.assertFormula(opposite);
            cvc5::Result res = solver.checkSat();
            solver.pop();
            if (res.isUnsat())
            {
                result.forced.emplace_back(v, modelValue);
            }
            else
            {
                candidates.push_back(v);
            }
        }
        solver.pop();
    }

    if (candidates.empty())
    {
        return result;
    }

    // ---- Pass 2: implicit definability (Padoa) over the candidates -------
    // Build phi(x) AND phi(z) where z is a fresh copy of every variable. A
    // candidate v is implicitly defined by the already-accepted support when,
    // forcing the support copies equal, v and its copy must agree -- i.e. the
    // formula admits no model with x and z agreeing on the support yet
    // differing on v.
    solver.push();

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

    cvc5::Term phiZ = formula.substitute(allVars, allCopies);
    solver.assertFormula(tm.mkTerm(cvc5::Kind::AND, {formula, phiZ}));

    std::vector<cvc5::Term> support;
    support.reserve(candidates.size());
    for (const cvc5::Term& v : candidates)
    {
        solver.push();
        for (const cvc5::Term& s : support)
        {
            cvc5::Term zs = copyMap[s.toString()];
            solver.assertFormula(tm.mkTerm(cvc5::Kind::EQUAL, {s, zs}));
        }
        cvc5::Term zv = copyMap[v.toString()];
        cvc5::Term neq = tm.mkTerm(cvc5::Kind::NOT,
                                   {tm.mkTerm(cvc5::Kind::EQUAL, {v, zv})});
        solver.assertFormula(neq);

        cvc5::Result res = solver.checkSat();
        solver.pop();

        if (res.isUnsat())
        {
            result.defined.push_back(v);
        }
        else
        {
            support.push_back(v);
        }
    }

    solver.pop();

    result.support = std::move(support);
    return result;
}

void minimizeProjectionSet(cvc5::Solver& solver,
                           const cvc5::Term& formula,
                           std::vector<cvc5::Term>& projVars)
{
    Reduction r = analyzeProjectionSet(solver, formula, projVars);
    projVars = std::move(r.support);
}

} // namespace arjun
