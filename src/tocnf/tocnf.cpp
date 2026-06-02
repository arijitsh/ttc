#include "tocnf.hpp"

#include "features.hpp"

#include <chrono>
#include <cstdint>
#include <fstream>
#include <utility>

namespace {

cvc5::Term preprocessTerm(cvc5::Solver& solver, const cvc5::Term& term)
{
    if (term.isNull())
    {
        return term;
    }
    return solver.simplify(term);
}

std::size_t decodeClauseData(const std::vector<uint32_t>& cnfData,
                             std::vector<std::vector<int>>* clauses,
                             int* maxVar)
{
    if (clauses)
    {
        clauses->clear();
    }

    std::vector<int> clause;
    std::size_t clauseCount = 0;
    bool clauseStarted = false;
    int localMax = 0;

    for (uint32_t raw : cnfData)
    {
        int lit = static_cast<int32_t>(raw);
        if (lit == 0)
        {
            if (clauses)
            {
                clauses->push_back(clause);
            }
            clause.clear();
            clauseStarted = false;
            ++clauseCount;
            continue;
        }

        clauseStarted = true;
        int absLit = lit < 0 ? -lit : lit;
        if (absLit > localMax)
        {
            localMax = absLit;
        }
        if (clauses)
        {
            clause.push_back(lit);
        }
    }

    if (clauseStarted)
    {
        if (clauses)
        {
            clauses->push_back(clause);
        }
        ++clauseCount;
    }

    if (maxVar)
    {
        *maxVar = localMax;
    }

    return clauseCount;
}

} // namespace

ToCNF::ToCNF(cvc5::Solver& solver, const std::vector<cvc5::Term>& assertions)
    : d_solver(solver)
{
    d_assertions.reserve(assertions.size());
    for (const auto& term : assertions)
    {
        d_assertions.push_back(preprocessTerm(d_solver, term));
    }
}

std::size_t ToCNF::countOps(const cvc5::Term& t)
{
    std::size_t sum = 1;
    for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
    {
        sum += countOps(t[i]);
    }
    return sum;
}

std::vector<ToCNFStats> ToCNF::convert(const std::string& filename)
{
    std::vector<ToCNFStats> stats;
    stats.reserve(d_assertions.size());

    for (const auto& assertion : d_assertions)
    {
        auto start = std::chrono::steady_clock::now();
        auto abstraction = ttc::getBooleanAbstraction(d_solver, assertion);
        auto end = std::chrono::steady_clock::now();
        double sec = std::chrono::duration<double>(end - start).count();

        std::size_t clauseCount = decodeClauseData(abstraction.clauses, nullptr, nullptr);
        ToCNFStats st{countOps(assertion), abstraction.mapping.size(), clauseCount, sec};
        stats.push_back(st);
    }

    CNFFormula cnf = build();
    std::ofstream out(filename);
    out << "p cnf " << cnf.varCount << " " << cnf.clauses.size() << "\n";
    out << "c p show";
    for (int i = 1; i <= cnf.varCount; ++i)
    {
        out << " " << i;
    }
    out << " 0\n";
    for (const auto& cl : cnf.clauses)
    {
        for (int lit : cl) out << lit << ' ';
        out << "0\n";
    }
    return stats;
}

ToCNF::CNFFormula ToCNF::build()
{
    cvc5::Term formula;
    auto& tm = ttc::getTermBuilder(d_solver);
    if (d_assertions.size() == 1)
    {
        formula = d_assertions[0];
    }
    else if (!d_assertions.empty())
    {
        formula = tm.mkTerm(cvc5::Kind::AND, d_assertions);
    }
    else
    {
        formula = tm.mkBoolean(true);
    }

    formula = preprocessTerm(d_solver, formula);

    auto abstraction = ttc::getBooleanAbstraction(d_solver, formula);

    int maxVar = 0;
    std::vector<std::vector<int>> clauses;
    decodeClauseData(abstraction.clauses, &clauses, &maxVar);

    CNFFormula res;
    res.varCount = maxVar;
    res.clauses = std::move(clauses);
    res.idxToTerm.assign(res.varCount + 1, cvc5::Term());

    for (const auto& [idx, term] : abstraction.mapping)
    {
        int id = static_cast<int>(idx);
        res.varToIdx[term.toString()] = id;
        res.termToIdx.emplace(term, id);

        if (id > res.varCount)
        {
            res.varCount = id;
            if (static_cast<int>(res.idxToTerm.size()) <= res.varCount)
            {
                res.idxToTerm.resize(res.varCount + 1);
            }
        }
        else if (static_cast<int>(res.idxToTerm.size()) <= id)
        {
            res.idxToTerm.resize(id + 1);
        }

        res.idxToTerm[id] = term;
    }

    if (static_cast<int>(res.idxToTerm.size()) <= res.varCount)
    {
        res.idxToTerm.resize(res.varCount + 1);
    }

    return res;
}

