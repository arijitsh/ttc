#include "parser.hpp"

#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "cache_mode.hpp"
#include "features.hpp"
#include "logger.hpp"
#include "stp_sat.hpp"
#include "volume/allsat_volume.hpp"

#ifdef TTC_ENABLE_DDNNF
#include "arjun.hpp"
#include "proj_ddnnf.hpp"
#include "var_order.hpp"
#endif
TTCParser::TTCParser()
    : d_tm(),
      d_solver(ttc::createSolverWithStorage<cvc5::Solver>(d_tm)),
      d_parser(nullptr),
      d_formula(),
      d_projVars(),
      d_boolVars(),
      d_intVars(),
      d_realVars(),
      d_numConstraints(0),
      d_polytopes()
{
    d_solver.setOption("print-success", "false");
    d_solver.setOption("incremental", "true");
    d_solver.setOption("produce-models", "true");
    d_solver.setOption("produce-learned-literals", "true");
}

TTCParser::~TTCParser() {
    if (d_parser != nullptr) {
        delete d_parser;
        d_parser = nullptr;
    }
}

void TTCParser::parseFormula(const std::string &smtFormula) {
    Log(3) << "Starting parse" << std::endl;
    std::istringstream iss(smtFormula);
    d_parser = new cvc5::parser::InputParser(&d_solver);
    d_parser->setStreamInput(cvc5::modes::InputLanguage::SMT_LIB_2_6, iss, "input_stream");

    cvc5::parser::Command cmd;
    while (true) {
        cmd = d_parser->nextCommand();
        if (cmd.isNull()) {
            break;
        }
        std::stringstream out;
        cmd.invoke(&d_solver, d_parser->getSymbolManager(), out);
    }
    Log(3) << "Commands processed" << std::endl;
    auto assertions = d_solver.getAssertions();
    d_assertions = assertions;
    d_projVars.clear();
    d_nonProjVars.clear();
    d_boolVars.clear();
    d_intVars.clear();
    d_realVars.clear();
    d_bvVars.clear();
    d_polytopes.clear();
    if (!assertions.empty())
    {
        auto& tm = ttc::getTermBuilder(d_solver);
        cvc5::Term rawFormula = assertions.size() == 1
                                    ? assertions[0]
                                    : tm.mkTerm(cvc5::Kind::AND, assertions);
        d_formula = d_solver.simplify(rawFormula);

        // Count top-level constraints recursively.
        std::function<std::size_t(const cvc5::Term&)> countConstraints = [&](const cvc5::Term& t) {
            if (t.getKind() == cvc5::Kind::AND)
            {
                std::size_t sum = 0;
                for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
                {
                    sum += countConstraints(t[i]);
                }
                return sum;
            }
            return static_cast<std::size_t>(1);
        };
        d_numConstraints = countConstraints(d_formula);

        std::unordered_set<std::string> seen;
        std::function<void(const cvc5::Term&)> collect = [&](const cvc5::Term& t) {
            if (t.getNumChildren() == 0 && t.hasSymbol())
            {
                std::string name = t.toString();
                if (name != "true" && name != "false" && seen.insert(name).second)
                {
                    if (t.getSort().isBoolean())
                    {
                        d_boolVars.insert(name);
                        Log(4) << "Boolean variable: " << name << std::endl;
                        if (name.rfind("proj_", 0) == 0)
                        {
                            d_projVars.push_back(t);
                        }
                        else
                        {
                            d_nonProjVars.push_back(t);
                        }
                    }
                    else if (t.getSort().isInteger())
                    {
                        d_intVars.insert(name);
                        d_nonProjVars.push_back(t);
                    }
                    else if (t.getSort().isReal())
                    {
                        d_realVars.insert(name);
                        d_nonProjVars.push_back(t);
                    }
                    else if (t.getSort().isBitVector())
                    {
                        d_bvVars.insert(name);
                        d_nonProjVars.push_back(t);
                    }
                }
            }
            for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
            {
                collect(t[i]);
            }
        };
        // Collect projection-variable candidates from the unsimplified formula
        // so that variables which survive in the asserted constraints but are
        // rewritten away by simplification (e.g. when the formula reduces to a
        // tautology over them) are still recognised as part of the support.
        collect(rawFormula);
    }
}

void TTCParser::promoteBooleanAndBvToProjection()
{
    if (d_nonProjVars.empty())
    {
        return;
    }
    std::vector<cvc5::Term> remaining;
    remaining.reserve(d_nonProjVars.size());
    for (const cvc5::Term& t : d_nonProjVars)
    {
        cvc5::Sort sort = t.getSort();
        if (sort.isBoolean() || sort.isBitVector())
        {
            d_projVars.push_back(t);
        }
        else
        {
            remaining.push_back(t);
        }
    }
    d_nonProjVars.swap(remaining);
}

#ifdef TTC_ENABLE_DDNNF
std::uint64_t TTCParser::projectedModelCount(CacheMode cacheMode,
                                             bool simplifyProp,
                                             bool contract,
                                             bool netrel,
                                             bool useAssumptions,
                                             bool useBitset,
                                             int propAt,
                                             bool useMono,
                                             bool monoTrue)
{
    Log(3) << "Starting model count with DDNNF" << std::endl;
    ProjDDNNF ddnnf(d_solver,
                    d_formula,
                    d_projVars,
                    cacheMode,
                    simplifyProp,
                    contract,
                    netrel,
                    useAssumptions,
                    useBitset,
                    propAt,
                    useMono,
                    monoTrue);
    return ddnnf.count();
}

void TTCParser::printTreeDecomposition(bool contract, bool netrel)
{
    computeProjVarOrder(d_formula, d_projVars, d_solver, true, contract, netrel);
}

void TTCParser::minimizeProjectionSet()
{
    d_projVarsBefore = d_projVars.size();
    Log(3) << "Minimizing projection set of size " << d_projVarsBefore << std::endl;
    arjun::minimizeProjectionSet(d_solver, d_formula, d_projVars);
    d_projVarsAfter = d_projVars.size();
    Log(3) << "Projection set reduced to " << d_projVarsAfter << std::endl;
}
#else
std::uint64_t TTCParser::projectedModelCount(CacheMode,
                                             bool,
                                             bool,
                                             bool,
                                             bool,
                                             bool,
                                             int,
                                             bool,
                                             bool)
{
    throw std::runtime_error("DDNNF-based counting is not enabled in this build");
}

void TTCParser::printTreeDecomposition(bool, bool)
{
    Log(2) << "Tree decomposition requested but DDNNF support is disabled" << std::endl;
}

void TTCParser::minimizeProjectionSet()
{
    Log(2) << "Projection minimization requested but DDNNF support is disabled" << std::endl;
    d_projVarsBefore = d_projVars.size();
    d_projVarsAfter = d_projVars.size();
}
#endif

std::optional<bool> TTCParser::checkSatWithSTP(const std::string& smtFormula)
{
    // Require only Boolean and bit-vector variables.
    if (!d_intVars.empty() || !d_realVars.empty())
    {
        return std::nullopt;
    }

    return csb::stpCheckSat(smtFormula);
}

void TTCParser::computePolytopes()
{
    d_polytopes.clear();
    if (d_formula.isNull())
    {
        return;
    }

    try
    {
        auto abstraction = ttc::getBooleanAbstractionAig(d_solver, d_formula);
        auto result = ttc::enumeratePolytopes(abstraction, d_solver);
        d_polytopes = std::move(result.polytopes);
    }
    catch (const std::exception&)
    {
        try
        {
            auto manual = ttc::buildBooleanAigFromTerm(d_formula);
            auto result = ttc::enumeratePolytopes(manual, d_solver);
            d_polytopes = std::move(result.polytopes);
        }
        catch (const std::exception&)
        {
            d_polytopes.clear();
        }
    }
}

std::vector<cvc5::Term> TTCParser::realVariables() const
{
    std::unordered_set<cvc5::Term> seen;
    std::vector<cvc5::Term> vars;
    vars.reserve(d_nonProjVars.size() + d_projVars.size());
    for (const auto& term : d_nonProjVars)
    {
        if (term.getSort().isReal() && seen.insert(term).second)
        {
            vars.push_back(term);
        }
    }
    for (const auto& term : d_projVars)
    {
        if (term.getSort().isReal() && seen.insert(term).second)
        {
            vars.push_back(term);
        }
    }
    return vars;
}
