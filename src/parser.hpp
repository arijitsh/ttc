#ifndef PARSER_HPP
#define PARSER_HPP

#include <cvc5/cvc5.h>
#include <cvc5/cvc5_parser.h>
#include <string>
#include <vector>
#include <cstdint>
#include <unordered_set>
#include <optional>
#include "cache_mode.hpp"
#include "features.hpp"
#include "volume/allsat_volume.hpp"

class TTCParser {
public:
    TTCParser();
    ~TTCParser();

    // Parse an SMT-LIB formula provided as a string.
    void parseFormula(const std::string &smtFormula);

    // Minimize the set of projection variables via implicit definability.
    void minimizeProjectionSet();

    // Count models projected on Boolean variables with prefix "proj_".
    std::uint64_t projectedModelCount(CacheMode cacheMode,
                                      bool simplifyProp,
                                      bool contract,
                                      bool netrel,
                                      bool useAssumptions,
                                      bool useBitset,
                                      int propAt,
                                      bool useMono,
                                      bool monoTrue);

    void printTreeDecomposition(bool contract, bool netrel);

    void computePolytopes();
    const std::vector<ttc::Polytope>& polytopes() const { return d_polytopes; }
    std::vector<cvc5::Term> realVariables() const;

    // Accessors for variable counts.
    std::size_t numBoolVars() const { return d_boolVars.size(); }
    std::size_t numIntVars() const { return d_intVars.size(); }
    std::size_t numRealVars() const { return d_realVars.size(); }
    std::size_t numBvVars() const { return d_bvVars.size(); }
    std::size_t numProjVars() const { return d_projVars.size(); }
    std::size_t projVarsBeforeMin() const { return d_projVarsBefore; }
    std::size_t projVarsAfterMin() const { return d_projVarsAfter; }
    std::size_t numConstraints() const { return d_numConstraints; }

    // Promote all Boolean and bit-vector variables to projection variables when
    // no explicit projection variables were provided.
    void promoteBooleanAndBvToProjection();

    // Access solver and projection variables for external algorithms.
    cvc5::Solver& solver() { return d_solver; }
    const std::vector<cvc5::Term>& projectionVars() const { return d_projVars; }
    const std::vector<cvc5::Term>& nonProjectionVars() const
    {
        return d_nonProjVars;
    }
    const std::vector<cvc5::Term>& assertions() const { return d_assertions; }
    const cvc5::Term& formula() const { return d_formula; }

    // If all variables are Boolean or bit-vector, check satisfiability with STP.
    // Returns std::nullopt if STP is unavailable or other sorts are present.
    std::optional<bool> checkSatWithSTP(const std::string& smtFormula);

private:
    // Term manager owning term data for the solver when available.
    ttc::TermBuilderHelper<cvc5::Solver>::storage_type d_tm;
    // The main solver instance.
    cvc5::Solver d_solver;

    // Pointer to the input parser (created when formula is provided).
    cvc5::parser::InputParser* d_parser;

    // Conjunction of all asserted formulas.
     cvc5::Term d_formula;

    // Individual assertions
    std::vector<cvc5::Term> d_assertions;

    // List of projection variables (Boolean terms starting with "proj").
    std::vector<cvc5::Term> d_projVars;
    std::vector<cvc5::Term> d_nonProjVars;
    std::size_t d_projVarsBefore = 0;
    std::size_t d_projVarsAfter = 0;

    // Sets of variables by sort.
    std::unordered_set<std::string> d_boolVars;
    std::unordered_set<std::string> d_intVars;
    std::unordered_set<std::string> d_realVars;
    std::unordered_set<std::string> d_bvVars;

    // Number of top-level constraints in the formula.
    std::size_t d_numConstraints = 0;

    std::vector<ttc::Polytope> d_polytopes;

};

#endif // TTC_PARSER_H
