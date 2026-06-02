#ifndef TOCNF_HPP
#define TOCNF_HPP

#include <cvc5/cvc5.h>
#include <vector>
#include <string>
#include <unordered_map>

struct ToCNFStats {
    std::size_t smtOps;
    std::size_t boolVars;
    std::size_t clauses;
    double timeSec;
};

class ToCNF {
public:
    ToCNF(cvc5::Solver& solver, const std::vector<cvc5::Term>& assertions);

    // Convert assertions to CNF and write to given file.
    // Returns statistics per assertion.
    std::vector<ToCNFStats> convert(const std::string& filename);

    // Build CNF representation of the assertions and return clauses
    // together with a mapping from variable name to index.  The returned
    // vector contains one clause per entry where positive integers
    // represent positive literals and negative integers represent
    // negated literals.
    struct CNFFormula
    {
        int varCount;                                  // number of variables
        std::unordered_map<std::string,int> varToIdx;  // mapping of variable names
        std::unordered_map<cvc5::Term,int> termToIdx;  // mapping of terms
        std::vector<cvc5::Term> idxToTerm;             // index to term mapping
        std::vector<std::vector<int>> clauses;         // CNF clauses
    };

    CNFFormula build();

private:
    cvc5::Solver& d_solver;
    std::vector<cvc5::Term> d_assertions;

    std::size_t countOps(const cvc5::Term& t);
};

#endif // TOCNF_HPP
