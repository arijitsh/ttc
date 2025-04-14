#ifndef PARSER_HPP
#define PARSER_HPP

#include <cvc5/cvc5.h>
#include <cvc5/cvc5_parser.h>
#include <string>

class TTCParser {
public:
    TTCParser();
    ~TTCParser();

    // Parse an SMT-LIB formula provided as a string.
    void parseFormula(const std::string &smtFormula);

    // Check the satisfiability of the current assertions and return the result as a string.
    std::string checkSatisfiability();

private:
    // The main solver instance.
    cvc5::Solver d_solver;

    // Pointer to the input parser (created when formula is provided).
    cvc5::parser::InputParser* d_parser;
};

#endif // TTC_PARSER_H
