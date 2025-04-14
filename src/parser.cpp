#include "parser.hpp"
#include <iostream>
#include <sstream>

TTCParser::TTCParser() : d_solver(), d_parser(nullptr) {
    // Set options on the solver:
    d_solver.setOption("print-success", "true");
    // d_solver.setLogic("ALL");  // You can change this as needed
}

TTCParser::~TTCParser() {
    if (d_parser != nullptr) {
        delete d_parser;
        d_parser = nullptr;
    }
}

void TTCParser::parseFormula(const std::string &smtFormula) {
    // Create an input stream from the provided formula string.
    std::istringstream iss(smtFormula);

    // Create a new parser associated with our solver.
    d_parser = new cvc5::parser::InputParser(&d_solver);

    // Set the input stream (using the SMT-LIB 2.6 language mode here)
    d_parser->setStreamInput(
        cvc5::modes::InputLanguage::SMT_LIB_2_6, iss, "input_stream");

    // Loop to fetch and execute commands from the input stream.
    cvc5::parser::Command cmd;
    while (true) {
        cmd = d_parser->nextCommand();
        if (cmd.isNull()) {
            break;
        }
        std::cout << "Executing command: " << cmd << std::endl;
        // Invoke the command on our solver; output is sent to std::cout.

        cmd.invoke(&d_solver, d_parser->getSymbolManager(), std::cout);
    }
    auto terms = d_parser->getSymbolManager()->getDeclaredTerms();

    for(auto& term : terms) {
        std::cout << "Term: " << term << " type: " << term.getSort() << std::endl;
    }

    std::cout << "Parsed terms: " << terms << std::endl;

    if (d_parser->getSolver()->getLogic() == "QF_LRA")
    {
      std::cout << "Logic is QF_LRA, turning to Volume mode" << std::endl;
    }
    else if (d_parser->getSolver()->getLogic() == "QF_LIA")
    {
      std::cout << "Logic is QF_LIA, turning to Counting mode" << std::endl;
    }
    else
    {
      std::cout << "Logic is not QF_LRA or QF_LIA" << std::endl;
    }
}

std::string TTCParser::checkSatisfiability() {
    cvc5::Result res = d_solver.checkSat();
    auto stat = d_solver.getAssertions();
    std::cout << "Assertions: " << stat << std::endl;
    std::stringstream ss;
    ss << res;
    return ss.str();
}
