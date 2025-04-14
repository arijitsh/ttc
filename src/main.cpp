#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "parser.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  // Define command-line options
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "Display help information")(
      "solve", po::value<std::string>(), "SMT2 file to solve");

  po::variables_map vm;
  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  }
  catch (const std::exception& ex)
  {
    std::cerr << "Error parsing command line: " << ex.what() << std::endl;
    return 1;
  }

  if (vm.count("help") || !vm.count("solve"))
  {
    std::cout << "Usage: " << argv[0] << " --solve filename.smt2" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  std::string filename = vm["solve"].as<std::string>();
  std::ifstream inFile(filename);
  if (!inFile)
  {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return 1;
  }

  // Read the entire file contents into a string
  std::stringstream buffer;
  buffer << inFile.rdbuf();
  std::string smtFormula = buffer.str();

  // Create our TTCParser instance and parse the formula
  TTCParser parser;
  parser.parseFormula(smtFormula);

  // Run check-sat (if not already done by the parsed script)
  std::string result = parser.checkSatisfiability();
  std::cout << "Satisfiability result: " << result << std::endl;

  return 0;
}
