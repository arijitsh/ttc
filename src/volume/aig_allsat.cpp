#include "volume/aig_allsat.hpp"

#include <iostream>

#include "Aiger/AigerParser.hpp"
#include "AllSatAlgo/Blocking/TseitinEnc/AllSatAlgoTseitinEnc.hpp"
#include "Globals/AllSatSolverGloblas.hpp"
#include "Globals/TernaryVal.hpp"
#include "Utilities/InputParser.hpp"

namespace ttc
{
namespace
{

// The library reads AIGs from AIGER files; this loads one from memory.
class InMemoryAigerParser : public AigerParser
{
 public:
  explicit InMemoryAigerParser(const BooleanAbstractionAigData& aig)
  {
    m_Inputs.assign(aig.inputs.begin(), aig.inputs.end());
    m_Outputs.assign(aig.outputs.begin(), aig.outputs.end());
    m_AndGates.reserve(aig.andGates.size());
    m_IsVarRef.assign(static_cast<std::size_t>(aig.maxVariable) + 1, false);
    if (!m_IsVarRef.empty()) m_IsVarRef[0] = true;
    auto markVar = [&](std::uint32_t literal) {
      const std::uint32_t var = literal >> 1U;
      if (var >= m_IsVarRef.size()) m_IsVarRef.resize(var + 1, false);
      m_IsVarRef[var] = true;
    };
    for (std::uint32_t literal : aig.inputs) markVar(literal);
    for (const auto& gate : aig.andGates)
    {
      markVar(gate.lhs);
      markVar(gate.rhs0);
      markVar(gate.rhs1);
      m_AndGates.emplace_back(gate.lhs, gate.rhs0, gate.rhs1);
    }
  }
};

// The Tseitin-encoding blocking algorithm ("roc" mode), driven step by step
// so the generalized models are returned instead of printed.
class CubeEnumerator : public AllSatAlgoTseitinEnc
{
 public:
  CubeEnumerator(const BooleanAbstractionAigData& aig, const InputParser& params)
      : AllSatAlgoTseitinEnc(params)
  {
    m_AigParser = InMemoryAigerParser(aig);
    delete m_CirSimulation;
    m_CirSimulation = nullptr;
    if (m_UseCirSim)
      m_CirSimulation = new CirSim(m_AigParser,
                                   m_UseTopToBotSim ? SimStrat::TopToBot : SimStrat::BotToTop);
    m_Inputs = m_AigParser.GetInputs();
    m_InputSize = m_Inputs.size();
    m_Solver->InitializeSolver(m_AigParser);
    if (m_UseDualSolver && m_DualSolver != nullptr) m_DualSolver->InitializeSolver(m_AigParser);
  }

  std::vector<INPUT_ASSIGNMENT> enumerate()
  {
    std::vector<INPUT_ASSIGNMENT> cubes;
    while (m_Solver->Solve() == SAT_RET_STATUS)
    {
      INPUT_ASSIGNMENT model = m_Solver->GetAssignmentForAIGLits(m_Inputs);
      INPUT_ASSIGNMENT cube = GeneralizeModel(model);
      if (m_IsTimeOut) break;
      cubes.push_back(cube);
      BlockModel(cube);
    }
    return cubes;
  }
};

InputCube toInputCube(const INPUT_ASSIGNMENT& assignment)
{
  InputCube cube;
  cube.reserve(assignment.size());
  for (const auto& [literal, value] : assignment)
  {
    const std::uint32_t var = literal >> 1U;
    switch (value)
    {
      case TVal::True: cube.emplace_back(var, true); break;
      case TVal::False: cube.emplace_back(var, false); break;
      case TVal::DontCare: cube.emplace_back(var, std::nullopt); break;
      default: break;  // never produced for a generalized model's inputs
    }
  }
  return cube;
}

}  // namespace

std::vector<InputCube> allSatInputCubes(const BooleanAbstractionAigData& aig)
{
  if (aig.outputs.empty()) return {};

  int argc = 1;
  char arg0[] = "allsat";
  char* argv[] = {arg0};
  InputParser params(argc, argv);
  params.AppendParams({"/general/print_enumer", "0",
                       "/mode", "roc",
                       "/alg/blocking/use_ucore", "1",
                       "/alg/blocking/use_ipasir_for_dual", "0"});

  std::cout << "c starting DNFization" << std::endl;
  CubeEnumerator enumerator(aig, params);
  const std::vector<INPUT_ASSIGNMENT> assignments = enumerator.enumerate();
  std::cout << "c done DNFization" << std::endl;

  std::vector<InputCube> cubes;
  cubes.reserve(assignments.size());
  for (const auto& assignment : assignments) cubes.push_back(toInputCube(assignment));
  return cubes;
}

}  // namespace ttc
