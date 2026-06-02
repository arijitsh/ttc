#include "OptionPartialOrderHeuristic.hpp"
#include "src/configurations/ConfigurationPartialOrderHeuristic.hpp"

namespace d4 {
OptionPartialOrderHeuristic::OptionPartialOrderHeuristic()
    : OptionPartialOrderHeuristic(ConfigurationPartialOrderHeuristic()) {}

OptionPartialOrderHeuristic::OptionPartialOrderHeuristic(
    const ConfigurationPartialOrderHeuristic& config) {
  partialOrderMethod = config.partialOrderMethod;
  givenOrder = config.givenOrder;
  scaleFactor = config.scaleFactor;
  optionTreeDecomposition = nullptr;
}

OptionPartialOrderHeuristic::~OptionPartialOrderHeuristic() {
  if (optionTreeDecomposition) delete optionTreeDecomposition;
}

std::ostream& operator<<(std::ostream& out,
                         const OptionPartialOrderHeuristic&) {
  out << " Option Partitioning Heuristic";
  return out;
}
} // namespace d4
