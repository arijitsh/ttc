#pragma once
#include <optional>
#include <string>

namespace csb {
std::optional<bool> stpCheckSat(const std::string& smtFormula);
}

