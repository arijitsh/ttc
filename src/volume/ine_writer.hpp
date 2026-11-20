#pragma once

// --dump-ine: a polytope in the cdd / Avis-Fukuda H-representation format,
// readable by external volume tools (vinci, lrs, cddlib).

#include <string>
#include <vector>

#include <Eigen/Dense>

namespace ttc
{

// Writes {x : A x <= b} to `path`; column j is the variable names[j].  A file
// that cannot be opened is skipped silently.
void writePolytopeIne(const std::string& path, const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                      const std::vector<std::string>& names);

}  // namespace ttc
