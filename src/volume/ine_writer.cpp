#include "volume/ine_writer.hpp"

#include <fstream>
#include <iomanip>

namespace ttc
{

void writePolytopeIne(const std::string& path, const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                      const std::vector<std::string>& names)
{
  std::ofstream out(path);
  if (!out) return;
  // A cdd row [b_i, -a_i1, ..., -a_in] encodes b_i - a_i.x >= 0, i.e. a_i.x <= b_i.
  out << "ttc_polytope\n";
  out << "* variables:";  // cdd treats '*' lines as comments
  for (const auto& name : names) out << " " << name;
  out << "\n";
  out << "H-representation\n";
  out << "begin\n";
  out << " " << A.rows() << " " << (A.cols() + 1) << " real\n";
  out << std::setprecision(17);
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    out << " " << b(i);
    for (Eigen::Index j = 0; j < A.cols(); ++j) out << " " << -A(i, j);
    out << "\n";
  }
  out << "end\n";
}

}  // namespace ttc
