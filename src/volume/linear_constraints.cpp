#include "volume/linear_constraints.hpp"

#include <optional>
#include <stdexcept>
#include <unordered_set>

#include "logger.hpp"

namespace ttc
{
namespace
{

constexpr double kZeroTolerance = 1e-9;

double parseRational(const std::string& value)
{
  if (value.empty()) return 0.0;
  const auto slash = value.find('/');
  if (slash == std::string::npos) return std::stod(value);
  const double numerator = std::stod(value.substr(0, slash));
  const double denominator = std::stod(value.substr(slash + 1));
  if (denominator == 0.0) throw std::runtime_error("encountered zero denominator in rational value");
  return numerator / denominator;
}

// The value of a variable-free arithmetic term, if it has one.
std::optional<double> constantValue(const cvc5::Term& term)
{
  if (term.isRealValue()) return parseRational(term.getRealValue());
  if (term.isIntegerValue()) return parseRational(term.getIntegerValue());
  switch (term.getKind())
  {
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER: return constantValue(term[0]);
    case cvc5::Kind::NEG:
      if (auto inner = constantValue(term[0])) return -*inner;
      return std::nullopt;
    case cvc5::Kind::SUB:
      if (term.getNumChildren() != 1) return std::nullopt;
      if (auto inner = constantValue(term[0])) return -*inner;
      return std::nullopt;
    case cvc5::Kind::DIVISION:
    {
      if (term.getNumChildren() != 2) return std::nullopt;
      const auto numerator = constantValue(term[0]);
      const auto denominator = constantValue(term[1]);
      if (numerator && denominator && *denominator != 0.0) return *numerator / *denominator;
      return std::nullopt;
    }
    case cvc5::Kind::MULT:
    {
      double product = 1.0;
      for (const auto& child : term)
      {
        const auto factor = constantValue(child);
        if (!factor) return std::nullopt;
        product *= *factor;
      }
      if (term.getNumChildren() == 0) return std::nullopt;
      return product;
    }
    default: return std::nullopt;
  }
}

// a.x + constant, the linear form of an arithmetic term.
struct LinearForm
{
  Eigen::VectorXd a;
  double constant = 0.0;
};

// Adds multiplier * term to `form`; false if term is not linear in the variables.
bool accumulate(const cvc5::Term& term, double multiplier, const RealVariables& variables,
                LinearForm& form)
{
  if (multiplier == 0.0) return true;
  if (term.isRealValue())
  {
    form.constant += multiplier * parseRational(term.getRealValue());
    return true;
  }
  if (term.isIntegerValue())
  {
    form.constant += multiplier * parseRational(term.getIntegerValue());
    return true;
  }
  if (term.getNumChildren() == 0)
  {
    if (!term.getSort().isReal() && !term.getSort().isInteger()) return false;
    const long column = variables.column(term);
    if (column < 0)
    {
      Log(2) << "Unknown variable in linear term: " << term << std::endl;
      return false;
    }
    form.a[column] += multiplier;
    return true;
  }

  switch (term.getKind())
  {
    case cvc5::Kind::ADD:
      for (const auto& child : term)
        if (!accumulate(child, multiplier, variables, form)) return false;
      return true;
    case cvc5::Kind::SUB:
      if (term.getNumChildren() == 1) return accumulate(term[0], -multiplier, variables, form);
      if (!accumulate(term[0], multiplier, variables, form)) return false;
      for (std::size_t i = 1; i < term.getNumChildren(); ++i)
        if (!accumulate(term[i], -multiplier, variables, form)) return false;
      return true;
    case cvc5::Kind::NEG: return accumulate(term[0], -multiplier, variables, form);
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER: return accumulate(term[0], multiplier, variables, form);
    case cvc5::Kind::DIVISION:
    {
      if (term.getNumChildren() != 2) return false;
      const auto divisor = constantValue(term[1]);
      if (!divisor || *divisor == 0.0) return false;
      return accumulate(term[0], multiplier / *divisor, variables, form);
    }
    case cvc5::Kind::MULT:
    {
      // All factors but (at most) one must be constant.
      double factor = 1.0;
      std::optional<cvc5::Term> linearFactor;
      for (const auto& child : term)
      {
        if (const auto value = constantValue(child))
          factor *= *value;
        else if (linearFactor)
          return false;
        else
          linearFactor = child;
      }
      if (!linearFactor)
      {
        form.constant += multiplier * factor;
        return true;
      }
      return accumulate(*linearFactor, multiplier * factor, variables, form);
    }
    default: return false;
  }
}

// The rows of A x <= b, gathered literal by literal.
class RowCollector
{
 public:
  explicit RowCollector(const RealVariables& variables) : d_variables(variables) {}

  // Adds the rows of `term` (or of its negation if !positive).  False if it is
  // not a linear comparison, a conjunction of them or a Boolean constant.
  bool add(const cvc5::Term& term, bool positive)
  {
    if (d_empty) return true;
    if (term.getKind() == cvc5::Kind::NOT) return add(term[0], !positive);
    if (term.isBooleanValue())
    {
      if (term.getBooleanValue() != positive) d_empty = true;
      return true;
    }
    if (term.getKind() == cvc5::Kind::AND)
    {
      for (const auto& child : term)
        if (!add(child, positive)) return false;
      return true;
    }
    return addComparison(term, positive);
  }

  bool empty() const { return d_empty; }
  bool hasRows() const { return !d_rows.empty(); }

  LinearConstraints matrices() const
  {
    const long n = static_cast<long>(d_variables.size());
    const long m = static_cast<long>(d_rows.size());
    LinearConstraints P;
    P.A = Eigen::MatrixXd(m, n);
    P.b = Eigen::VectorXd(m);
    for (long i = 0; i < m; ++i)
    {
      P.A.row(i) = d_rows[i].a.transpose();
      P.b[i] = d_rows[i].bound;
    }
    return P;
  }

 private:
  struct Row
  {
    Eigen::VectorXd a;
    double bound = 0.0;  // a.x <= bound
  };

  bool addComparison(const cvc5::Term& atom, bool positive)
  {
    cvc5::Kind kind = atom.getKind();
    if (!positive)
    {
      switch (kind)
      {
        case cvc5::Kind::LEQ: kind = cvc5::Kind::GEQ; break;
        case cvc5::Kind::LT: kind = cvc5::Kind::GT; break;
        case cvc5::Kind::GEQ: kind = cvc5::Kind::LEQ; break;
        case cvc5::Kind::GT: kind = cvc5::Kind::LT; break;
        case cvc5::Kind::EQUAL: return true;  // a != b: a null set, no row
        default: break;
      }
    }
    switch (kind)
    {
      case cvc5::Kind::LEQ:
      case cvc5::Kind::LT: return addLessEqual(atom[0], atom[1]);
      case cvc5::Kind::GEQ:
      case cvc5::Kind::GT: return addLessEqual(atom[1], atom[0]);
      case cvc5::Kind::EQUAL: return addLessEqual(atom[0], atom[1]) && addLessEqual(atom[1], atom[0]);
      default: return false;
    }
  }

  // lhs <= rhs as (lhs - rhs).x <= -(constant of lhs - rhs).
  bool addLessEqual(const cvc5::Term& lhs, const cvc5::Term& rhs)
  {
    LinearForm form{Eigen::VectorXd::Zero(static_cast<long>(d_variables.size())), 0.0};
    if (!accumulate(lhs, 1.0, d_variables, form) || !accumulate(rhs, -1.0, d_variables, form))
      return false;
    if (form.a.norm() <= kZeroTolerance)
    {
      // A comparison of constants: true adds nothing, false empties P(c).
      if (form.constant > kZeroTolerance) d_empty = true;
      return true;
    }
    d_rows.push_back({std::move(form.a), -form.constant});
    return true;
  }

  const RealVariables& d_variables;
  std::vector<Row> d_rows;
  bool d_empty = false;
};

}  // namespace

RealVariables::RealVariables(const std::vector<cvc5::Term>& variables)
{
  std::unordered_set<cvc5::Term> seen;
  for (const auto& term : variables)
  {
    if (!term.getSort().isReal() || !seen.insert(term).second) continue;
    d_column.emplace(term, d_terms.size());
    d_terms.push_back(term);
    d_names.push_back(term.hasSymbol() ? term.getSymbol() : term.toString());
  }
}

long RealVariables::column(const cvc5::Term& term) const
{
  const auto it = d_column.find(term);
  return it == d_column.end() ? -1 : static_cast<long>(it->second);
}

LinearConstraints cubeConstraints(const Cube& cube, const RealVariables& variables)
{
  RowCollector rows(variables);
  for (const TermLiteral& literal : cube.literals)
  {
    if (!rows.add(literal.term, literal.positive))
      throw std::runtime_error("unsupported literal encountered during LRA volume conversion");
    if (rows.empty()) break;
  }
  if (rows.empty())
  {
    LinearConstraints P;
    P.A = Eigen::MatrixXd::Zero(0, static_cast<long>(variables.size()));
    P.b = Eigen::VectorXd::Zero(0);
    P.empty = true;
    return P;
  }
  if (!rows.hasRows()) throw std::runtime_error("polytope has no linear constraints after abstraction");
  return rows.matrices();
}

}  // namespace ttc
