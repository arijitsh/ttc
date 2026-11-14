#include "volume/lra_volume.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <list>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

#include <random>

// GMP arbitrary-precision floats for the sampling / union bookkeeping. The
// default (--gmp) path stores points and evaluates membership tests in mpf at a
// arithmetic precision derived from lattice precision b; --fullgmp additionally runs
// the selected walk itself in this type (see gmp_walks.hpp).
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/eigen.hpp>

#include "logger.hpp"
#include "volume/gmp_walks.hpp"
#include "volume/polytope_sampler.hpp"
#include "volume/polytope_volume.hpp"

// cddlib (double precision) for polytope canonicalization, mirroring the
// Python prototype's use of pycddlib (src/polytope_operations.py::canonicalize).
#define set_int cdd_set_int
extern "C" {
#include <cddlib/setoper.h>
#include <cddlib/cdd.h>
}
#undef set_int

namespace ttc
{
namespace
{

using TermIndexMap = std::unordered_map<cvc5::Term, std::size_t>;

// MpFloat (GMP float for the sampling / union bookkeeping) is defined in
// gmp_walks.hpp; its runtime precision is set globally via
// MpFloat::default_precision(digits) before sampling starts.

constexpr double kZeroTolerance = 1e-9;

// Write the polytope {x : A x <= b} to `path` in the cdd / Avis-Fukuda
// H-representation format (readable by vinci, lrs, cddlib). Each cdd row is
// [b_i, -a_i0, -a_i1, ...], encoding b_i - a_i . x >= 0, i.e. a_i . x <= b_i.
void writePolytopeIne(const std::string& path, const Eigen::MatrixXd& A,
                      const Eigen::VectorXd& b)
{
  std::ofstream out(path);
  if (!out)
  {
    return;
  }
  const Eigen::Index rows = A.rows();
  const Eigen::Index cols = A.cols();
  out << "ttc_polytope\n";
  out << "H-representation\n";
  out << "begin\n";
  out << " " << rows << " " << (cols + 1) << " real\n";
  out << std::setprecision(17);
  for (Eigen::Index i = 0; i < rows; ++i)
  {
    out << " " << b(i);
    for (Eigen::Index j = 0; j < cols; ++j)
    {
      out << " " << -A(i, j);
    }
    out << "\n";
  }
  out << "end\n";
}

double parseRationalString(const std::string& value)
{
  if (value.empty())
  {
    return 0.0;
  }
  auto slash = value.find('/');
  if (slash == std::string::npos)
  {
    return std::stod(value);
  }
  double numerator = std::stod(value.substr(0, slash));
  double denominator = std::stod(value.substr(slash + 1));
  if (denominator == 0.0)
  {
    throw std::runtime_error("encountered zero denominator in rational value");
  }
  return numerator / denominator;
}

std::optional<double> getConstantValue(const cvc5::Term& term)
{
  if (term.isRealValue())
  {
    return parseRationalString(term.getRealValue());
  }
  if (term.isIntegerValue())
  {
    return parseRationalString(term.getIntegerValue());
  }
  switch (term.getKind())
  {
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER:
      return getConstantValue(term[0]);
    case cvc5::Kind::NEG:
      if (auto inner = getConstantValue(term[0]))
      {
        return -(*inner);
      }
      return std::nullopt;
    case cvc5::Kind::SUB:
      if (term.getNumChildren() == 1)
      {
        if (auto inner = getConstantValue(term[0]))
        {
          return -(*inner);
        }
        return std::nullopt;
      }
      break;
    case cvc5::Kind::DIVISION:
      if (term.getNumChildren() == 2)
      {
        auto numerator = getConstantValue(term[0]);
        auto denominator = getConstantValue(term[1]);
        if (numerator && denominator && *denominator != 0.0)
        {
          return *numerator / *denominator;
        }
        return std::nullopt;
      }
      break;
    case cvc5::Kind::MULT:
    {
      double product = 1.0;
      bool hasFactor = false;
      for (const auto& child : term)
      {
        auto factor = getConstantValue(child);
        if (!factor)
        {
          return std::nullopt;
        }
        product *= *factor;
        hasFactor = true;
      }
      if (hasFactor)
      {
        return product;
      }
      return std::nullopt;
    }
    default:
      break;
  }
  return std::nullopt;
}

bool accumulateLinear(const cvc5::Term& term,
                      double multiplier,
                      const TermIndexMap& index,
                      Eigen::VectorXd& coefficients,
                      double& constant)
{
  if (multiplier == 0.0)
  {
    return true;
  }

  if (term.isRealValue())
  {
    constant += multiplier * parseRationalString(term.getRealValue());
    return true;
  }
  if (term.isIntegerValue())
  {
    constant += multiplier * parseRationalString(term.getIntegerValue());
    return true;
  }

  if (term.getNumChildren() == 0)
  {
    if (term.getSort().isReal() || term.getSort().isInteger())
    {
      auto it = index.find(term);
      if (it == index.end())
      {
        Log(2) << "Unknown variable in linear term: " << term << std::endl;
        return false;
      }
      coefficients[it->second] += multiplier;
      return true;
    }
    return false;
  }

  switch (term.getKind())
  {
    case cvc5::Kind::ADD:
      for (const auto& child : term)
      {
        if (!accumulateLinear(child, multiplier, index, coefficients, constant))
        {
          return false;
        }
      }
      return true;
    case cvc5::Kind::SUB:
      if (term.getNumChildren() == 1)
      {
        return accumulateLinear(term[0], -multiplier, index, coefficients, constant);
      }
      if (!accumulateLinear(term[0], multiplier, index, coefficients, constant))
      {
        return false;
      }
      for (std::size_t i = 1, n = term.getNumChildren(); i < n; ++i)
      {
        if (!accumulateLinear(term[i], -multiplier, index, coefficients, constant))
        {
          return false;
        }
      }
      return true;
    case cvc5::Kind::NEG:
      return accumulateLinear(term[0], -multiplier, index, coefficients, constant);
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER:
      return accumulateLinear(term[0], multiplier, index, coefficients, constant);
    case cvc5::Kind::DIVISION:
      if (term.getNumChildren() == 2)
      {
        auto denom = getConstantValue(term[1]);
        if (!denom || *denom == 0.0)
        {
          return false;
        }
        return accumulateLinear(term[0], multiplier / *denom, index, coefficients,
                                constant);
      }
      return false;
    case cvc5::Kind::MULT:
    {
      double factor = 1.0;
      std::optional<cvc5::Term> expr;
      for (const auto& child : term)
      {
        auto constantValue = getConstantValue(child);
        if (constantValue)
        {
          factor *= *constantValue;
          continue;
        }
        if (expr.has_value())
        {
          return false;
        }
        expr = child;
      }
      if (!expr.has_value())
      {
        constant += multiplier * factor;
        return true;
      }
      return accumulateLinear(expr.value(), multiplier * factor, index,
                              coefficients, constant);
    }
    default:
      break;
  }
  return false;
}

struct Inequality
{
  Eigen::VectorXd coefficients;
  double bound = 0.0;  // Ax <= bound
};

struct PolytopeMatrices
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  bool empty = false;
};

bool appendLessEqualConstraint(const cvc5::Term& lhs,
                               const cvc5::Term& rhs,
                               const TermIndexMap& index,
                               std::size_t dimension,
                               std::vector<Inequality>& inequalities,
                               bool& empty)
{
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(static_cast<long>(dimension));
  double constant = 0.0;
  if (!accumulateLinear(lhs, 1.0, index, coeffs, constant) ||
      !accumulateLinear(rhs, -1.0, index, coeffs, constant))
  {
    return false;
  }
  double norm = coeffs.norm();
  if (norm <= kZeroTolerance)
  {
    if (constant <= kZeroTolerance)
    {
      return true;
    }
    empty = true;
    return true;
  }
  inequalities.push_back({std::move(coeffs), -constant});
  return true;
}

bool appendComparisonConstraints(const cvc5::Term& atom,
                                 bool positive,
                                 const TermIndexMap& index,
                                 std::size_t dimension,
                                 std::vector<Inequality>& inequalities,
                                 bool& empty)
{
  cvc5::Kind kind = atom.getKind();
  if (!positive)
  {
    switch (kind)
    {
      case cvc5::Kind::LEQ:
        kind = cvc5::Kind::GEQ;
        break;
      case cvc5::Kind::LT:
        kind = cvc5::Kind::GT;
        break;
      case cvc5::Kind::GEQ:
        kind = cvc5::Kind::LEQ;
        break;
      case cvc5::Kind::GT:
        kind = cvc5::Kind::LT;
        break;
      case cvc5::Kind::EQUAL:
        throw std::runtime_error("negated equality is not supported for LRA volume computation");
      default:
        break;
    }
  }

  switch (kind)
  {
    case cvc5::Kind::LEQ:
    case cvc5::Kind::LT:
      return appendLessEqualConstraint(atom[0], atom[1], index, dimension,
                                       inequalities, empty);
    case cvc5::Kind::GEQ:
    case cvc5::Kind::GT:
      return appendLessEqualConstraint(atom[1], atom[0], index, dimension,
                                       inequalities, empty);
    case cvc5::Kind::EQUAL:
      if (!appendLessEqualConstraint(atom[0], atom[1], index, dimension,
                                     inequalities, empty))
      {
        return false;
      }
      return appendLessEqualConstraint(atom[1], atom[0], index, dimension,
                                       inequalities, empty);
    default:
      break;
  }
  return false;
}

bool collectConstraints(const cvc5::Term& term,
                        bool positive,
                        const TermIndexMap& index,
                        std::size_t dimension,
                        std::vector<Inequality>& inequalities,
                        bool& empty)
{
  if (empty)
  {
    return true;
  }

  if (term.getKind() == cvc5::Kind::NOT)
  {
    return collectConstraints(term[0], !positive, index, dimension, inequalities,
                              empty);
  }

  if (term.isBooleanValue())
  {
    bool value = term.getBooleanValue();
    bool satisfied = positive ? value : !value;
    if (!satisfied)
    {
      empty = true;
    }
    return true;
  }

  if (term.getKind() == cvc5::Kind::AND)
  {
    for (const auto& child : term)
    {
      if (!collectConstraints(child, positive, index, dimension, inequalities,
                              empty))
      {
        return false;
      }
    }
    return true;
  }

  return appendComparisonConstraints(term, positive, index, dimension,
                                     inequalities, empty);
}

PolytopeMatrices buildPolytopeMatrices(const Polytope& poly,
                                       const TermIndexMap& index,
                                       std::size_t dimension)
{
  PolytopeMatrices matrices;
  std::vector<Inequality> inequalities;
  inequalities.reserve(poly.termLiterals.size());
  bool empty = false;

  for (const auto& literal : poly.termLiterals)
  {
    if (!literal.value.has_value())
    {
      continue;
    }
    if (!collectConstraints(literal.term, literal.value.value(), index, dimension,
                            inequalities, empty))
    {
      throw std::runtime_error("unsupported literal encountered during LRA volume conversion");
    }
    if (empty)
    {
      break;
    }
  }

  if (empty)
  {
    matrices.empty = true;
    matrices.A = Eigen::MatrixXd::Zero(0, static_cast<long>(dimension));
    matrices.b = Eigen::VectorXd::Zero(0);
    return matrices;
  }

  if (inequalities.empty())
  {
    throw std::runtime_error("polytope has no linear constraints after abstraction");
  }

  matrices.A = Eigen::MatrixXd(static_cast<long>(inequalities.size()),
                               static_cast<long>(dimension));
  matrices.b = Eigen::VectorXd(static_cast<long>(inequalities.size()));

  for (std::size_t i = 0; i < inequalities.size(); ++i)
  {
    matrices.A.row(static_cast<long>(i)) = inequalities[i].coefficients.transpose();
    matrices.b[static_cast<long>(i)] = inequalities[i].bound;
  }

  return matrices;
}

// Canonicalizes the polytope {x : A x <= b} with cddlib, mirroring the Python
// prototype (src/polytope_operations.py::canonicalize -> cdd.matrix_canonicalize).
// Redundant inequalities are removed in place.  Returns false when the polytope
// is empty or has implicit equalities (i.e. it is lower-dimensional, so its
// full-dimensional volume is zero) -- the same cases where the Python code
// returns -1 / skips the polytope.  On any cddlib error the matrices are left
// unchanged and true is returned (fall back to the raw representation).
bool canonicalizePolytope(Eigen::MatrixXd& A, Eigen::VectorXd& b)
{
  const long m = A.rows();
  const long n = A.cols();
  if (m == 0)
  {
    return true;
  }

  dd_set_global_constants();

  // cddlib H-representation: row (g_0, g_1, ..., g_n) means
  // g_0 + g_1 x_1 + ... + g_n x_n >= 0.  For A x <= b that is b - A x >= 0,
  // so the row is (b_i, -A_i1, ..., -A_in).
  dd_MatrixPtr M = dd_CreateMatrix(m, n + 1);
  M->representation = dd_Inequality;
  for (long i = 0; i < m; ++i)
  {
    dd_set_d(M->matrix[i][0], b(i));
    for (long j = 0; j < n; ++j)
    {
      dd_set_d(M->matrix[i][j + 1], -A(i, j));
    }
  }

  dd_rowset implLin = nullptr;
  dd_rowset redset = nullptr;
  dd_rowindex newpos = nullptr;
  dd_ErrorType err = dd_NoError;
  dd_MatrixCanonicalize(&M, &implLin, &redset, &newpos, &err);

  bool keep = true;
  bool replace = false;
  Eigen::MatrixXd newA;
  Eigen::VectorXd newB;

  if (err != dd_NoError || M == nullptr)
  {
    // Leave the raw matrices untouched.
    keep = true;
  }
  else if (set_card(M->linset) > 0)
  {
    // Implicit equalities => measure-zero in the full dimension.
    keep = false;
  }
  else
  {
    long rows = M->rowsize;
    newA = Eigen::MatrixXd(rows, n);
    newB = Eigen::VectorXd(rows);
    for (long i = 0; i < rows; ++i)
    {
      newB(i) = dd_get_d(M->matrix[i][0]);
      for (long j = 0; j < n; ++j)
      {
        newA(i, j) = -dd_get_d(M->matrix[i][j + 1]);
      }
    }
    replace = true;
  }

  if (implLin != nullptr) set_free(implLin);
  if (redset != nullptr) set_free(redset);
  if (newpos != nullptr) free(newpos);
  if (M != nullptr) dd_FreeMatrix(M);
  dd_free_global_constants();

  if (replace)
  {
    A = std::move(newA);
    b = std::move(newB);
  }
  return keep;
}

// Stores the running set of sample points "X" used by the union-volume
// (Karp-Luby style) estimator.  Mirrors the behaviour of the reference Python
// implementation in src/cube_processor_nondis.py.  Templated on the coordinate
// type: `double` for the --nogmp path, `MpFloat` for the GMP paths (the
// membership test A x <= b is then evaluated in GMP, avoiding the double
// cancellation that can misclassify points near a facet).
template <typename Coord>
class SampleStore
{
 public:
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;
  using Mat = Eigen::Matrix<Coord, Eigen::Dynamic, Eigen::Dynamic>;

  explicit SampleStore(double tolerance) : d_tolerance(tolerance) {}

  void addSample(const Vec& sample)
  {
    d_samples.push_back(sample);
  }

  // Removes every stored point that lies inside the polytope {x : A x <= b}.
  // Equivalent to the Python:  X = [s for s in X if not poly.is_in(s)].  The
  // half-space matrices are double (built once); they are cast to the
  // coordinate type so the dot products run at the stored precision.
  std::size_t removeInside(const Eigen::MatrixXd& Ad, const Eigen::VectorXd& bd)
  {
    if (Ad.rows() == 0)
    {
      return 0;
    }
    Mat A = Ad.template cast<Coord>();
    Vec b = bd.template cast<Coord>();
    Coord tol(d_tolerance);
    auto oldSize = d_samples.size();
    auto it = std::remove_if(d_samples.begin(), d_samples.end(),
                             [&](const Vec& sample) {
                               Vec diff = A * sample - b;
                               return (diff.array() <= tol).all();
                             });
    d_samples.erase(it, d_samples.end());
    return oldSize - d_samples.size();
  }

  // Keeps each stored point independently with probability 1/2.  Mirrors the
  // Python:  X = [s for s in X if random.random() >= 0.5].
  template <typename URNG>
  std::size_t downsampleHalf(URNG& rng)
  {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto oldSize = d_samples.size();
    auto it = std::remove_if(d_samples.begin(), d_samples.end(),
                             [&](const Vec&) {
                               return unit(rng) < 0.5;
                             });
    d_samples.erase(it, d_samples.end());
    return oldSize - d_samples.size();
  }

  std::size_t size() const
  {
    return d_samples.size();
  }

 private:
  double d_tolerance;
  std::vector<Vec> d_samples;
};

// GetPrecision in the revised manuscript: eta = epsilon/96. Work in logs
// so gamma/minRadius need not be representable as a double.
int precisionFromCubes(std::size_t dimension, std::size_t facetCount,
                       std::size_t numCubes, double minRadius, double epsilon)
{
  const double eta = epsilon / 96.0;
  const double logGamma = std::log(64.0) + std::log(double(numCubes))
      + std::log(double(dimension)) - std::log(eta)
      + 0.5 * std::log(std::log(4.0) + std::log(double(facetCount)) - std::log(eta));
  const double prec = std::max(0.0, std::ceil((logGamma - std::log(minRadius))
                                           / std::log(10.0)));
  if (!std::isfinite(prec) || prec > std::numeric_limits<int>::max() - 32)
  {
    throw std::runtime_error("lattice precision exceeds supported range");
  }
  return static_cast<int>(prec);
}

// A polytope kept after phase 1, with its (double) half-space matrices and the
// volesti volume estimate.
struct PolytopeVolume
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  double volume = 0.0;
};

// Phase 2: streaming union-volume (Karp-Luby style) estimation, templated on the
// coordinate type used to store the surviving sample set `X`.  `p` is the
// running sampling probability; the union volume is |X| / (p*C), C=T*10^(b*n).
// `generate(A, b, n, samplingTime)` draws `n` fresh points inside {x : A x <= b}
// (the only step that differs between the double and GMP walk backends).
template <typename Coord, typename GenFn>
void runUnionAlgorithm(const std::vector<PolytopeVolume>& kept, double thresh,
                       unsigned seed, double log2Scale, GenFn&& generate,
                       VolumeComputationResult& result,
                       const std::function<void(const VolumeComputationRow&)>& onRow)
{
  std::mt19937 auxRng(seed);  // drives Poisson draws and down-sampling
  SampleStore<Coord> store(0.0);

  // Track log2(p*C), avoiding overflow of C and underflow of p.
  double log2Density = log2Scale;
  std::size_t rowIndex = 0;
  std::size_t totalGenerated = 0;
  std::size_t totalDeleted = 0;
  double totalSamplingTime = 0.0;
  double totalFilterTime = 0.0;

  for (const auto& entry : kept)
  {
    const double volume = entry.volume;

    // Drop stored points that lie inside the current polytope.
    const double filterStart = Log.elapsed();
    std::size_t deleted = store.removeInside(entry.A, entry.b);

    // Shrink the sampling probability until p*C*volume < thresh.
    while (log2Density + std::log2(volume) >= std::log2(thresh))
    {
      store.downsampleHalf(auxRng);
      log2Density -= 1.0;
    }

    // Draw the number of fresh samples; keep |X| + N within the threshold.
    std::poisson_distribution<long long> poisson(std::exp2(log2Density + std::log2(volume)));
    long long n = poisson(auxRng);
    while (static_cast<double>(n) + static_cast<double>(store.size()) > thresh)
    {
      store.downsampleHalf(auxRng);
      log2Density -= 1.0;
      n = std::binomial_distribution<long long>(n, 0.5)(auxRng);
    }

    totalFilterTime += Log.elapsed() - filterStart;

    std::size_t sampleCount = 0;
    double samplingTime = 0.0;
    if (n > 0)
    {
      auto samples = generate(entry.A, entry.b, n, samplingTime);
      const double storeStart = Log.elapsed();
      for (auto& sample : samples)
      {
        store.addSample(sample);
      }
      totalFilterTime += Log.elapsed() - storeStart;
      sampleCount = samples.size();
    }

    totalGenerated += sampleCount;
    totalDeleted += deleted;
    totalSamplingTime += samplingTime;
    VolumeComputationRow row;
    row.index = ++rowIndex;
    row.volume = volume;
    row.samplesDeleted = deleted;
    row.totalSamples = store.size();
    row.samplesGenerated = sampleCount;
    result.rows.push_back(row);
    if (onRow)
    {
      onRow(row);
    }
  }

  result.finalSampleCount = store.size();
  result.volumeEstimate =
      store.size() ? std::exp2(std::log2(double(store.size())) - log2Density) : 0.0;
  result.totalSamplesGenerated = totalGenerated;
  result.totalSamplesDeleted = totalDeleted;
  result.samplingTime = totalSamplingTime;
  result.filterTime = totalFilterTime;
}

std::vector<cvc5::Term> deduplicateRealVariables(
    const std::vector<cvc5::Term>& variables)
{
  std::vector<cvc5::Term> unique;
  unique.reserve(variables.size());
  std::unordered_set<cvc5::Term> seen;
  for (const auto& term : variables)
  {
    if (!term.getSort().isReal())
    {
      continue;
    }
    if (seen.insert(term).second)
    {
      unique.push_back(term);
    }
  }
  return unique;
}

}  // namespace

double lraSampleThreshold(double epsilon, double delta, std::size_t numPolytopes)
{
  const double mvcEps = epsilon / 12.0;
  const double m = std::max(1.0, static_cast<double>(numPolytopes));
  return std::max(
      24.0 * (std::log(24.0) - std::log(delta)) / ((1.0 - mvcEps) * mvcEps * mvcEps),
      15.0 * (std::log(6.0) - std::log(delta) + std::log(m)));
}

unsigned lraDefaultWalkLength(std::size_t dimension)
{
  return static_cast<unsigned>(10 + dimension * dimension);
}

VolumeComputationResult computeLraVolume(
    const std::vector<Polytope>& polytopes,
    const std::vector<cvc5::Term>& realVariables,
    const VolumeOptions& options,
    const std::function<void(const VolumeComputationRow&)>& onRow,
    const std::function<void(const VolumeSetup&)>& onSetup)
{
  VolumeComputationResult result;
  if (!std::isfinite(options.epsilon) || options.epsilon <= 0.0 || options.epsilon > 1.0
      || !std::isfinite(options.delta) || options.delta <= 0.0 || options.delta >= 1.0)
    throw std::invalid_argument("LRA volume requires 0 < epsilon <= 1 and 0 < delta < 1");
  if (options.sampleWalkLength == 0 || options.volumeWalkLength == 0 || options.precision < 0)
    throw std::invalid_argument("volume walk lengths and arithmetic precision must be positive or automatic");
  auto uniqueVars = deduplicateRealVariables(realVariables);
  if (uniqueVars.empty())
  {
    throw std::runtime_error("LRA volume computation requires at least one real variable");
  }

  std::size_t dimension = uniqueVars.size();
  TermIndexMap index;
  index.reserve(dimension);
  for (std::size_t i = 0; i < uniqueVars.size(); ++i)
  {
    index.emplace(uniqueVars[i], i);
  }

  // Algorithm 1 in the revised manuscript.
  const double kEpsilon = options.epsilon;
  const double kDelta = options.delta;
  const unsigned kSeed = options.seed;
  // Per-polytope sampling streams.  `generate` is invoked once per polytope by
  // runUnionAlgorithm, and each invocation seeds its walk from scratch, so a
  // single base seed would replay the identical direction sequence for every
  // polytope and correlate the sample sets the union estimator combines.
  // Derive a distinct stream per call instead, mixing with splitmix64 so
  // consecutive polytopes do not get near-adjacent mt19937 states.
  unsigned streamIndex = 0;
  auto streamSeed = [&kSeed, &streamIndex]() {
    std::uint64_t z = kSeed + 0x9E3779B97F4A7C15ULL * (++streamIndex);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return static_cast<unsigned>(z ^ (z >> 31));
  };
  const double mvcEps = kEpsilon / 12.0;
  const double numCubes = static_cast<double>(polytopes.size());
  const double oracleDelta = kDelta / (4.0 * std::max(1.0, numCubes));
  const double thresh = lraSampleThreshold(kEpsilon, kDelta, polytopes.size());
  if (!std::isfinite(thresh) || thresh >= double(std::numeric_limits<long long>::max() / 2)
      || oracleDelta <= 0.0)
    throw std::invalid_argument("epsilon/delta require an unsupported sample threshold");

  // The oracle walk defaults per oracle (defaultOracleWalkLength: 10 + dim^2
  // for the Gaussian ball walk) and the sampling walk to 10 + dim^2
  // (lraDefaultWalkLength); --walklen-vol / --walklen-samp override them (a
  // negative option value keeps the default).
  const unsigned volumeWalkLength =
      options.volumeWalkLength >= 0
          ? static_cast<unsigned>(options.volumeWalkLength)
          : defaultOracleWalkLength(options.volumeOracle, dimension);
  const unsigned sampleWalkLength =
      options.sampleWalkLength >= 0
          ? static_cast<unsigned>(options.sampleWalkLength)
          : defaultSampleWalkLength(options.samplerWalk, dimension);

  Log(2) << "Starting union algorithm, threshold: " << thresh << std::endl;
  // Volesti's practical cooling API accepts epsilon but no confidence delta.
  Log(2) << "Per-polytope oracle targets: epsilon=" << mvcEps
         << ", delta=" << oracleDelta << " (practical Volesti oracle)" << std::endl;

  // Phase 1a: build and canonicalize the half-space matrices and find every
  // inner ball and outer radius, dropping polytopes whose (full-dimensional)
  // volume is zero.  This fixes the lattice precision before the volume
  // oracles run (phase 1b), so the caller can report it up front.
  struct PreparedPolytope
  {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    InnerBall ball;  // certified Chebyshev ball
    double innerRadius;
    double outerRadius;
    long rawFacets;
  };
  std::vector<PreparedPolytope> prepared;
  prepared.reserve(polytopes.size());

  std::size_t numZeroVolume = 0;
  double totalVolumeTime = 0.0;
  double ballTime = 0.0;         // phase 1a: inner ball + outer-radius LPs
  double samplerBallTime = 0.0;  // phase 2: inner ball per sampler restart
  std::size_t facetTotal = 0;
  double minRadius = std::numeric_limits<double>::infinity();
  std::size_t dumpIndex = 0;  // running counter for --dump-ine file names
  for (const auto& poly : polytopes)
  {
    PolytopeMatrices matrices = buildPolytopeMatrices(poly, index, dimension);
    if (matrices.empty)
    {
      ++numZeroVolume;
      continue;
    }
    if (!matrices.A.allFinite() || !matrices.b.allFinite())
      throw std::runtime_error("polytope coefficients exceed double geometry range");
    // Canonicalize (remove redundant constraints, detect equalities) with
    // cddlib, exactly as the Python prototype does before volume/sampling.
    // --no-cdd-simp skips this and feeds volesti the raw half-spaces.
    long rawFacets = matrices.A.rows();
    if (options.cddSimplify && !canonicalizePolytope(matrices.A, matrices.b))
    {
      ++numZeroVolume;
      continue;
    }
    if (!options.dumpInePrefix.empty())
    {
      std::string path = options.dumpInePrefix + "_cube" +
                         std::to_string(++dumpIndex) + ".ine";
      writePolytopeIne(path, matrices.A, matrices.b);
      Log(1) << "wrote polytope to " << path << std::endl;
    }
    const double ballStart = Log.elapsed();
    InnerBall innerBall = certifiedInnerBall(matrices.A, matrices.b);
    ballTime += Log.elapsed() - ballStart;
    if (innerBall.radius <= 0.0)
    {
      ++numZeroVolume;
      continue;
    }
    if (!std::isfinite(innerBall.radius))
      throw std::runtime_error("polytope has no finite inner radius");
    const double radiusStart = Log.elapsed();
    const double radius = outerRadius(matrices.A, matrices.b, innerBall.center);
    ballTime += Log.elapsed() - radiusStart;
    minRadius = std::min(minRadius, innerBall.radius);
    facetTotal += static_cast<std::size_t>(matrices.A.rows());
    const double innerRadius = innerBall.radius;
    prepared.push_back({std::move(matrices.A), std::move(matrices.b), std::move(innerBall),
                        innerRadius, radius, rawFacets});
  }

  if (numZeroVolume > 0)
  {
    Log(2) << "Skipping " << numZeroVolume << " polytopes out of "
           << polytopes.size() << " where volume is zero" << std::endl;
  }

  // Lattice spacing and GMP arithmetic precision are distinct quantities.
  const int precision =
      prepared.empty() ? 0
                       : precisionFromCubes(dimension, facetTotal, polytopes.size(),
                                            minRadius, kEpsilon);
  result.latticePrecision = precision;
  if (options.gmpMode != GmpMode::None && !prepared.empty())
  {
    result.samplingPrecision = std::max(options.precision, precision + 32);
    MpFloat::default_precision(static_cast<unsigned>(result.samplingPrecision));
  }
  if (onSetup)
  {
    VolumeSetup setup;
    setup.dimension = dimension;
    setup.keptPolytopes = prepared.size();
    setup.volumeWalkLength = volumeWalkLength;
    setup.sampleWalkLength = sampleWalkLength;
    setup.latticePrecision = result.latticePrecision;
    setup.samplingPrecision = result.samplingPrecision;
    onSetup(setup);
  }

  // Phase 1b: run the volume oracle on every polytope kept by phase 1a.
  std::vector<PolytopeVolume> kept;
  kept.reserve(prepared.size());
  for (auto& entry : prepared)
  {
    double volumeComputeStart = Log.elapsed();
    // Round first (Cousins-Vempala, Fig. 2 step 1): every volesti walk steps in
    // units of the inner radius, so on a skinny body it never crosses the body
    // and the estimate collapses.  vol(P) = det * vol(rounded P), folded into
    // RoundedPolytope::volume.  The original body keeps its radii for clipping,
    // GetPrecision and sampling.
    const RoundedPolytope body(entry.A, entry.b, entry.ball, options.roundPolytopes);
    if (body.rounded())
      Log(2) << "rounded polytope: |det T| = " << body.det() << ", inner radius "
             << entry.innerRadius << " -> " << body.innerRadius() << std::endl;
    else if (options.roundPolytopes)
      Log(1) << "WARNING: polytope rounding skipped (" << body.roundingFailure()
             << "); the volume estimate may be poor if the polytope is skinny" << std::endl;
    double volume = body.volume(options.volumeOracle, mvcEps, volumeWalkLength);
    // A practical oracle run can still end non-finite or non-positive, which is
    // a detected failure.  Rerun it on fresh independent streams before giving up.
    constexpr unsigned kOracleAttempts = 3;
    for (unsigned attempt = 1;
         attempt < kOracleAttempts && !(std::isfinite(volume) && volume > 0.0);
         ++attempt)
    {
      Log(1) << "WARNING: volume oracle returned " << volume
             << " (inner radius " << entry.innerRadius << ", outer radius "
             << entry.outerRadius << "); retrying on a fresh random stream" << std::endl;
      volume = body.volumeWithSeed(options.volumeOracle, mvcEps, volumeWalkLength,
                                   kSeed + 0x51ED2701U * attempt);
    }
    double volumeComputeElapsed = Log.elapsed() - volumeComputeStart;
    totalVolumeTime += volumeComputeElapsed;
    Log(2) << "polytope volume " << volume << " (" << entry.A.rows()
           << " facets, was " << entry.rawFacets << ", dim " << dimension << ") in "
           << volumeComputeElapsed << "s" << std::endl;
    if (!std::isfinite(volume) || volume <= 0.0)
      throw std::runtime_error("volume oracle failed on a positive-dimensional polytope");
    const double lower = (1.0 - mvcEps) * std::pow(2.0 * entry.innerRadius / std::sqrt(double(dimension)), double(dimension));
    const double upper = (1.0 + mvcEps) * std::pow(2.0 * entry.outerRadius, double(dimension));
    if (std::isnan(lower) || std::isnan(upper) || lower > upper)
      throw std::runtime_error("invalid polytope volume clipping bounds");
    volume = std::clamp(volume, lower, upper);
    if (!std::isfinite(volume) || volume <= 0.0)
      throw std::runtime_error("clipped polytope volume exceeds double range");
    kept.push_back({std::move(entry.A), std::move(entry.b), volume});
  }
  prepared.clear();

  result.volumeComputationTime = totalVolumeTime;
  result.ballTime = ballTime;
  if (kept.empty()) return result;
  const double log2Scale = std::log2(thresh) + double(precision) * double(dimension) * std::log2(10.0);

  // Restart the selected walk at the certified centre before every proposal.
  auto generateDouble =
      [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
          double& samplingTime) -> std::vector<Eigen::VectorXd> {
    double ballSeconds = 0.0;
    double samplingStart = Log.elapsed();
    auto out = sampleRestartWalk(A, b, n, options.samplerWalk, sampleWalkLength, streamSeed(),
                                 ballSeconds, options.roundPolytopes);
    samplerBallTime += ballSeconds;
    samplingTime = Log.elapsed() - samplingStart - ballSeconds;
    return out;
  };

  if (options.gmpMode == GmpMode::Full &&
      options.samplerWalk == SamplerWalk::AcceleratedBilliard)
  {
    Log(0) << "WARNING: --fullgmp has no accelerated-billiard implementation; "
              "using the GMP billiard walk instead"
           << std::endl;
  }

  switch (options.gmpMode)
  {
    case GmpMode::None:
    {
      auto generate = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                          long long n, double& time) {
        return sampleLattice<double>(A, b, n, precision, mvcEps, oracleDelta,
                                     streamSeed(), generateDouble, time);
      };
      runUnionAlgorithm<double>(kept, thresh, kSeed, log2Scale, generate, result,
                                onRow);
      break;
    }
    case GmpMode::PointRepr:
    {
      // Double proposals, with lattice rounding and membership in GMP.
      auto generate =
          [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
              double& samplingTime) {
            return sampleLattice<MpFloat>(A, b, n, precision, mvcEps, oracleDelta,
                                          streamSeed(), generateDouble, samplingTime);
          };
      runUnionAlgorithm<MpFloat>(kept, thresh, kSeed, log2Scale, generate, result, onRow);
      break;
    }
    case GmpMode::Full:
    {
      // Run the walk itself in GMP.  The Chebyshev centre and inner radius
      // (used to size the walk step) come from the double inner ball.  Billiard,
      // ball, RDHR and CDHR have GMP implementations; accelerated-billiard was
      // warned about above and falls back to the GMP billiard walk.
      const SamplerWalk walk = options.samplerWalk;
      auto generate =
          [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
              double& samplingTime) {
            const double ballStart = Log.elapsed();
            const InnerBall ball = certifiedInnerBall(A, b);
            // As in sampleRestartWalk: walk in the rounded body and map back.
            PolytopeRounding rounding;
            if (options.roundPolytopes) rounding = roundPolytope(A, b, ball);
            const Eigen::MatrixXd& Aw = rounding.applied ? rounding.A : A;
            const Eigen::VectorXd& bw = rounding.applied ? rounding.b : b;
            samplerBallTime += Log.elapsed() - ballStart;
            const Eigen::VectorXd& center = rounding.applied ? rounding.ball.center : ball.center;
            const double radius = rounding.applied ? rounding.ball.radius : ball.radius;
            double samplingStart = Log.elapsed();
            // One stream advance per polytope, drawn before the switch so the
            // sequence of seeds does not depend on which walk is selected.
            const unsigned walkSeed = streamSeed();
            std::vector<MpVector> pts;
            switch (walk)
            {
              case SamplerWalk::Ball:
                pts = sampleGmpBallWalk(Aw, bw, center, radius, n,
                                        sampleWalkLength, walkSeed, true);
                break;
              case SamplerWalk::RDHR:
                pts = sampleGmpRDHR(Aw, bw, center, radius, n, sampleWalkLength,
                                    walkSeed, true);
                break;
              case SamplerWalk::CDHR:
                pts = sampleGmpCDHR(Aw, bw, center, radius, n, sampleWalkLength,
                                    walkSeed, true);
                break;
              case SamplerWalk::AcceleratedBilliard:
              case SamplerWalk::Billiard:
              default:
                pts = sampleGmpBilliard(Aw, bw, center, radius, n,
                                        sampleWalkLength, walkSeed, true);
                break;
            }
            if (rounding.applied)
            {
              const MpVector shift = rounding.shift.cast<MpFloat>();
              const auto T = rounding.T.cast<MpFloat>().eval();
              for (auto& y : pts) y = shift + T * y;
            }
            samplingTime = Log.elapsed() - samplingStart;
            return pts;
          };
      auto latticeGenerate = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                 long long n, double& time) {
        return sampleLattice<MpFloat>(A, b, n, precision, mvcEps, oracleDelta,
                                      streamSeed(), generate, time);
      };
      runUnionAlgorithm<MpFloat>(kept, thresh, kSeed, log2Scale, latticeGenerate, result, onRow);
      break;
    }
  }

  // sampleLattice times its whole call, including the sampler's inner-ball
  // computation; report that part under ball time only.
  result.samplingTime = std::max(0.0, result.samplingTime - samplerBallTime);
  result.ballTime = ballTime + samplerBallTime;
  result.volumeComputationTime = totalVolumeTime;
  return result;
}

}  // namespace ttc

