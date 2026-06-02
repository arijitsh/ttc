#include "volume/lra_volume.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

#include <random>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
#include "logger.hpp"
#include "sampling/sphere.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/sampling_policies.hpp"
#include "volume/volume_cooling_gaussians.hpp"

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

constexpr double kZeroTolerance = 1e-9;

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
// implementation in src/cube_processor_nondis.py.
class SampleStore
{
 public:
  explicit SampleStore(double tolerance) : d_tolerance(tolerance) {}

  void addSample(const Eigen::VectorXd& sample)
  {
    d_samples.push_back(sample);
  }

  // Removes every stored point that lies inside the polytope {x : A x <= b}.
  // Equivalent to the Python:  X = [s for s in X if not poly.is_in(s)].
  std::size_t removeInside(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
  {
    if (A.rows() == 0)
    {
      return 0;
    }
    auto oldSize = d_samples.size();
    auto it = std::remove_if(d_samples.begin(), d_samples.end(),
                             [&](const Eigen::VectorXd& sample) {
                               Eigen::VectorXd diff = A * sample - b;
                               return (diff.array() <= d_tolerance).all();
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
                             [&](const Eigen::VectorXd&) {
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
  std::vector<Eigen::VectorXd> d_samples;
};

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

VolumeComputationResult computeLraVolume(
    const std::vector<Polytope>& polytopes,
    const std::vector<cvc5::Term>& realVariables,
    const std::function<void(const VolumeComputationRow&)>& onRow)
{
  VolumeComputationResult result;
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

  using Kernel = Cartesian<double>;
  using Point = Kernel::Point;
  using HPolytope = HPolytope<Point>;
  // Volume RNG matches the volesti `volume` tool: a fixed compile-time seed (3),
  // so the cooling-gaussians estimate is deterministic.
  using VolumeRNG = BoostRandomNumberGenerator<boost::mt19937, double, 3>;
  // Sampling RNG matches the volesti `sample` tool: explicitly seeded.
  using SampleRNG = BoostRandomNumberGenerator<boost::mt19937, double>;
  using SampleWalk = BilliardWalk::template Walk<HPolytope, SampleRNG>;

  // Algorithm parameters; defaults taken from the reference Python prototype
  // (src/cube_processor_nondis.py + src/global_storage.py).
  constexpr double kEpsilon = 0.8;
  constexpr double kDelta = 0.2;
  constexpr unsigned kSeed = 123;
  const double mvcEps = kEpsilon / 2.0;
  const double numCubes = static_cast<double>(polytopes.size());
  const double thresh = std::max(
      12.0 * std::log(24.0 / kDelta) / (mvcEps * mvcEps),
      6.0 * (std::log(6.0 / kDelta) + std::log(std::max(1.0, numCubes))));

  // The volesti `volume` tool uses walk length 10 + dim/10; the `sample` tool
  // uses a fixed walk length of 10.
  const unsigned volumeWalkLength = static_cast<unsigned>(10 + dimension / 10);
  const unsigned sampleWalkLength = 10;

  Log(2) << "Starting union algorithm, threshold: " << thresh << std::endl;

  // Phase 1: build the half-space matrices and compute every volume up front,
  // dropping polytopes whose (full-dimensional) volume is zero.
  struct PolytopeVolume
  {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    double volume = 0.0;
  };
  std::vector<PolytopeVolume> kept;
  kept.reserve(polytopes.size());

  std::size_t numZeroVolume = 0;
  double totalVolumeTime = 0.0;
  for (const auto& poly : polytopes)
  {
    PolytopeMatrices matrices = buildPolytopeMatrices(poly, index, dimension);
    if (matrices.empty)
    {
      ++numZeroVolume;
      continue;
    }
    // Canonicalize (remove redundant constraints, detect equalities) with
    // cddlib, exactly as the Python prototype does before volume/sampling.
    long rawFacets = matrices.A.rows();
    if (!canonicalizePolytope(matrices.A, matrices.b))
    {
      ++numZeroVolume;
      continue;
    }
    HPolytope hpoly(static_cast<unsigned>(dimension), matrices.A, matrices.b);
    if (hpoly.ComputeInnerBall().second <= 0.0)
    {
      ++numZeroVolume;
      continue;
    }
    double volumeComputeStart = Log.elapsed();
    double volume = volume_cooling_gaussians<GaussianBallWalk, VolumeRNG>(
        hpoly, 0.1, volumeWalkLength);
    double volumeComputeElapsed = Log.elapsed() - volumeComputeStart;
    totalVolumeTime += volumeComputeElapsed;
    Log(2) << "polytope volume " << volume << " (" << matrices.A.rows()
           << " facets, was " << rawFacets << ", dim " << dimension << ") in "
           << volumeComputeElapsed << "s" << std::endl;
    if (!std::isfinite(volume) || volume <= 0.0)
    {
      ++numZeroVolume;
      continue;
    }
    kept.push_back({std::move(matrices.A), std::move(matrices.b), volume});
  }

  if (numZeroVolume > 0)
  {
    Log(2) << "Skipping " << numZeroVolume << " polytopes out of "
           << polytopes.size() << " where volume is zero" << std::endl;
  }

  // Phase 2: streaming union-volume (Karp-Luby style) estimation.  `p` is the
  // running sampling probability and `store` holds the surviving samples `X`;
  // the union volume is estimated as |X| / p.
  std::mt19937 auxRng(kSeed);  // drives Poisson draws and down-sampling
  SampleRNG sampleRng(static_cast<int>(dimension));
  SampleStore store(1e-6);

  double p = 1.0;
  std::size_t rowIndex = 0;
  std::size_t totalGenerated = 0;
  std::size_t totalDeleted = 0;
  double totalSamplingTime = 0.0;

  for (const auto& entry : kept)
  {
    const double volume = entry.volume;

    // Drop stored points that lie inside the current polytope.
    std::size_t deleted = store.removeInside(entry.A, entry.b);

    // Shrink the sampling probability so that p*volume stays within thresh.
    while (p * volume > thresh)
    {
      store.downsampleHalf(auxRng);
      p /= 2.0;
    }

    // Draw the number of fresh samples; keep |X| + N within the threshold.
    std::poisson_distribution<long long> poisson(p * volume);
    long long n = poisson(auxRng);
    while (static_cast<double>(n) + static_cast<double>(store.size()) > thresh)
    {
      store.downsampleHalf(auxRng);
      p /= 2.0;
      poisson = std::poisson_distribution<long long>(p * volume);
      n = poisson(auxRng);
    }

    std::size_t sampleCount = 0;
    double samplingTime = 0.0;
    if (n > 0)
    {
      HPolytope hpoly(static_cast<unsigned>(dimension), entry.A, entry.b);
      Point start = hpoly.ComputeInnerBall().first;
      PushBackWalkPolicy policy;
      std::list<Point> samples;
      sampleRng.set_seed(kSeed);
      double samplingStart = Log.elapsed();
      RandomPointGenerator<SampleWalk>::apply(hpoly, start,
                                              static_cast<unsigned>(n),
                                              sampleWalkLength, samples, policy,
                                              sampleRng);
      samplingTime = Log.elapsed() - samplingStart;
      for (const Point& sample : samples)
      {
        store.addSample(sample.getCoefficients());
      }
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
      (p > 0.0) ? static_cast<double>(store.size()) / p : 0.0;
  result.totalSamplesGenerated = totalGenerated;
  result.totalSamplesDeleted = totalDeleted;
  result.volumeComputationTime = totalVolumeTime;
  result.samplingTime = totalSamplingTime;
  return result;
}

}  // namespace ttc

