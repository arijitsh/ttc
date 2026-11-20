#include "volume/unionvolume.hpp"

// Algorithm 1, ttc(F, eps, delta), and the file that implements each line:
//
//    1  F_bc, M <- BoolAbs(F)                         boolean_abstraction.hpp \ run by
//    2  F_bd <- ToDNF(F_bc)                           dnf_enumeration.hpp     / TTCParser
//       P_i = P(c_i) for each cube c_i of F_bd         prepared_polytope.hpp
//    3  b <- GetPrecision(F_bd, M, eps/96)            get_precision.hpp
//    4  eps' <- eps/12,  delta' <- delta/(4m)          union_parameters.hpp
//    5  Thresh <- max(...)                             union_parameters.hpp
//    6  p <- 1;  X <- {}                               union_estimator.hpp
//    7  for i = 1 to m
//    8    t <- Thresh 10^{bn} VolEst(P_i, eps', delta')          vol_est.hpp
//    9-18 thin X and draw N_i                         union_estimator.hpp
//   19    S <- GenerateSamples(P_i, N_i, b, eps', delta')        generate_samples.hpp
//   20    X.Append(S)                                 union_estimator.hpp
//   22  output |X| / (p Thresh 10^{bn})               union_estimator.hpp
//
// m counts every cube, but only polytopes of positive volume enter the loop.

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "logger.hpp"
#include "volume/generate_samples.hpp"
#include "volume/get_precision.hpp"
#include "volume/linear_constraints.hpp"
#include "volume/mp_float.hpp"
#include "volume/polytope_sampler.hpp"
#include "volume/prepared_polytope.hpp"
#include "volume/seed_streams.hpp"
#include "volume/union_estimator.hpp"
#include "volume/union_parameters.hpp"
#include "volume/vol_est.hpp"

namespace ttc
{
namespace
{

void checkOptions(const VolumeOptions& options)
{
  if (!std::isfinite(options.epsilon) || options.epsilon <= 0.0 || options.epsilon > 1.0
      || !std::isfinite(options.delta) || options.delta <= 0.0 || options.delta >= 1.0)
    throw std::invalid_argument("LRA volume requires 0 < epsilon <= 1 and 0 < delta < 1");
  if (options.sampleWalkLength == 0 || options.volumeWalkLength == 0 || options.precision < 0)
    throw std::invalid_argument(
        "volume walk lengths and arithmetic precision must be positive or automatic");
}

}  // namespace

VolumeComputationResult computeLraVolume(const std::vector<Cube>& cubes,
                                         const std::vector<cvc5::Term>& realVariables,
                                         const VolumeOptions& options, const RowCallback& onRow,
                                         const SetupCallback& onSetup)
{
  checkOptions(options);
  const RealVariables variables(realVariables);
  if (variables.size() == 0)
    throw std::runtime_error("LRA volume computation requires at least one real variable");
  const std::size_t n = variables.size();
  const std::size_t m = cubes.size();

  // Lines 1-2 ran in the parser, so `cubes` is F_bd.  P_i for each cube,
  // without the polytopes of zero volume.
  const PreparedPolytopes prepared =
      preparePolytopes(cubes, variables, {options.cddSimplify, options.dumpInePrefix});
  const std::vector<PreparedPolytope>& P = prepared.kept;

  const int b = getPrecision(P, n, m, options.epsilon / 96.0);                               // line 3
  const UnionParameters params = unionParameters(options.epsilon, options.delta, m, n, b);  // lines 4-5
  Log(2) << "Starting union algorithm, threshold: " << params.thresh << std::endl;
  Log(2) << "Per-polytope oracle targets: epsilon=" << params.epsPrime
         << ", delta=" << params.deltaPrime << " (practical Volesti oracle)" << std::endl;

  // Arithmetic of X and GenerateSamples, and the walk lengths.
  const double digits = latticeDigits(b, maxAbsCoordinate(P));
  const GmpMode gmpMode = resolveGmpMode(options.gmpMode, digits);
  if (gmpMode != GmpMode::Full && digits > 15.9)
    Log(1) << "WARNING: lattice coordinates need " << digits
           << " digits; the double walk cannot resolve the 10^-" << b
           << " lattice (use --fullgmp for a faithful lattice sampler)" << std::endl;
  VolumeComputationResult result;
  result.latticePrecision = b;
  result.ballTime = prepared.ballSeconds;
  if (gmpMode != GmpMode::None && !P.empty())
  {
    result.samplingPrecision = gmpDigits(options.precision, b);
    MpFloat::default_precision(static_cast<unsigned>(result.samplingPrecision));
  }
  const unsigned volumeWalkLength = options.volumeWalkLength >= 0
                                        ? static_cast<unsigned>(options.volumeWalkLength)
                                        : defaultOracleWalkLength(options.volumeOracle, n);
  const unsigned sampleWalkLength = options.sampleWalkLength >= 0
                                        ? static_cast<unsigned>(options.sampleWalkLength)
                                        : defaultSampleWalkLength(options.samplerWalk, n);
  if (onSetup)
    onSetup({n, P.size(), volumeWalkLength, sampleWalkLength, b, result.samplingPrecision,
             gmpMode, digits});
  if (P.empty()) return result;
  if (gmpMode == GmpMode::Full && options.samplerWalk == SamplerWalk::AcceleratedBilliard)
    Log(0) << "WARNING: --fullgmp has no accelerated-billiard implementation; "
              "using the GMP billiard walk instead"
           << std::endl;

  // Lines 6-22, with VolEst (line 8) and GenerateSamples (line 19) as oracles.
  VolEst volEst(options.volumeOracle, volumeWalkLength, options.roundPolytopes, params.epsPrime,
                options.seed);
  auto volEstOf = [&](std::size_t i) { return volEst(P[i], i); };
  // Seeds of the sampling side, in call order: one for each GenerateSamples
  // call (its randomized rounding), then one per batch of walk proposals.
  SeedStream seeds(options.seed);
  double walkBallSeconds = 0.0;  // inner balls of the walks' bodies
  // Runs the loop with X and GenerateSamples in `Coord` arithmetic and the
  // proposals of `restartWalk` (sampleRestartWalk or sampleRestartWalkGmp).
  auto runLoop = [&](auto coordTag, auto restartWalk) {
    using Coord = decltype(coordTag);
    auto proposals = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& rhs, long long count,
                         double&) {
      double ballSeconds = 0.0;
      auto points = restartWalk(A, rhs, count, options.samplerWalk, sampleWalkLength, seeds.next(),
                                ballSeconds, options.roundPolytopes);
      walkBallSeconds += ballSeconds;
      return points;
    };
    auto generate = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& rhs, long long N,
                        double& seconds) {
      return generateSamples<Coord>(A, rhs, N, b, params.epsPrime, params.deltaPrime, seeds.next(),
                                    proposals, seconds);
    };
    estimateUnionVolume<Coord>(P, params.thresh, params.log2c, options.seed, volEstOf, generate,
                               result, onRow);
  };
  switch (gmpMode)
  {
    case GmpMode::None: runLoop(double{}, sampleRestartWalk); break;
    case GmpMode::PointRepr: runLoop(MpFloat{}, sampleRestartWalk); break;  // GMP X, double walk
    case GmpMode::Full: runLoop(MpFloat{}, sampleRestartWalkGmp); break;    // the walk in GMP too
    case GmpMode::Auto: throw std::logic_error("GmpMode::Auto must be resolved before sampling");
  }

  // GenerateSamples times its whole call; report the walks' inner balls under ball time only.
  result.samplingTime = std::max(0.0, result.samplingTime - walkBallSeconds);
  result.ballTime += walkBallSeconds;
  result.volumeComputationTime = volEst.seconds();
  return result;
}

}  // namespace ttc
