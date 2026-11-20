#include "volume/prepared_polytope.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "logger.hpp"
#include "volume/canonicalize.hpp"
#include "volume/cube.hpp"
#include "volume/linear_constraints.hpp"
#include "volume/ine_writer.hpp"

namespace ttc
{

PreparedPolytopes preparePolytopes(const std::vector<Cube>& cubes, const RealVariables& variables,
                                   const PrepareOptions& options)
{
  PreparedPolytopes out;
  std::size_t dumped = 0;
  auto timed = [&out](auto&& lp) {
    const double start = Log.elapsed();
    auto value = lp();
    out.ballSeconds += Log.elapsed() - start;
    return value;
  };

  for (const Cube& cube : cubes)
  {
    LinearConstraints P = cubeConstraints(cube, variables);
    if (P.empty)
    {
      ++out.zeroVolume;
      continue;
    }
    if (!P.A.allFinite() || !P.b.allFinite())
      throw std::runtime_error("polytope coefficients exceed double geometry range");
    const long rawFacets = P.A.rows();
    if (options.cddSimplify && !canonicalizePolytope(P.A, P.b))
    {
      ++out.zeroVolume;
      continue;
    }
    if (!options.dumpInePrefix.empty())
    {
      const std::string path = options.dumpInePrefix + "_cube" + std::to_string(++dumped) + ".ine";
      writePolytopeIne(path, P.A, P.b, variables.names());
      Log(1) << "wrote polytope to " << path << std::endl;
    }

    InnerBall ball = timed([&] { return certifiedInnerBall(P.A, P.b); });
    if (ball.radius <= 0.0)
    {
      ++out.zeroVolume;
      continue;
    }
    if (!std::isfinite(ball.radius)) throw std::runtime_error("polytope has no finite inner radius");
    const double R = timed([&] { return outerRadius(P.A, P.b, ball.center); });
    out.kept.push_back({std::move(P.A), std::move(P.b), std::move(ball), R, rawFacets});
  }

  if (out.zeroVolume > 0)
    Log(2) << "Skipping " << out.zeroVolume << " polytopes out of " << cubes.size()
           << " where volume is zero" << std::endl;
  return out;
}

double maxAbsCoordinate(const std::vector<PreparedPolytope>& P)
{
  double bound = 0.0;
  for (const PreparedPolytope& Pi : P)
    bound = std::max(bound, Pi.ball.center.cwiseAbs().maxCoeff() + Pi.outerRadius);
  return bound;
}

}  // namespace ttc
