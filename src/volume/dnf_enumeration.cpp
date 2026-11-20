#include "volume/dnf_enumeration.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "volume/aig_allsat.hpp"

namespace ttc
{
namespace
{

using AtomMap = std::unordered_map<std::uint32_t, cvc5::Term>;  // M: input variable -> atom

// The cube over theory atoms that `input` stands for, or nullopt when two
// inputs that M maps to the same atom disagree.  The literals come in the
// iteration order of the per-atom map below, and that order becomes the row
// order of A: changing the container gives different (equally valid)
// estimates on the same seed.
std::optional<Cube> theoryCube(const InputCube& input, const AtomMap& M)
{
  std::unordered_map<cvc5::Term, std::optional<bool>> valueOf;
  valueOf.reserve(input.size());
  for (const auto& [var, value] : input)
  {
    const auto atom = M.find(var);
    if (atom == M.end()) continue;
    const auto [it, inserted] = valueOf.emplace(atom->second, value);
    if (inserted || !value) continue;
    if (!it->second)
      it->second = value;
    else if (*it->second != *value)
      return std::nullopt;
  }
  Cube cube;
  for (const auto& [atom, value] : valueOf)
    if (value) cube.literals.push_back({atom, *value});
  return cube;
}

// The literals of a cube as sorted "atom=1" / "atom=0" strings: equal keys
// are equal cubes, and key inclusion is cube inclusion.
std::vector<std::string> literalKey(const Cube& cube)
{
  std::vector<std::string> key;
  key.reserve(cube.literals.size());
  for (const TermLiteral& literal : cube.literals)
    key.push_back(literal.term.toString() + (literal.positive ? "=1" : "=0"));
  std::sort(key.begin(), key.end());
  return key;
}

bool strictlyContains(const std::vector<std::string>& larger,
                      const std::vector<std::string>& smaller)
{
  return smaller.size() < larger.size()
         && std::includes(larger.begin(), larger.end(), smaller.begin(), smaller.end());
}

}  // namespace

std::vector<Cube> enumerateCubes(const BooleanAbstractionAigData& aig)
{
  AtomMap M;
  M.reserve(aig.mapping.size());
  for (const auto& [var, atom] : aig.mapping) M.emplace(var, atom);

  std::vector<Cube> cubes;
  std::vector<std::vector<std::string>> keys;
  std::set<std::vector<std::string>> seen;
  for (const InputCube& input : allSatInputCubes(aig))
  {
    std::optional<Cube> cube = theoryCube(input, M);
    if (!cube) continue;
    std::vector<std::string> key = literalKey(*cube);
    if (!seen.insert(key).second) continue;
    cubes.push_back(std::move(*cube));
    keys.push_back(std::move(key));
  }

  std::vector<Cube> minimal;
  for (std::size_t i = 0; i < cubes.size(); ++i)
  {
    const bool subsumed = std::any_of(keys.begin(), keys.end(), [&](const auto& other) {
      return strictlyContains(keys[i], other);
    });
    if (!subsumed) minimal.push_back(std::move(cubes[i]));
  }
  return minimal;
}

}  // namespace ttc
