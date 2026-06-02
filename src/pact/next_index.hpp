#pragma once

#include <cstddef>
#include <optional>

// Progression strategy for the number of hash constraints. The strategy
// performs an exponential search followed by binary search to identify the
// smallest number of hashes that yields an exact model count. It also keeps
// track of information between rounds so that the search in round n can start
// from the satisfying index discovered in round n-1.
class NextIndex
{
 public:
  NextIndex();

  std::size_t initial() const;

  void startRound();

  std::size_t nextCandidate();

  std::optional<std::size_t> previewNextCandidate() const;

  void updateOnSaturation(std::size_t level);
  void updateOnSatisfying(std::size_t level);

  void finishRound();

  bool needsRefinement() const;

  bool hasSatisfying() const;
  std::size_t currentSatisfying() const;
  std::size_t currentLower() const;

 private:
  std::size_t d_lastSatisfying;
  bool d_hasLastSatisfying;
  std::size_t d_lastSaturating;
  bool d_hasLastSaturating;

  bool d_roundActive;
  bool d_firstCandidate;
  bool d_haveLow;
  std::size_t d_low;
  bool d_lowValidated;
  bool d_haveHigh;
  std::size_t d_high;
  std::size_t d_step;
};

