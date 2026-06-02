#include "pact/next_index.hpp"

#include <algorithm>
#include <limits>

NextIndex::NextIndex()
    : d_lastSatisfying(0),
      d_hasLastSatisfying(false),
      d_lastSaturating(0),
      d_hasLastSaturating(false),
      d_roundActive(false),
      d_firstCandidate(true),
      d_haveLow(false),
      d_low(0),
      d_lowValidated(true),
      d_haveHigh(false),
      d_high(0),
      d_step(1)
{}

std::size_t NextIndex::initial() const { return 0; }

void NextIndex::startRound()
{
  d_roundActive = true;
  d_firstCandidate = true;
  if (d_hasLastSaturating)
  {
    d_haveLow = true;
    d_low = d_lastSaturating;
    d_lowValidated = (d_low == 0);
  }
  else
  {
    d_haveLow = false;
    d_low = 0;
    d_lowValidated = true;
  }
  if (d_hasLastSatisfying)
  {
    d_haveHigh = true;
    d_high = d_lastSatisfying;
  }
  else
  {
    d_haveHigh = false;
    d_high = 0;
  }
  if (d_haveHigh && d_haveLow && d_high <= d_low)
  {
    d_haveHigh = false;
  }
  d_step = 1;
}

std::size_t NextIndex::nextCandidate()
{
  if (!d_roundActive)
  {
    startRound();
  }
  if (d_firstCandidate)
  {
    d_firstCandidate = false;
    if (d_haveHigh)
    {
      return d_high;
    }
  }

  if (d_haveHigh && d_haveLow && !d_lowValidated)
  {
    return d_low;
  }

  if (d_haveHigh && d_haveLow)
  {
    if (d_low + 1 >= d_high)
    {
      return d_high;
    }
    std::size_t mid = d_low + (d_high - d_low) / 2;
    if (mid <= d_low)
    {
      mid = d_low + 1;
    }
    if (mid >= d_high)
    {
      mid = d_high - 1;
    }
    return mid;
  }

  std::size_t base = d_haveLow ? d_low : 0;
  std::size_t candidate = base + d_step;
  if (candidate <= base)
  {
    candidate = base + 1;
  }
  if (candidate < base)  // overflow guard
  {
    candidate = std::numeric_limits<std::size_t>::max();
  }
  if (d_step > std::numeric_limits<std::size_t>::max() / 2)
  {
    d_step = std::numeric_limits<std::size_t>::max();
  }
  else
  {
    d_step *= 2;
    if (d_step == 0)
    {
      d_step = 1;
    }
  }
  return candidate;
}

std::optional<std::size_t> NextIndex::previewNextCandidate() const
{
  bool roundActive = d_roundActive;
  bool firstCandidate = d_firstCandidate;
  bool haveLow = d_haveLow;
  std::size_t low = d_low;
  bool lowValidated = d_lowValidated;
  bool haveHigh = d_haveHigh;
  std::size_t high = d_high;
  std::size_t step = d_step;

  if (!roundActive)
  {
    roundActive = true;
    firstCandidate = true;
    if (d_hasLastSaturating)
    {
      haveLow = true;
      low = d_lastSaturating;
      lowValidated = (low == 0);
    }
    else
    {
      haveLow = false;
      low = 0;
      lowValidated = true;
    }
    if (d_hasLastSatisfying)
    {
      haveHigh = true;
      high = d_lastSatisfying;
    }
    else
    {
      haveHigh = false;
      high = 0;
    }
    if (haveHigh && haveLow && high <= low)
    {
      haveHigh = false;
    }
    step = 1;
  }

  if (firstCandidate)
  {
    if (haveHigh)
    {
      return high;
    }
    firstCandidate = false;
  }

  if (haveHigh && haveLow && !lowValidated)
  {
    return low;
  }

  if (haveHigh && haveLow)
  {
    if (low + 1 >= high)
    {
      return high;
    }
    std::size_t mid = low + (high - low) / 2;
    if (mid <= low)
    {
      mid = low + 1;
    }
    if (mid >= high)
    {
      mid = high > 0 ? high - 1 : 0;
    }
    return mid;
  }

  std::size_t base = haveLow ? low : 0;
  std::size_t candidate = base + step;
  if (candidate <= base)
  {
    candidate = base + 1;
  }
  if (candidate < base)
  {
    candidate = std::numeric_limits<std::size_t>::max();
  }
  return candidate;
}

void NextIndex::updateOnSaturation(std::size_t level)
{
  if (!d_hasLastSaturating || level > d_lastSaturating)
  {
    d_lastSaturating = level;
  }
  d_hasLastSaturating = true;

  if (!d_roundActive)
  {
    return;
  }

  if (!d_haveLow || level > d_low)
  {
    d_low = level;
    d_haveLow = true;
  }
  if (level == d_low)
  {
    d_lowValidated = true;
  }
  if (d_haveHigh && level >= d_high)
  {
    d_haveHigh = false;
  }
}

void NextIndex::updateOnSatisfying(std::size_t level)
{
  if (!d_hasLastSatisfying || level < d_lastSatisfying)
  {
    d_lastSatisfying = level;
  }
  d_hasLastSatisfying = true;

  if (!d_roundActive)
  {
    return;
  }

  if (!d_haveHigh || level < d_high)
  {
    d_high = level;
  }
  else
  {
    d_high = level;
  }
  d_haveHigh = true;
  d_step = 1;

  if (!d_haveLow || d_low >= d_high)
  {
    if (level == 0)
    {
      d_haveLow = false;
      d_low = 0;
      d_lowValidated = true;
    }
    else
    {
      d_low = level - 1;
      d_haveLow = true;
      d_lowValidated = (d_low == 0);
    }
  }
}

void NextIndex::finishRound() { d_roundActive = false; }

bool NextIndex::needsRefinement() const
{
  if (!d_roundActive)
  {
    return false;
  }
  if (!d_haveHigh)
  {
    return true;
  }
  if (!d_haveLow)
  {
    return true;
  }
  if (!d_lowValidated)
  {
    return true;
  }
  return d_low + 1 < d_high;
}

bool NextIndex::hasSatisfying() const { return d_haveHigh; }

std::size_t NextIndex::currentSatisfying() const { return d_high; }

std::size_t NextIndex::currentLower() const { return d_low; }

