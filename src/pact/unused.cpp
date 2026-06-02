// Legacy helpers that were used in previous prototype implementations of the
// approximate counting procedure. They are intentionally excluded from the
// build but kept in the repository for reference.

#if 0
#include <cvc5/cvc5.h>
#include <random>
#include <vector>

cvc5::Term legacy_generate_xor_hash(cvc5::Solver& solver,
                                    const std::vector<cvc5::Term>& vars,
                                    std::mt19937& rng)
{
  using namespace cvc5;
  Term parity = solver.mkBitVector(1, 0u);
  std::uniform_int_distribution<int> bit(0, 1);

  for (const Term& v : vars)
  {
    if (v.getSort().isBoolean())
    {
      if (bit(rng))
      {
        Term bv = solver.mkTerm(Kind::ITE,
                                {v,
                                 solver.mkBitVector(1, 1u),
                                 solver.mkBitVector(1, 0u)});
        parity = solver.mkTerm(Kind::BITVECTOR_XOR, {parity, bv});
      }
    }
    else if (v.getSort().isBitVector())
    {
      uint32_t w = v.getSort().getBitVectorSize();
      for (uint32_t i = 0; i < w; ++i)
      {
        if (bit(rng))
        {
          std::vector<uint32_t> params = {i, i};
          Op ext = solver.mkOp(Kind::BITVECTOR_EXTRACT, params);
          Term bt = solver.mkTerm(ext, {v});
          parity = solver.mkTerm(Kind::BITVECTOR_XOR, {parity, bt});
        }
      }
    }
  }

  Term rhs = solver.mkBitVector(1, bit(rng));
  return solver.mkTerm(Kind::EQUAL, {parity, rhs});
}


// Return next index to look at for count, negative denotes that we have found
// the count in the last iteration, done
int64_t SmtApproxMc::getNextIndex(uint64_t prev_index,
                                  uint64_t prev_prev_index,
                                  uint64_t count_i,
                                  bool start_of_iter)
{
  hash_cnt = prev_index;
  hash_prev = prev_prev_index;
  bool is_first_iter = false;

  if (count_i == threshold + 1)
  {
    threshold_sols[hash_cnt] = 1;
    Trace("pact-getnext") << "Threshold solutions found at " << hash_cnt
                          << "\n";
  }
  else
  {
    threshold_sols[hash_cnt] = 0;
    Trace("pact-getnext") << "Threshold solutions not found at " << hash_cnt
                          << "\n";
  }

  // We are doing a galloping search here (see our IJCAI-16 paper for more
  // details). lowerFib is referred to as loIndex and upperFib is referred to as
  // hiIndex The key idea is that we first do an exponential search and then do
  // binary search This is implemented by using two sentinels: lowerFib and
  // upperFib. The correct answer
  //  is always between lowFib and upperFib. We do exponential search until
  //  upperFib < lowerFib/2 Once upperFib < lowerFib/2; we do a binary search.
  // while (num_explored < total_max_hashes)
  // {

  if (lower_fib == 0 && upper_fib == total_max_hashes && hash_cnt == 0)
  {
    Assert(prev_index == 0);
    Trace("pact-getnext") << "getNextIndex returning 0 as first iteration\n";
    upper_fib = total_max_hashes - 1;
    is_first_iter = true;
    if (prev_measure != 0)
    {
      Assert(prev_measure != 0) << "prev_measure is 0, which is unlikely";
      Trace("pact-getnext")
          << "Returning prev_measure " << prev_measure << "\n";
    }
  }

  uint64_t cur_hash_cnt = hash_cnt;
  // prev_measure = prev_index;  // TODO AS : not sure

  const uint64_t num_sols = std::min<uint64_t>(count, threshold + 1);
  Assert(num_sols <= threshold + 1);
  // bool found_full = (num_sols == threshold + 1);

  if (num_sols < threshold + 1)
  {
    num_explored = lower_fib + total_max_hashes - hash_cnt;
    Trace("pact-getnext") << "Found " << num_sols << " solutions at "
                          << hash_cnt << "threshold " << threshold << "\n";

    // one less hash count had threshold solutions
    // this one has less than threshold
    // so this is the real deal!
    if (hash_cnt == 0
        || (threshold_sols.find(hash_cnt - 1) != threshold_sols.end()
            && threshold_sols[hash_cnt - 1] == 1))
    {
      num_hash_list.push_back(hash_cnt);
      num_count_list.push_back(num_sols);
      prev_measure = hash_cnt;
      Trace("pact-getnext")
          << "Found winner (595)" << num_sols << " solutions at " << hash_cnt
          << prev_measure << "\n";
      return -(hash_cnt);
    }

    threshold_sols[hash_cnt] = 0;
    sols_for_hash[hash_cnt] = num_sols;

    Trace("pact-getnext") << "At this point " << prev_measure << " " << hash_cnt
                          << " " << upper_fib << " " << iter << "\n";
    if (std::abs(hash_cnt - prev_measure) <= 2)
    {
      // Doing linear, this is a re-count
      upper_fib = hash_cnt;
      hash_cnt--;
      Trace("pact-getnext") << "Linear search recount : " << lower_fib << " "
                            << hash_cnt << " " << upper_fib << "\n";
    }
    else
    {
      if (hash_prev > hash_cnt) hash_prev = 0;
      upper_fib = hash_cnt;
      if (hash_prev > lower_fib) lower_fib = hash_prev;
      hash_cnt = (upper_fib + lower_fib) / 2;
      Trace("pact-getnext")
          << "Binary search (found < threshold): " << lower_fib << " "
          << hash_cnt << " " << upper_fib << "\n";
    }
  }
  else
  {
    Assert(num_sols == threshold + 1);
    Trace("pact-getnext") << "Found full solutions at " << hash_cnt << "\n";
    num_explored = hash_cnt + total_max_hashes - upper_fib;

    // success record for +1 hashcount exists and is 0
    // so one-above hashcount was below threshold, this is above
    // we have a winner -- the one above!
    if (threshold_sols.find(hash_cnt + 1) != threshold_sols.end()
        && threshold_sols[hash_cnt + 1] == 0 && !is_first_iter)
    {
      num_hash_list.push_back(hash_cnt + 1);
      num_count_list.push_back(sols_for_hash[hash_cnt + 1]);
      prev_measure = hash_cnt + 1;
      Trace("pact-getnext")
          << "Found winner (617)" << sols_for_hash[hash_cnt + 1]
          << " solutions at " << hash_cnt + 1 << "\n";
      return -(hash_cnt + 1);
    }

    threshold_sols[hash_cnt] = 1;
    sols_for_hash[hash_cnt] = threshold + 1;
    if (iter > 0 && std ::abs(hash_cnt - prev_measure) < 2)
    {
      // Doing linear, this is a re-count
      lower_fib = hash_cnt;
      hash_cnt++;
      Trace("pact-getnext") << "Linear search: " << lower_fib << " " << hash_cnt
                            << " " << upper_fib << "\n";
    }
    else if (lower_fib + (hash_cnt - lower_fib) * 2 >= upper_fib - 1)
    {
      // Whenever the above condition is satisfied, we are in binary search
      // mode
      lower_fib = hash_cnt;
      hash_cnt = (lower_fib + upper_fib) / 2;
      Trace("pact-getnext") << "Binary search: " << lower_fib << " " << hash_cnt
                            << " " << upper_fib << "\n";
    }
    else
    {
      // We are in exponential search mode.
      const auto old_hash_cnt = hash_cnt;
      hash_cnt = lower_fib + (hash_cnt - lower_fib) * 2;
      if (old_hash_cnt == hash_cnt) hash_cnt++;
      Trace("pact-getnext")
          << "Exponential search: " << lower_fib << " " << old_hash_cnt << " "
          << upper_fib << " " << hash_cnt << "\n";
    }
  }
  Trace("pact-getnext") << "Returning hash count: " << hash_cnt << "\n";
  hash_prev = cur_hash_cnt;
  if (is_first_iter)
  {
    if (prev_measure != 0)
      return prev_measure;
    else
      return 0;
  }
  return hash_cnt;
}
void SmtApproxMc::init_iteration_data()
{
  num_hash_list.clear();
  num_count_list.clear();
  threshold_sols.clear();
  sols_for_hash.clear();
  // TODO AS max hashes is not the right term here
  total_max_hashes = projection_vars.size()
                     * 32;  // should be getmaxBW or something complicated
  lower_fib = 0;
  upper_fib = total_max_hashes;
  hash_cnt = 0;
  hash_prev = 0;

  threshold = getPivot();
  count = threshold + 1;
  // TODO AS remember to add the actual constant generation routine
  threshold_sols[total_max_hashes] = 0;
  sols_for_hash[total_max_hashes] = 1;
}
#endif

