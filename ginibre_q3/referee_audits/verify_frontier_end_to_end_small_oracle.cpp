// Referee-side end-to-end oracle for the H25--H28 B-frontier computation.
//
// The production program starts at doubled terminal partitions and propagates
// backward with its optimized predecessor generator.  This audit instead
// enumerates every partition of 0,2,...,2*MAX_J and propagates forward using
// the defining interlacing criterion for a horizontal two-strip.  It thereby
// checks terminal generation, the doubled-partition boundary condition,
// iteration count, path multiplicities, and final aggregation in addition to
// the already-audited local predecessor kernel.

// Compile against any one of the four live sources, for example:
//   g++ -O2 -std=c++20 -fopenmp [FRONTIER_SOURCE definition]
//     verify_frontier_end_to_end_small_oracle.cpp -lgmpxx -lgmp

#ifndef FRONTIER_SOURCE
#define FRONTIER_SOURCE \
  "../character_ring_iter/post_m29_bc_twentyfifth_frontier_gmp.cpp"
#endif

#define main full_q3_frontier_certificate_main_disabled
#include FRONTIER_SOURCE
#undef main

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace independent {

using Partition = std::vector<int>;

void generate_partitions(
    int remaining,
    int maximum,
    Partition& prefix,
    std::vector<Partition>& output
) {
  if (remaining == 0) {
    output.push_back(prefix);
    return;
  }
  for (int part = std::min(remaining, maximum); part >= 1; --part) {
    prefix.push_back(part);
    generate_partitions(remaining - part, part, prefix, output);
    prefix.pop_back();
  }
}

std::vector<Partition> partitions(int size) {
  std::vector<Partition> output;
  Partition prefix;
  generate_partitions(size, size, prefix, output);
  return output;
}

bool is_horizontal_two_strip(
    const Partition& larger,
    const Partition& smaller
) {
  int larger_size = 0;
  int smaller_size = 0;
  for (int part : larger) larger_size += part;
  for (int part : smaller) smaller_size += part;
  if (larger_size - smaller_size != 2) return false;

  const std::size_t rows = std::max(larger.size(), smaller.size());
  for (std::size_t row = 0; row < rows; ++row) {
    const int larger_row = row < larger.size() ? larger[row] : 0;
    const int smaller_row = row < smaller.size() ? smaller[row] : 0;
    const int larger_next =
        row + 1 < larger.size() ? larger[row + 1] : 0;
    if (!(larger_row >= smaller_row && smaller_row >= larger_next)) {
      return false;
    }
  }
  return true;
}

Partition doubled_rows(const Partition& mu) {
  Partition lambda;
  lambda.reserve(2 * mu.size());
  for (int part : mu) {
    lambda.push_back(part);
    lambda.push_back(part);
  }
  return lambda;
}

}  // namespace independent

int main() {
  constexpr int MAX_J = 14;
  std::cout << std::unitbuf;

  std::map<independent::Partition, mpz_class> forward;
  forward[independent::Partition{}] = 1;
  int failures = 0;
  int cases_checked = 0;
  int endpoint_shapes_checked = 0;
  int transition_pairs_checked = 0;

  for (int j = 1; j <= MAX_J; ++j) {
    const std::vector<independent::Partition> next_shapes =
        independent::partitions(2 * j);
    std::map<independent::Partition, mpz_class> next;
    for (const independent::Partition& larger : next_shapes) {
      mpz_class count = 0;
      for (const auto& [smaller, multiplicity] : forward) {
        ++transition_pairs_checked;
        if (independent::is_horizontal_two_strip(larger, smaller)) {
          count += multiplicity;
        }
      }
      if (count != 0) next.emplace(larger, count);
    }
    forward = std::move(next);

    const std::vector<independent::Partition> mus =
        independent::partitions(j);
    for (int q = 1; q <= j / 2 + 1; ++q) {
      mpz_class expected = 0;
      int endpoints = 0;
      for (const independent::Partition& mu : mus) {
        if (mu.front() < 2 * q) continue;
        const independent::Partition lambda = independent::doubled_rows(mu);
        const auto found = forward.find(lambda);
        if (found != forward.end()) expected += found->second;
        ++endpoints;
      }
      endpoint_shapes_checked += endpoints;

      // Case and compute_b_width_bad_count come from the selected production
      // source above.  The production route is deliberately the opposite
      // direction from the independent forward enumeration.
      const mpz_class actual = compute_b_width_bad_count(Case{q, j});
      ++cases_checked;
      if (actual != expected) {
        ++failures;
        std::cout << "FRONTIER_END_TO_END_FAIL j=" << j
                  << " q=" << q
                  << " actual=" << actual
                  << " expected=" << expected << "\n";
      }
    }
  }

  if (cases_checked != 63) ++failures;
  std::cout << "FRONTIER_END_TO_END_SMALL_ORACLE"
            << " source=" << FRONTIER_SOURCE
            << " max_j=" << MAX_J
            << " cases_checked=" << cases_checked
            << " endpoint_shapes_checked=" << endpoint_shapes_checked
            << " transition_pairs_checked=" << transition_pairs_checked
            << " failures=" << failures << "\n";
  if (failures == 0) {
    std::cout << "FRONTIER_END_TO_END_SMALL_ORACLE VERIFICATION: ALL PASS\n";
  }
  return failures == 0 ? 0 : 1;
}
