#define main full_q3_frontier_certificate_main_disabled
#include "post_m29_bc_twentyfifth_frontier_gmp.cpp"
#undef main

#include <algorithm>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

namespace oracle {

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

bool is_horizontal_two_strip(const Partition& lambda, const Partition& nu) {
  int lambda_size = 0;
  int nu_size = 0;
  for (int part : lambda) lambda_size += part;
  for (int part : nu) nu_size += part;
  if (lambda_size - nu_size != 2) return false;
  const std::size_t rows = std::max(lambda.size(), nu.size());
  for (std::size_t row = 0; row < rows; ++row) {
    const int lambda_row = row < lambda.size() ? lambda[row] : 0;
    const int nu_row = row < nu.size() ? nu[row] : 0;
    const int lambda_next = row + 1 < lambda.size() ? lambda[row + 1] : 0;
    if (!(lambda_row >= nu_row && nu_row >= lambda_next)) return false;
  }
  return true;
}

}  // namespace oracle

int main() {
  std::cout << std::unitbuf;
  int failures = 0;
  int partitions_checked = 0;
  int predecessor_sets_checked = 0;

  for (int size = 2; size <= 18; ++size) {
    const std::vector<oracle::Partition> larger = oracle::partitions(size);
    const std::vector<oracle::Partition> smaller = oracle::partitions(size - 2);
    for (const oracle::Partition& lambda : larger) {
      ++partitions_checked;
      std::set<oracle::Partition> expected;
      for (const oracle::Partition& nu : smaller) {
        if (oracle::is_horizontal_two_strip(lambda, nu)) expected.insert(nu);
      }
      std::set<oracle::Partition> actual;
      for (const std::string& key : previous_horizontal_two_strip_shapes(encode(lambda))) {
        actual.insert(decode(key));
      }
      ++predecessor_sets_checked;
      if (actual != expected) {
        ++failures;
        std::cout << "KERNEL_FAIL size=" << size << " lambda=";
        for (int part : lambda) std::cout << part << ',';
        std::cout << " actual=" << actual.size()
                  << " expected=" << expected.size() << "\n";
      }
    }
  }

  const oracle::Partition largest_live_label = {121, 120, 3, 1};
  if (decode(encode(largest_live_label)) != largest_live_label) ++failures;
  if (partitions_checked != 1595 || predecessor_sets_checked != 1595) ++failures;

  std::cout << "FRONTIER_HORIZONTAL_STRIP_KERNEL partitions_checked="
            << partitions_checked
            << " predecessor_sets_checked=" << predecessor_sets_checked
            << " maximum_live_part=121 encoding_limit=255 failures=" << failures
            << "\n";
  if (failures == 0) {
    std::cout << "FRONTIER_HORIZONTAL_STRIP_KERNEL VERIFICATION: ALL PASS\n";
  }
  std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << "\n";
  return failures == 0 ? 0 : 1;
}
