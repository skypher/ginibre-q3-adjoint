#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using Character = std::map<int, cpp_int>;

static Character multiply(const Character& input, int p) {
  Character output;
  for (const auto& [r, coefficient] : input)
    for (int s = std::abs(r - p); s <= r + p; s += 2)
      output[s] += coefficient;
  return output;
}

static cpp_int coefficient(const Character& character, int label) {
  const auto it = character.find(label);
  return it == character.end() ? cpp_int(0) : it->second;
}

static std::string vector_string(const std::vector<int>& values) {
  std::ostringstream out;
  out << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) out << ",";
    out << values[i];
  }
  out << "]";
  return out.str();
}

static bool is_primitive_invariant(std::uint64_t mask,
                                   const std::vector<Character>& products) {
  if (mask == 0) return true;
  if (coefficient(products[static_cast<std::size_t>(mask)], 0) == 0)
    return false;
  for (std::uint64_t proper = (mask - 1) & mask; proper != 0;
       proper = (proper - 1) & mask)
    if (coefficient(products[static_cast<std::size_t>(proper)], 0) > 0)
      return false;
  return true;
}

struct Witness {
  bool found = false;
  std::vector<int> labels;
  int q = 0;
  int a = 0;
  cpp_int demand = 0;
  cpp_int primitive_capacity = 0;
  cpp_int full_capacity = 0;
};

static bool increment_labels(std::vector<int>& labels, int max_label) {
  for (int i = static_cast<int>(labels.size()) - 1; i >= 0; --i) {
    if (labels[static_cast<std::size_t>(i)] < max_label) {
      const int next = labels[static_cast<std::size_t>(i)] + 1;
      for (int j = i; j < static_cast<int>(labels.size()); ++j)
        labels[static_cast<std::size_t>(j)] = next;
      return true;
    }
  }
  return false;
}

int main(int argc, char** argv) {
  const int max_factors = argc > 1 ? std::stoi(argv[1]) : 8;
  const int max_label = argc > 2 ? std::stoi(argv[2]) : 5;
  std::uint64_t lists = 0;
  std::uint64_t entries = 0;
  Witness witness;
  for (int n = 1; n <= max_factors && !witness.found; ++n) {
    std::vector<int> labels(static_cast<std::size_t>(n), 1);
    do {
      ++lists;
      const std::uint64_t subset_count = std::uint64_t{1} << n;
      const std::uint64_t full_mask = subset_count - 1;
      std::vector<Character> products(static_cast<std::size_t>(subset_count));
      products[0][0] = 1;
      for (std::uint64_t mask = 1; mask < subset_count; ++mask) {
        const int bit = __builtin_ctzll(mask);
        const std::uint64_t previous = mask & (mask - 1);
        products[static_cast<std::size_t>(mask)] =
            multiply(products[static_cast<std::size_t>(previous)],
                     labels[static_cast<std::size_t>(bit)]);
      }
      int total = 0;
      for (int p : labels) total += p;
      for (int q = 1; q <= total && !witness.found; ++q) {
        for (int a = 0; a <= total && !witness.found; ++a) {
          cpp_int demand = 0;
          cpp_int primitive_capacity = 0;
          cpp_int full_capacity = 0;
          for (std::uint64_t mask = 0; mask < subset_count; ++mask) {
            const std::uint64_t complement = full_mask ^ mask;
            demand += coefficient(products[static_cast<std::size_t>(mask)], q)
                      * coefficient(
                          products[static_cast<std::size_t>(complement)], a);
            const cpp_int invariant =
                coefficient(products[static_cast<std::size_t>(mask)], 0);
            if (invariant == 0) continue;
            cpp_int endpoint_capacity = 0;
            for (int t = std::abs(a - q); t <= a + q; t += 2)
              endpoint_capacity += coefficient(
                  products[static_cast<std::size_t>(complement)], t);
            const cpp_int contribution = invariant * endpoint_capacity;
            full_capacity += contribution;
            if (is_primitive_invariant(mask, products))
              primitive_capacity += contribution;
          }
          ++entries;
          if (primitive_capacity < demand && demand <= full_capacity) {
            witness = {true, labels, q, a, demand, primitive_capacity,
                       full_capacity};
          }
        }
      }
    } while (!witness.found && increment_labels(labels, max_label));
  }
  if (witness.found) {
    std::cout << "BOUNDARY_PRIMITIVE_SUPPLIERS OBSTRUCTION"
              << " labels=" << vector_string(witness.labels)
              << " q=" << witness.q
              << " a=" << witness.a
              << " demand=" << witness.demand
              << " primitive_capacity=" << witness.primitive_capacity
              << " full_capacity=" << witness.full_capacity
              << " lists=" << lists
              << " entries=" << entries << "\n";
  } else {
    std::cout << "BOUNDARY_PRIMITIVE_SUPPLIERS NO_OBSTRUCTION"
              << " max_factors=" << max_factors
              << " max_label=" << max_label
              << " lists=" << lists
              << " entries=" << entries << "\n";
  }
}
