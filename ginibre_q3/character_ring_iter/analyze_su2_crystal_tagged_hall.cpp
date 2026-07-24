#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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

struct HighestPath {
  int weight;
  std::vector<int> elements;
  std::vector<int> heights;
};

static void enumerate_highest_rec(const std::vector<int>& labels,
                                  std::size_t at, int height,
                                  std::vector<int>& elements,
                                  std::vector<int>& heights,
                                  std::vector<HighestPath>& out) {
  if (at == labels.size()) {
    out.push_back({height, elements, heights});
    return;
  }
  const int p = labels[at];
  for (int element = 0; element <= std::min(p, height); ++element) {
    const int next = height + p - 2 * element;
    elements.push_back(element);
    heights.push_back(next);
    enumerate_highest_rec(labels, at + 1, next, elements, heights, out);
    heights.pop_back();
    elements.pop_back();
  }
}

static std::vector<HighestPath> enumerate_highest(
    const std::vector<int>& labels) {
  std::vector<HighestPath> out;
  std::vector<int> elements;
  std::vector<int> heights = {0};
  enumerate_highest_rec(labels, 0, 0, elements, heights, out);
  return out;
}

struct Bin {
  std::uint64_t retained_mask;
  std::vector<int> retained_elements;
  int x_weight;

  bool operator<(const Bin& other) const {
    return std::tie(retained_mask, retained_elements, x_weight)
           < std::tie(other.retained_mask, other.retained_elements,
                      other.x_weight);
  }
};

static std::vector<int> labels_on_mask(const std::vector<int>& labels,
                                       std::uint64_t mask) {
  std::vector<int> out;
  for (int i = 0; i < static_cast<int>(labels.size()); ++i)
    if (mask & (std::uint64_t{1} << i))
      out.push_back(labels[static_cast<std::size_t>(i)]);
  return out;
}

struct HallResult {
  bool overload = false;
  Bin bin{};
  std::uint64_t demand = 0;
  std::uint64_t capacity = 0;
  std::uint64_t global_demand = 0;
  std::uint64_t global_capacity = 0;
  std::uint64_t sources = 0;
  std::uint64_t bins = 0;
};

static HallResult check_instance(const std::vector<int>& labels, int q) {
  const int n = static_cast<int>(labels.size());
  const std::uint64_t full_mask = (std::uint64_t{1} << n) - 1;
  std::map<Bin, std::uint64_t> demand;
  std::map<int, std::uint64_t> demand_by_weight;
  std::map<int, std::uint64_t> capacity_by_weight;
  HallResult result;
  for (std::uint64_t y_mask = 0; y_mask <= full_mask; ++y_mask) {
    const std::uint64_t x_mask = full_mask ^ y_mask;
    const auto y_paths = enumerate_highest(labels_on_mask(labels, y_mask));
    const auto x_paths = enumerate_highest(labels_on_mask(labels, x_mask));
    std::vector<int> y_ids;
    for (int i = 0; i < n; ++i)
      if (y_mask & (std::uint64_t{1} << i)) y_ids.push_back(i);
    for (const auto& y_path : y_paths) {
      if (y_path.weight != q) continue;
      int last_zero = 0;
      for (int i = 0; i < static_cast<int>(y_path.heights.size()); ++i)
        if (y_path.heights[static_cast<std::size_t>(i)] == 0)
          last_zero = i;
      std::uint64_t retained_mask = 0;
      std::vector<int> retained_elements;
      for (int i = 0; i < last_zero; ++i) {
        retained_mask |= std::uint64_t{1}
                         << y_ids[static_cast<std::size_t>(i)];
        retained_elements.push_back(
            y_path.elements[static_cast<std::size_t>(i)]);
      }
      for (const auto& x_path : x_paths) {
        ++demand[{retained_mask, retained_elements, x_path.weight}];
        ++demand_by_weight[x_path.weight];
        ++result.sources;
      }
    }
  }
  for (std::uint64_t y_mask = 0; y_mask <= full_mask; ++y_mask) {
    std::uint64_t invariant_paths = 0;
    for (const auto& path :
         enumerate_highest(labels_on_mask(labels, y_mask)))
      if (path.weight == 0) ++invariant_paths;
    if (invariant_paths == 0) continue;
    std::vector<int> target_labels = {q};
    const auto x_labels = labels_on_mask(labels, full_mask ^ y_mask);
    target_labels.insert(target_labels.end(), x_labels.begin(),
                         x_labels.end());
    for (const auto& path : enumerate_highest(target_labels))
      capacity_by_weight[path.weight] += invariant_paths;
  }
  result.bins = demand.size();
  for (const auto& [bin, count] : demand) {
    std::vector<int> target_labels = {q};
    const std::uint64_t target_x_mask = full_mask ^ bin.retained_mask;
    const auto suffix = labels_on_mask(labels, target_x_mask);
    target_labels.insert(target_labels.end(), suffix.begin(), suffix.end());
    std::uint64_t capacity = 0;
    for (const auto& path : enumerate_highest(target_labels))
      if (path.weight == bin.x_weight) ++capacity;
    if (count > capacity) {
      result.overload = true;
      result.bin = bin;
      result.demand = count;
      result.capacity = capacity;
      result.global_demand = demand_by_weight[bin.x_weight];
      result.global_capacity = capacity_by_weight[bin.x_weight];
      return result;
    }
  }
  return result;
}

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
  const int max_factors = argc > 1 ? std::stoi(argv[1]) : 7;
  const int max_label = argc > 2 ? std::stoi(argv[2]) : 4;
  std::uint64_t instances = 0;
  std::uint64_t sources = 0;
  std::uint64_t bins = 0;
  for (int n = 1; n <= max_factors; ++n) {
    std::vector<int> labels(static_cast<std::size_t>(n), 1);
    do {
      int total = 0;
      for (int p : labels) total += p;
      for (int q = 1; q <= total; ++q) {
        const HallResult result = check_instance(labels, q);
        ++instances;
        sources += result.sources;
        bins += result.bins;
        if (result.overload) {
          std::cout << "CRYSTAL_TAGGED_HALL OVERLOAD"
                    << " labels=" << vector_string(labels)
                    << " q=" << q
                    << " retained_mask=" << result.bin.retained_mask
                    << " retained_elements="
                    << vector_string(result.bin.retained_elements)
                    << " x_weight=" << result.bin.x_weight
                    << " demand=" << result.demand
                    << " capacity=" << result.capacity
                    << " global_demand=" << result.global_demand
                    << " global_capacity=" << result.global_capacity
                    << " instances=" << instances << "\n";
          return 0;
        }
      }
    } while (increment_labels(labels, max_label));
  }
  std::cout << "CRYSTAL_TAGGED_HALL NO_OVERLOAD"
            << " max_factors=" << max_factors
            << " max_label=" << max_label
            << " instances=" << instances
            << " sources=" << sources
            << " bins=" << bins << "\n";
}
