#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

struct Item {
  int id;
  int label;
  int element;
};

static int tensor_epsilon(int a, int i, int j) {
  return i + std::max(0, j - (a - i));
}

static int tensor_weight(int a, int i, int b, int j) {
  return a - 2 * i + b - 2 * j;
}

static std::pair<int, int> commute_pair(int a, int i, int b, int j) {
  const int epsilon = tensor_epsilon(a, i, j);
  const int weight = tensor_weight(a, i, b, j);
  std::optional<std::pair<int, int>> answer;
  for (int jp = 0; jp <= b; ++jp) {
    for (int ip = 0; ip <= a; ++ip) {
      if (tensor_epsilon(b, jp, ip) != epsilon
          || tensor_weight(b, jp, a, ip) != weight)
        continue;
      if (answer) {
        std::cerr << "PAIR_COMMUTOR_NOT_UNIQUE\n";
        std::exit(2);
      }
      answer = {jp, ip};
    }
  }
  if (!answer) {
    std::cerr << "PAIR_COMMUTOR_MISSING\n";
    std::exit(2);
  }
  return *answer;
}

static void swap_adjacent(std::vector<Item>& items, int left) {
  Item& x = items[static_cast<std::size_t>(left)];
  Item& y = items[static_cast<std::size_t>(left + 1)];
  const auto [new_y, new_x] =
      commute_pair(x.label, x.element, y.label, y.element);
  const Item old_x = x;
  const Item old_y = y;
  x = {old_y.id, old_y.label, new_y};
  y = {old_x.id, old_x.label, new_x};
}

static void reorder(std::vector<Item>& items,
                    const std::vector<int>& desired_ids) {
  for (int target = 0; target < static_cast<int>(desired_ids.size());
       ++target) {
    int found = target;
    while (found < static_cast<int>(items.size())
           && items[static_cast<std::size_t>(found)].id
                  != desired_ids[static_cast<std::size_t>(target)])
      ++found;
    if (found == static_cast<int>(items.size())) {
      std::cerr << "PERMUTATION_ID_MISSING\n";
      std::exit(2);
    }
    while (found > target) {
      swap_adjacent(items, found - 1);
      --found;
    }
  }
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

struct Source {
  std::uint64_t y_mask;
  std::vector<int> y_elements;
  std::vector<int> x_elements;
};

static std::string source_string(const Source& source) {
  std::ostringstream out;
  out << "ymask=" << source.y_mask
      << " y=" << vector_string(source.y_elements)
      << " x=" << vector_string(source.x_elements);
  return out.str();
}

static std::string target_key(std::uint64_t y_mask,
                              const std::vector<int>& y_elements,
                              const std::vector<Item>& x_items) {
  std::ostringstream out;
  out << y_mask << ":";
  for (int value : y_elements) out << value << ",";
  out << ":";
  for (const Item& item : x_items)
    out << item.id << "/" << item.element << ",";
  return out.str();
}

struct Search {
  std::vector<int> labels;
  int q;
  std::map<std::string, Source> images;
  std::uint64_t sources = 0;
  bool collision = false;
  Source collision_first;
  Source collision_second;
  std::string collision_key;

  Search(std::vector<int> input_labels, int input_q)
      : labels(std::move(input_labels)), q(input_q) {}

  void enumerate_x(const std::vector<int>& x_ids, std::size_t at,
                   std::vector<int>& x_elements, std::uint64_t y_mask,
                   const std::vector<int>& y_ids,
                   const std::vector<int>& y_elements,
                   const std::vector<int>& heights) {
    if (collision) return;
    if (at != x_ids.size()) {
      const int p = labels[static_cast<std::size_t>(
          x_ids[at])];
      for (int element = 0; element <= p; ++element) {
        x_elements.push_back(element);
        enumerate_x(x_ids, at + 1, x_elements, y_mask, y_ids, y_elements,
                    heights);
        x_elements.pop_back();
      }
      return;
    }

    int last_zero = 0;
    for (int i = 0; i < static_cast<int>(heights.size()); ++i)
      if (heights[static_cast<std::size_t>(i)] == 0) last_zero = i;

    std::uint64_t target_y_mask = 0;
    std::vector<int> target_y_elements;
    for (int i = 0; i < last_zero; ++i) {
      target_y_mask |= std::uint64_t{1}
                       << y_ids[static_cast<std::size_t>(i)];
      target_y_elements.push_back(
          y_elements[static_cast<std::size_t>(i)]);
    }

    std::vector<Item> items;
    items.push_back({-1, q, 0});
    for (int i = static_cast<int>(y_ids.size()) - 1; i >= last_zero; --i) {
      const int previous_height =
          heights[static_cast<std::size_t>(i + 1)];
      const int next_height = heights[static_cast<std::size_t>(i)];
      const int id = y_ids[static_cast<std::size_t>(i)];
      const int p = labels[static_cast<std::size_t>(id)];
      const int element = (previous_height + p - next_height) / 2;
      items.push_back({id, p, element});
    }
    for (std::size_t i = 0; i < x_ids.size(); ++i) {
      const int id = x_ids[i];
      items.push_back(
          {id, labels[static_cast<std::size_t>(id)], x_elements[i]});
    }

    std::vector<int> desired_ids = {-1};
    for (int id = 0; id < static_cast<int>(labels.size()); ++id)
      if ((target_y_mask & (std::uint64_t{1} << id)) == 0)
        desired_ids.push_back(id);
    reorder(items, desired_ids);

    const Source source = {y_mask, y_elements, x_elements};
    const std::string key =
        target_key(target_y_mask, target_y_elements, items);
    const auto [it, inserted] = images.emplace(key, source);
    ++sources;
    if (!inserted) {
      collision = true;
      collision_first = it->second;
      collision_second = source;
      collision_key = key;
    }
  }

  void enumerate_y_elements(std::uint64_t y_mask,
                            const std::vector<int>& y_ids,
                            const std::vector<int>& x_ids, std::size_t at,
                            int height, std::vector<int>& y_elements,
                            std::vector<int>& heights) {
    if (collision) return;
    if (at == y_ids.size()) {
      if (height != q) return;
      std::vector<int> x_elements;
      enumerate_x(x_ids, 0, x_elements, y_mask, y_ids, y_elements, heights);
      return;
    }
    const int p =
        labels[static_cast<std::size_t>(y_ids[at])];
    for (int element = 0; element <= std::min(p, height); ++element) {
      const int next_height = height + p - 2 * element;
      y_elements.push_back(element);
      heights.push_back(next_height);
      enumerate_y_elements(y_mask, y_ids, x_ids, at + 1, next_height,
                           y_elements, heights);
      heights.pop_back();
      y_elements.pop_back();
    }
  }

  void run() {
    const int n = static_cast<int>(labels.size());
    for (std::uint64_t y_mask = 0;
         y_mask < (std::uint64_t{1} << n) && !collision; ++y_mask) {
      std::vector<int> y_ids;
      std::vector<int> x_ids;
      for (int id = 0; id < n; ++id) {
        if (y_mask & (std::uint64_t{1} << id)) y_ids.push_back(id);
        else x_ids.push_back(id);
      }
      std::vector<int> y_elements;
      std::vector<int> heights = {0};
      enumerate_y_elements(y_mask, y_ids, x_ids, 0, 0, y_elements, heights);
    }
  }
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
  const int max_factors = argc > 1 ? std::stoi(argv[1]) : 6;
  const int max_label = argc > 2 ? std::stoi(argv[2]) : 4;
  std::uint64_t instances = 0;
  std::uint64_t sources = 0;
  for (int n = 1; n <= max_factors; ++n) {
    std::vector<int> labels(static_cast<std::size_t>(n), 1);
    do {
      int total = 0;
      for (int p : labels) total += p;
      for (int q = 1; q <= total; ++q) {
        Search search(labels, q);
        search.run();
        ++instances;
        sources += search.sources;
        if (search.collision) {
          std::cout << "CRYSTAL_LAST_RETURN COLLISION"
                    << " labels=" << vector_string(labels)
                    << " q=" << q
                    << " first={" << source_string(search.collision_first)
                    << "} second={" << source_string(search.collision_second)
                    << "} target=" << search.collision_key << "\n";
          return 0;
        }
      }
    } while (increment_labels(labels, max_label));
  }
  std::cout << "CRYSTAL_LAST_RETURN NO_COLLISION"
            << " max_factors=" << max_factors
            << " max_label=" << max_label
            << " instances=" << instances
            << " sources=" << sources << "\n";
}
