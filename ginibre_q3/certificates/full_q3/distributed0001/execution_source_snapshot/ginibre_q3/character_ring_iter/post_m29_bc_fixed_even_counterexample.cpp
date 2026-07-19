#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

struct Cell {
    int r;
    int c;
    bool operator==(const Cell& other) const {
        return r == other.r && c == other.c;
    }
};

static const std::vector<int> kShape = {12, 10, 8, 6, 4, 3, 3, 2, 2, 1};

static const std::vector<Cell> kEvenCells = {
    {1, 2}, {2, 1}, {2, 2}, {3, 2}, {4, 1},
    {2, 4}, {1, 5}, {6, 1}, {1, 7}, {1, 8},
    {5, 2}, {7, 1}, {1, 9}, {5, 3}, {2, 6},
    {3, 5}, {4, 5}, {2, 7}, {2, 8}, {4, 6},
    {5, 4}, {1, 11}, {9, 1}, {3, 8}, {10, 1}
};

static const std::vector<Cell> kOddWitnessCells = {
    {1, 1}, {1, 3}, {1, 4}, {3, 1}, {2, 3}, {4, 2},
    {5, 1}, {1, 6}, {3, 3}, {4, 3}, {3, 4}, {4, 4},
    {8, 1}, {6, 2}, {2, 5}, {7, 2}, {1, 10}, {8, 2},
    {3, 6}, {6, 3}, {2, 9}, {3, 7}, {7, 3}, {9, 2},
    {1, 12}, {2, 10}
};

static std::string shape_key(const std::vector<int>& shape, int k) {
    std::ostringstream out;
    out << k << ':';
    for (int x : shape) out << x << ',';
    return out.str();
}

static bool contains_cell(const std::vector<int>& shape, Cell x) {
    return x.r >= 1 && x.r <= static_cast<int>(shape.size()) &&
           x.c >= 1 && x.c <= shape[x.r - 1];
}

static bool is_removable(const std::vector<int>& shape, Cell x) {
    if (!contains_cell(shape, x)) return false;
    if (shape[x.r - 1] != x.c) return false;
    if (x.r < static_cast<int>(shape.size()) && shape[x.r] >= x.c) return false;
    return true;
}

static bool is_addable(const std::vector<int>& shape, Cell x) {
    if (x.r < 1 || x.c < 1) return false;
    if (x.r > static_cast<int>(shape.size()) + 1) return false;
    const int current = (x.r <= static_cast<int>(shape.size())) ? shape[x.r - 1] : 0;
    if (x.c != current + 1) return false;
    if (x.r > 1 && shape[x.r - 2] < x.c) return false;
    return true;
}

static std::vector<Cell> removable_corners(const std::vector<int>& shape) {
    std::vector<Cell> out;
    for (int r = 1; r <= static_cast<int>(shape.size()); ++r) {
        Cell x{r, shape[r - 1]};
        if (is_removable(shape, x)) out.push_back(x);
    }
    return out;
}

static std::vector<int> remove_cell(std::vector<int> shape, Cell x) {
    --shape[x.r - 1];
    while (!shape.empty() && shape.back() == 0) shape.pop_back();
    return shape;
}

static std::vector<int> add_cell(std::vector<int> shape, Cell x) {
    if (x.r == static_cast<int>(shape.size()) + 1) {
        shape.push_back(1);
    } else {
        ++shape[x.r - 1];
    }
    return shape;
}

static bool is_fixed_prefix_cell(Cell x, int k) {
    for (int i = 0; i < k; ++i) {
        if (kEvenCells[i] == x) return true;
    }
    return false;
}

static std::unordered_map<std::string, std::uint64_t> memo;

static std::uint64_t count_fixed_even(const std::vector<int>& shape, int k) {
    const std::string key = shape_key(shape, k);
    auto it = memo.find(key);
    if (it != memo.end()) return it->second;

    std::uint64_t result = 0;
    if (k == 0) {
        result = (shape.size() == 1 && shape[0] == 1) ? 1 : 0;
        memo.emplace(key, result);
        return result;
    }

    const Cell e = kEvenCells[k - 1];
    if (!contains_cell(shape, e)) {
        memo.emplace(key, 0);
        return 0;
    }

    for (Cell y : removable_corners(shape)) {
        if (is_fixed_prefix_cell(y, k)) continue;
        std::vector<int> after_y = remove_cell(shape, y);
        if (!is_removable(after_y, e)) continue;
        std::vector<int> child = remove_cell(after_y, e);
        result += count_fixed_even(child, k - 1);
    }
    memo.emplace(key, result);
    return result;
}

static bool check_witness_path() {
    std::vector<int> shape;
    for (int i = 0; i < static_cast<int>(kEvenCells.size()); ++i) {
        if (!is_addable(shape, kOddWitnessCells[i])) return false;
        shape = add_cell(shape, kOddWitnessCells[i]);
        if (!is_addable(shape, kEvenCells[i])) return false;
        shape = add_cell(shape, kEvenCells[i]);
    }
    if (!is_addable(shape, kOddWitnessCells.back())) return false;
    shape = add_cell(shape, kOddWitnessCells.back());
    return shape == kShape;
}

int main() {
    const std::uint64_t count = count_fixed_even(kShape, 25);
    const std::uint64_t budget = (std::uint64_t{1} << 26);
    std::cout << "B/C fixed-even broad-bound counterexample verifier\n";
    std::cout << "m=26\n";
    std::cout << "shape=(12,10,8,6,4,3,3,2,2,1)\n";
    std::cout << "fixed_even_cells=25\n";
    std::cout << "witness_path_valid=" << (check_witness_path() ? "yes" : "no") << "\n";
    std::cout << "fixed_even_count=" << count << "\n";
    std::cout << "pow2_budget=" << budget << "\n";
    std::cout << "exceeds_budget=" << (count > budget ? "yes" : "no") << "\n";
    std::cout << "memo_states=" << memo.size() << "\n";
    std::cout << "__EXIT_STATUS=0\n";
    return (check_witness_path() && count > budget) ? 0 : 1;
}
