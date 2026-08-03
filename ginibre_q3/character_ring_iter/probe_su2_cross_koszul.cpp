#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::int64_t kPrime = 1000003;
constexpr int kOuterDegree = 6;

struct MinusState {
  int exterior = 0;
  int x = 0;
  int y_plus = 0;
  int y_minus = 0;
  int z = 0;

  bool operator<(const MinusState& other) const {
    return std::array<int, 5>{exterior, x, y_plus, y_minus, z}
           < std::array<int, 5>{other.exterior, other.x, other.y_plus,
                                other.y_minus, other.z};
  }
};

struct PlusState {
  int nontrivial = 0;
  int index = 0;

  bool operator<(const PlusState& other) const {
    return std::array<int, 2>{nontrivial, index}
           < std::array<int, 2>{other.nontrivial, other.index};
  }
};

struct State {
  std::array<int, 2> minus{};
  std::array<int, 3> plus{};

  bool operator<(const State& other) const {
    if (minus != other.minus) {
      return minus < other.minus;
    }
    return plus < other.plus;
  }
};

struct Contribution {
  State state;
  std::int64_t coefficient = 0;
};

using Matrix = std::vector<std::vector<std::int64_t>>;

std::int64_t normalize(std::int64_t value) {
  value %= kPrime;
  return value < 0 ? value + kPrime : value;
}

std::int64_t inverse(std::int64_t value) {
  std::int64_t result = 1;
  std::int64_t base = normalize(value);
  std::int64_t exponent = kPrime - 2;
  while (exponent > 0) {
    if ((exponent & 1) != 0) {
      result = normalize(result * base);
    }
    base = normalize(base * base);
    exponent >>= 1;
  }
  return result;
}

int popcount(const int value) {
  int result = 0;
  for (int current = value; current != 0; current >>= 1) {
    result += current & 1;
  }
  return result;
}

int exterior_weight(const int exterior) {
  return ((exterior & 1) != 0 ? 1 : 0)
         + ((exterior & 2) != 0 ? -1 : 0);
}

int degree(const MinusState& state) {
  return popcount(state.exterior) + state.y_plus + state.y_minus
         + 2 * state.z;
}

int weight(const MinusState& state) {
  return exterior_weight(state.exterior) + state.y_plus - state.y_minus;
}

int degree(const PlusState& state, const int label) {
  return state.nontrivial != 0 ? label : 2 * state.index;
}

int weight(const PlusState& state, const int label) {
  return state.nontrivial != 0 ? label - 2 * state.index : 0;
}

std::vector<MinusState> make_minus_states(const int label) {
  const int symmetric_degree = label - 1;
  std::vector<MinusState> result;
  for (int exterior = 0; exterior < 4; ++exterior) {
    for (int x = 0; x <= symmetric_degree; ++x) {
      for (int y_plus = 0; y_plus + x <= symmetric_degree; ++y_plus) {
        for (int y_minus = 0;
             y_minus + y_plus + x <= symmetric_degree; ++y_minus) {
          const int z = symmetric_degree - x - y_plus - y_minus;
          result.push_back({exterior, x, y_plus, y_minus, z});
        }
      }
    }
  }
  return result;
}

std::vector<PlusState> make_plus_states(const int label) {
  std::vector<PlusState> result;
  for (int scalar = 0; scalar <= label; ++scalar) {
    result.push_back({0, scalar});
  }
  for (int index = 0; index <= label; ++index) {
    result.push_back({1, index});
  }
  return result;
}

int total_degree(const State& state,
                 const std::array<std::vector<MinusState>, 2>& minus,
                 const std::array<std::vector<PlusState>, 3>& plus) {
  constexpr std::array<int, 3> labels{1, 2, 2};
  int result = 0;
  for (int factor = 0; factor < 2; ++factor) {
    result += degree(minus[static_cast<std::size_t>(factor)]
                         [static_cast<std::size_t>(state.minus[
                             static_cast<std::size_t>(factor)])]);
  }
  for (int factor = 0; factor < 3; ++factor) {
    result += degree(plus[static_cast<std::size_t>(factor)]
                         [static_cast<std::size_t>(state.plus[
                             static_cast<std::size_t>(factor)])],
                     labels[static_cast<std::size_t>(factor)]);
  }
  return result;
}

int total_weight(const State& state,
                 const std::array<std::vector<MinusState>, 2>& minus,
                 const std::array<std::vector<PlusState>, 3>& plus) {
  constexpr std::array<int, 3> labels{1, 2, 2};
  int result = 0;
  for (int factor = 0; factor < 2; ++factor) {
    result += weight(minus[static_cast<std::size_t>(factor)]
                         [static_cast<std::size_t>(state.minus[
                             static_cast<std::size_t>(factor)])]);
  }
  for (int factor = 0; factor < 3; ++factor) {
    result += weight(plus[static_cast<std::size_t>(factor)]
                         [static_cast<std::size_t>(state.plus[
                             static_cast<std::size_t>(factor)])],
                     labels[static_cast<std::size_t>(factor)]);
  }
  return result;
}

int parity(const State& state,
           const std::array<std::vector<MinusState>, 2>& minus) {
  return (popcount(minus[0][static_cast<std::size_t>(state.minus[0])]
                       .exterior)
          + popcount(minus[1][static_cast<std::size_t>(state.minus[1])]
                       .exterior))
         & 1;
}

std::map<MinusState, int> index_minus(
    const std::vector<MinusState>& states) {
  std::map<MinusState, int> result;
  for (std::size_t index = 0; index < states.size(); ++index) {
    result.emplace(states[index], static_cast<int>(index));
  }
  return result;
}

std::map<PlusState, int> index_plus(const std::vector<PlusState>& states) {
  std::map<PlusState, int> result;
  for (std::size_t index = 0; index < states.size(); ++index) {
    result.emplace(states[index], static_cast<int>(index));
  }
  return result;
}

void append(std::map<State, std::int64_t>& target, const State& state,
            const std::int64_t coefficient) {
  target[state] = normalize(target[state] + coefficient);
  if (target[state] == 0) {
    target.erase(state);
  }
}

std::map<State, std::int64_t> apply_raising(
    const State& source,
    const std::array<std::vector<MinusState>, 2>& minus,
    const std::array<std::vector<PlusState>, 3>& plus,
    const std::array<std::map<MinusState, int>, 2>& minus_index,
    const std::array<std::map<PlusState, int>, 3>& plus_index) {
  std::map<State, std::int64_t> result;
  for (int factor = 0; factor < 2; ++factor) {
    const MinusState current = minus[static_cast<std::size_t>(factor)]
                                     [static_cast<std::size_t>(source.minus[
                                         static_cast<std::size_t>(factor)])];
    if ((current.exterior & 2) != 0 && (current.exterior & 1) == 0) {
      MinusState next = current;
      next.exterior = (next.exterior & ~2) | 1;
      State target = source;
      target.minus[static_cast<std::size_t>(factor)] =
          minus_index[static_cast<std::size_t>(factor)].at(next);
      append(result, target, 1);
    }
    if (current.y_minus > 0) {
      MinusState next = current;
      --next.y_minus;
      ++next.y_plus;
      State target = source;
      target.minus[static_cast<std::size_t>(factor)] =
          minus_index[static_cast<std::size_t>(factor)].at(next);
      append(result, target, current.y_minus);
    }
  }
  for (int factor = 0; factor < 3; ++factor) {
    const PlusState current = plus[static_cast<std::size_t>(factor)]
                                  [static_cast<std::size_t>(source.plus[
                                      static_cast<std::size_t>(factor)])];
    if (current.nontrivial != 0 && current.index > 0) {
      PlusState next = current;
      --next.index;
      State target = source;
      target.plus[static_cast<std::size_t>(factor)] =
          plus_index[static_cast<std::size_t>(factor)].at(next);
      append(result, target, current.index);
    }
  }
  return result;
}

std::map<State, std::int64_t> apply_differential(
    const State& source, const std::array<std::array<int, 2>, 2>& coupling,
    const std::array<std::vector<MinusState>, 2>& minus,
    const std::array<std::map<MinusState, int>, 2>& minus_index) {
  std::map<State, std::int64_t> result;
  for (int exterior_factor = 0; exterior_factor < 2; ++exterior_factor) {
    const int prior_parity = exterior_factor == 0
        ? 0
        : popcount(minus[0][static_cast<std::size_t>(source.minus[0])]
                         .exterior)
              & 1;
    const std::int64_t tensor_sign = prior_parity == 0 ? 1 : -1;
    const MinusState exterior_state =
        minus[static_cast<std::size_t>(exterior_factor)]
             [static_cast<std::size_t>(source.minus[
                 static_cast<std::size_t>(exterior_factor)])];
    for (int symmetric_factor = 0; symmetric_factor < 2;
         ++symmetric_factor) {
      const int scalar = coupling[static_cast<std::size_t>(exterior_factor)]
                                 [static_cast<std::size_t>(symmetric_factor)];
      if (scalar == 0) {
        continue;
      }
      const MinusState symmetric_state =
          minus[static_cast<std::size_t>(symmetric_factor)]
               [static_cast<std::size_t>(source.minus[
                   static_cast<std::size_t>(symmetric_factor)])];
      if ((exterior_state.exterior & 2) != 0 && symmetric_state.x > 0) {
        MinusState next_exterior = exterior_state;
        MinusState next_symmetric = symmetric_state;
        const int contraction_sign = (next_exterior.exterior & 1) != 0
            ? -1
            : 1;
        next_exterior.exterior &= ~2;
        --next_symmetric.x;
        ++next_symmetric.y_minus;
        State target = source;
        if (exterior_factor == symmetric_factor) {
          MinusState next = exterior_state;
          next.exterior &= ~2;
          --next.x;
          ++next.y_minus;
          target.minus[static_cast<std::size_t>(exterior_factor)] =
              minus_index[static_cast<std::size_t>(exterior_factor)].at(next);
        } else {
          target.minus[static_cast<std::size_t>(exterior_factor)] =
              minus_index[static_cast<std::size_t>(exterior_factor)]
                  .at(next_exterior);
          target.minus[static_cast<std::size_t>(symmetric_factor)] =
              minus_index[static_cast<std::size_t>(symmetric_factor)]
                  .at(next_symmetric);
        }
        append(result, target, tensor_sign * contraction_sign * scalar
                                   * symmetric_state.x);
      }
      if ((exterior_state.exterior & 1) != 0 && symmetric_state.x > 0) {
        MinusState next_exterior = exterior_state;
        MinusState next_symmetric = symmetric_state;
        next_exterior.exterior &= ~1;
        --next_symmetric.x;
        ++next_symmetric.y_plus;
        State target = source;
        if (exterior_factor == symmetric_factor) {
          MinusState next = exterior_state;
          next.exterior &= ~1;
          --next.x;
          ++next.y_plus;
          target.minus[static_cast<std::size_t>(exterior_factor)] =
              minus_index[static_cast<std::size_t>(exterior_factor)].at(next);
        } else {
          target.minus[static_cast<std::size_t>(exterior_factor)] =
              minus_index[static_cast<std::size_t>(exterior_factor)]
                  .at(next_exterior);
          target.minus[static_cast<std::size_t>(symmetric_factor)] =
              minus_index[static_cast<std::size_t>(symmetric_factor)]
                  .at(next_symmetric);
        }
        append(result, target, tensor_sign * scalar * symmetric_state.x);
      }
    }
  }
  return result;
}

Matrix make_operator(
    const std::vector<State>& source, const std::vector<State>& target,
    const std::map<State, int>& target_index,
    const std::vector<std::map<State, std::int64_t>>& images) {
  Matrix result(target.size(),
                std::vector<std::int64_t>(source.size(), 0));
  for (std::size_t column = 0; column < source.size(); ++column) {
    for (const auto& [image, coefficient] : images[column]) {
      const auto found = target_index.find(image);
      if (found == target_index.end()) {
        throw std::runtime_error("operator image has wrong weight or parity");
      }
      result[static_cast<std::size_t>(found->second)][column] = coefficient;
    }
  }
  return result;
}

int rank(Matrix matrix) {
  if (matrix.empty()) {
    return 0;
  }
  const std::size_t columns = matrix.front().size();
  std::size_t pivot_row = 0;
  for (std::size_t column = 0; column < columns
                                && pivot_row < matrix.size(); ++column) {
    std::size_t pivot = pivot_row;
    while (pivot < matrix.size() && matrix[pivot][column] == 0) {
      ++pivot;
    }
    if (pivot == matrix.size()) {
      continue;
    }
    std::swap(matrix[pivot_row], matrix[pivot]);
    const std::int64_t scale = inverse(matrix[pivot_row][column]);
    for (std::size_t entry = column; entry < columns; ++entry) {
      matrix[pivot_row][entry] = normalize(matrix[pivot_row][entry] * scale);
    }
    for (std::size_t row = 0; row < matrix.size(); ++row) {
      if (row == pivot_row || matrix[row][column] == 0) {
        continue;
      }
      const std::int64_t factor = matrix[row][column];
      for (std::size_t entry = column; entry < columns; ++entry) {
        matrix[row][entry] = normalize(matrix[row][entry]
                                       - factor * matrix[pivot_row][entry]);
      }
    }
    ++pivot_row;
  }
  return static_cast<int>(pivot_row);
}

Matrix nullspace(Matrix matrix) {
  const std::size_t rows = matrix.size();
  const std::size_t columns = rows == 0 ? 0 : matrix.front().size();
  std::vector<int> pivot_column;
  std::size_t pivot_row = 0;
  for (std::size_t column = 0; column < columns && pivot_row < rows;
       ++column) {
    std::size_t pivot = pivot_row;
    while (pivot < rows && matrix[pivot][column] == 0) {
      ++pivot;
    }
    if (pivot == rows) {
      continue;
    }
    std::swap(matrix[pivot_row], matrix[pivot]);
    const std::int64_t scale = inverse(matrix[pivot_row][column]);
    for (std::size_t entry = column; entry < columns; ++entry) {
      matrix[pivot_row][entry] = normalize(matrix[pivot_row][entry] * scale);
    }
    for (std::size_t row = 0; row < rows; ++row) {
      if (row == pivot_row || matrix[row][column] == 0) {
        continue;
      }
      const std::int64_t factor = matrix[row][column];
      for (std::size_t entry = column; entry < columns; ++entry) {
        matrix[row][entry] = normalize(matrix[row][entry]
                                       - factor * matrix[pivot_row][entry]);
      }
    }
    pivot_column.push_back(static_cast<int>(column));
    ++pivot_row;
  }
  std::vector<bool> is_pivot(columns, false);
  for (const int column : pivot_column) {
    is_pivot[static_cast<std::size_t>(column)] = true;
  }
  Matrix result(columns);
  for (std::size_t free_column = 0; free_column < columns; ++free_column) {
    if (is_pivot[free_column]) {
      continue;
    }
    std::vector<std::int64_t> vector(columns, 0);
    vector[free_column] = 1;
    for (std::size_t row = 0; row < pivot_column.size(); ++row) {
      vector[static_cast<std::size_t>(pivot_column[row])] =
          normalize(-matrix[row][free_column]);
    }
    for (std::size_t entry = 0; entry < columns; ++entry) {
      result[entry].push_back(vector[entry]);
    }
  }
  return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
  if (left.empty() || right.empty()) {
    return Matrix(left.size(),
                  std::vector<std::int64_t>(right.empty() ? 0U
                                                           : right.front().size(),
                                            0));
  }
  if (left.front().size() != right.size()) {
    throw std::runtime_error("incompatible matrices");
  }
  Matrix result(left.size(),
                std::vector<std::int64_t>(right.front().size(), 0));
  for (std::size_t row = 0; row < left.size(); ++row) {
    for (std::size_t middle = 0; middle < right.size(); ++middle) {
      if (left[row][middle] == 0) {
        continue;
      }
      for (std::size_t column = 0; column < right.front().size(); ++column) {
        result[row][column] = normalize(result[row][column]
                                         + left[row][middle]
                                               * right[middle][column]);
      }
    }
  }
  return result;
}

bool zero(const Matrix& matrix) {
  for (const auto& row : matrix) {
    for (const std::int64_t entry : row) {
      if (entry != 0) {
        return false;
      }
    }
  }
  return true;
}

struct WeightSpaces {
  std::array<std::vector<State>, 2> zero;
  std::array<std::vector<State>, 2> two;
  std::array<std::map<State, int>, 2> zero_index;
  std::array<std::map<State, int>, 2> two_index;
};

WeightSpaces make_weight_spaces(
    const std::array<std::vector<MinusState>, 2>& minus,
    const std::array<std::vector<PlusState>, 3>& plus) {
  WeightSpaces result;
  for (std::size_t first = 0; first < minus[0].size(); ++first) {
    for (std::size_t second = 0; second < minus[1].size(); ++second) {
      for (std::size_t first_plus = 0; first_plus < plus[0].size();
           ++first_plus) {
        for (std::size_t second_plus = 0; second_plus < plus[1].size();
             ++second_plus) {
          for (std::size_t third_plus = 0; third_plus < plus[2].size();
               ++third_plus) {
            const State state{{static_cast<int>(first), static_cast<int>(second)},
                              {static_cast<int>(first_plus),
                               static_cast<int>(second_plus),
                               static_cast<int>(third_plus)}};
            if (total_degree(state, minus, plus) != kOuterDegree) {
              continue;
            }
            const int state_parity = parity(state, minus);
            const int state_weight = total_weight(state, minus, plus);
            if (state_weight == 0) {
              result.zero[static_cast<std::size_t>(state_parity)].push_back(state);
            } else if (state_weight == 2) {
              result.two[static_cast<std::size_t>(state_parity)].push_back(state);
            }
          }
        }
      }
    }
  }
  for (int state_parity = 0; state_parity < 2; ++state_parity) {
    for (std::size_t index = 0;
         index < result.zero[static_cast<std::size_t>(state_parity)].size();
         ++index) {
      result.zero_index[static_cast<std::size_t>(state_parity)].emplace(
          result.zero[static_cast<std::size_t>(state_parity)][index],
          static_cast<int>(index));
    }
    for (std::size_t index = 0;
         index < result.two[static_cast<std::size_t>(state_parity)].size();
         ++index) {
      result.two_index[static_cast<std::size_t>(state_parity)].emplace(
          result.two[static_cast<std::size_t>(state_parity)][index],
          static_cast<int>(index));
    }
  }
  return result;
}

void test_case(const std::string& name,
               const std::array<std::array<int, 2>, 2>& coupling,
               const std::array<std::vector<MinusState>, 2>& minus,
               const std::array<std::vector<PlusState>, 3>& plus,
               const std::array<std::map<MinusState, int>, 2>& minus_index,
               const std::array<std::map<PlusState, int>, 3>& plus_index,
               const WeightSpaces& spaces) {
  std::array<Matrix, 2> raising;
  std::array<Matrix, 2> differential_zero;
  std::array<Matrix, 2> differential_two;
  for (int state_parity = 0; state_parity < 2; ++state_parity) {
    const auto& zero_space = spaces.zero[static_cast<std::size_t>(state_parity)];
    const auto& two_space = spaces.two[static_cast<std::size_t>(state_parity)];
    std::vector<std::map<State, std::int64_t>> raising_images;
    std::vector<std::map<State, std::int64_t>> differential_images;
    raising_images.reserve(zero_space.size());
    differential_images.reserve(zero_space.size());
    for (const State& state : zero_space) {
      raising_images.push_back(
          apply_raising(state, minus, plus, minus_index, plus_index));
      differential_images.push_back(
          apply_differential(state, coupling, minus, minus_index));
    }
    raising[static_cast<std::size_t>(state_parity)] = make_operator(
        zero_space, two_space,
        spaces.two_index[static_cast<std::size_t>(state_parity)],
        raising_images);
    differential_zero[static_cast<std::size_t>(state_parity)] = make_operator(
        zero_space, spaces.zero[static_cast<std::size_t>(1 - state_parity)],
        spaces.zero_index[static_cast<std::size_t>(1 - state_parity)],
        differential_images);

    std::vector<std::map<State, std::int64_t>> differential_two_images;
    differential_two_images.reserve(two_space.size());
    for (const State& state : two_space) {
      differential_two_images.push_back(
          apply_differential(state, coupling, minus, minus_index));
    }
    differential_two[static_cast<std::size_t>(state_parity)] = make_operator(
        two_space, spaces.two[static_cast<std::size_t>(1 - state_parity)],
        spaces.two_index[static_cast<std::size_t>(1 - state_parity)],
        differential_two_images);
  }

  const Matrix kernel_even = nullspace(raising[0]);
  const Matrix kernel_odd = nullspace(raising[1]);
  const int invariant_even = static_cast<int>(kernel_even.empty()
      ? 0U
      : kernel_even.front().size());
  const int invariant_odd = static_cast<int>(kernel_odd.empty()
      ? 0U
      : kernel_odd.front().size());
  const int rank_even_to_odd = rank(multiply(differential_zero[0], kernel_even));
  const int rank_odd_to_even = rank(multiply(differential_zero[1], kernel_odd));

  const bool square_zero = zero(multiply(differential_zero[1],
                                         differential_zero[0]))
      && zero(multiply(differential_zero[0], differential_zero[1]));
  const bool equivariant = multiply(raising[1], differential_zero[0])
      == multiply(differential_two[0], raising[0])
      && multiply(raising[0], differential_zero[1])
      == multiply(differential_two[1], raising[1]);
  const int homology_odd = invariant_odd - rank_even_to_odd - rank_odd_to_even;
  const int homology_even = invariant_even - rank_even_to_odd - rank_odd_to_even;

  std::cout << "CROSS_KOSZUL"
            << " case=" << name
            << " degree=" << kOuterDegree
            << " weight0_even=" << spaces.zero[0].size()
            << " weight0_odd=" << spaces.zero[1].size()
            << " invariant_even=" << invariant_even
            << " invariant_odd=" << invariant_odd
            << " rank_even_to_odd=" << rank_even_to_odd
            << " rank_odd_to_even=" << rank_odd_to_even
            << " homology_even=" << homology_even
            << " homology_odd=" << homology_odd
            << " square_zero=" << (square_zero ? "true" : "false")
            << " equivariant=" << (equivariant ? "true" : "false")
            << " result="
            << (square_zero && equivariant ? "PASS" : "FAIL") << '\n';
}

}  // namespace

int main() {
  const std::array<std::vector<MinusState>, 2> minus{
      make_minus_states(3), make_minus_states(3)};
  const std::array<std::vector<PlusState>, 3> plus{
      make_plus_states(1), make_plus_states(2), make_plus_states(2)};
  const std::array<std::map<MinusState, int>, 2> minus_index{
      index_minus(minus[0]), index_minus(minus[1])};
  const std::array<std::map<PlusState, int>, 3> plus_index{
      index_plus(plus[0]), index_plus(plus[1]), index_plus(plus[2])};
  const WeightSpaces spaces = make_weight_spaces(minus, plus);
  test_case("factorwise", {{{1, 0}, {0, 1}}}, minus, plus, minus_index,
            plus_index, spaces);
  test_case("generic_full", {{{1, 1}, {1, 2}}}, minus, plus, minus_index,
            plus_index, spaces);
  test_case("rank_one", {{{1, 1}, {1, 1}}}, minus, plus, minus_index,
            plus_index, spaces);
  return 0;
}
