#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Matrix = std::vector<std::vector<cpp_int>>;

int parse_nonnegative(const char* text, const char* name) {
    const std::string value{text};
    std::size_t used = 0U;
    const long parsed = std::stol(value, &used, 10);
    if (used != value.size() || parsed < 0L) {
        throw std::invalid_argument(std::string{name} + " must be nonnegative");
    }
    return static_cast<int>(parsed);
}

Matrix zero_matrix(const int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
}

template <typename Function>
void for_each_output(
    const int level, const int left, const int right, Function function
) {
    const int upper = std::min(left + right, 2 * level - left - right);
    for (int target = std::abs(left - right); target <= upper; target += 2) {
        function(target);
    }
}

Matrix additive_update(const int level, const Matrix& input, const int label) {
    Matrix output = zero_matrix(level + 1);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const cpp_int& value =
                input[static_cast<std::size_t>(left)][static_cast<std::size_t>(right)];
            if (value == 0) {
                continue;
            }
            for_each_output(level, left, label, [&](const int target) {
                output[static_cast<std::size_t>(target)]
                      [static_cast<std::size_t>(right)] += value;
            });
            for_each_output(level, right, label, [&](const int target) {
                output[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(target)] += value;
            });
        }
    }
    return output;
}

Matrix first_update(const int level, const Matrix& input, const int label) {
    Matrix output = zero_matrix(level + 1);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const cpp_int& value =
                input[static_cast<std::size_t>(left)][static_cast<std::size_t>(right)];
            if (value == 0) {
                continue;
            }
            for_each_output(level, left, label, [&](const int target) {
                output[static_cast<std::size_t>(target)]
                      [static_cast<std::size_t>(right)] += value;
            });
        }
    }
    return output;
}

Matrix second_update(const int level, const Matrix& input, const int label) {
    Matrix output = zero_matrix(level + 1);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const cpp_int& value =
                input[static_cast<std::size_t>(left)][static_cast<std::size_t>(right)];
            if (value == 0) {
                continue;
            }
            for_each_output(level, right, label, [&](const int target) {
                output[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(target)] += value;
            });
        }
    }
    return output;
}

Matrix carrier(
    const int level, const Matrix& prefix, const std::vector<int>& suffix
) {
    Matrix first = prefix;
    Matrix second = prefix;
    for (const int label : suffix) {
        first = first_update(level, first, label);
        second = second_update(level, second, label);
    }
    Matrix result = zero_matrix(level + 1);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            result[static_cast<std::size_t>(left)][static_cast<std::size_t>(right)]
                = first[static_cast<std::size_t>(left)]
                       [static_cast<std::size_t>(right)]
                + second[static_cast<std::size_t>(left)]
                         [static_cast<std::size_t>(right)];
        }
    }
    return result;
}

void render(const std::vector<int>& labels) {
    std::cout << '[';
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << labels[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::invalid_argument(
                "usage: verify_su2_level_one_odd_prefix_carrier "
                "MAXIMUM_LEVEL MAXIMUM_EVEN_FACTORS [ODD_COUNT]"
            );
        }
        const int maximum_level = parse_nonnegative(argv[1], "MAXIMUM_LEVEL");
        const int maximum_even_factors =
            parse_nonnegative(argv[2], "MAXIMUM_EVEN_FACTORS");
        const int odd_count = argc == 4
            ? parse_nonnegative(argv[3], "ODD_COUNT")
            : 1;
        if (maximum_level < 2) {
            throw std::invalid_argument("MAXIMUM_LEVEL must be at least two");
        }
        if (odd_count != 1 && odd_count != 2) {
            throw std::invalid_argument("ODD_COUNT must be one or two");
        }

        std::uint64_t cuts = 0U;
        std::uint64_t coefficients = 0U;
        for (int level = 2; level <= maximum_level; ++level) {
            std::vector<int> prefix_even;
            const auto visit_prefix = [&](const auto& self, const int next_even,
                                          const int remaining) -> void {
                std::vector<int> suffix_even;
                const auto visit_suffix = [&](
                    const auto& recurse, const int suffix_next,
                    const int suffix_remaining
                ) -> void {
                    Matrix base = zero_matrix(level + 1);
                    base[1U][0U] = 1;
                    base[0U][1U] = -1;
                    for (const int label : prefix_even) {
                        base = additive_update(level, base, label);
                    }
                    for (int first_odd = 1; first_odd <= level;
                         first_odd += 2) {
                        const int last_odd = odd_count == 1 ? first_odd : level;
                        for (int second_odd = first_odd; second_odd <= last_odd;
                             second_odd += 2) {
                            const std::vector<int> odds = odd_count == 1
                                ? std::vector<int>{first_odd}
                                : std::vector<int>{first_odd, second_odd};
                            const unsigned int assignments =
                                UINT32_C(1) << static_cast<unsigned int>(odd_count);
                            for (unsigned int mask = 0U; mask < assignments;
                                 ++mask) {
                                Matrix prefix = base;
                                std::vector<int> suffix = suffix_even;
                                for (int index = 0; index < odd_count; ++index) {
                                    const int odd =
                                        odds[static_cast<std::size_t>(index)];
                                    const unsigned int bit =
                                        UINT32_C(1) << static_cast<unsigned int>(index);
                                    if ((mask & bit) != 0U) {
                                        prefix = additive_update(level, prefix, odd);
                                    } else {
                                        suffix.insert(suffix.begin(), odd);
                                    }
                                }
                                const Matrix value = carrier(level, prefix, suffix);
                                ++cuts;
                                for (int target = 1; target <= level; ++target) {
                                    ++coefficients;
                                    const cpp_int& coefficient =
                                        value[static_cast<std::size_t>(target)][0U];
                                    if (coefficient < 0) {
                                        std::cout << "FAIL level=" << level
                                                  << " odds=";
                                        render(odds);
                                        std::cout << " odd_prefix_mask=" << mask
                                                  << " prefix_even=";
                                        render(prefix_even);
                                        std::cout << " suffix_even=";
                                        render(suffix_even);
                                        std::cout << " target=" << target
                                                  << " coefficient=" << coefficient
                                                  << '\n';
                                        std::exit(EXIT_FAILURE);
                                    }
                                }
                            }
                        }
                    }
                    if (suffix_remaining == 0) {
                        return;
                    }
                    for (int label = suffix_next; label <= level; label += 2) {
                        suffix_even.push_back(label);
                        recurse(recurse, label, suffix_remaining - 1);
                        suffix_even.pop_back();
                    }
                };
                visit_suffix(visit_suffix, 2, remaining);
                if (remaining == 0) {
                    return;
                }
                for (int label = next_even; label <= level; label += 2) {
                    prefix_even.push_back(label);
                    self(self, label, remaining - 1);
                    prefix_even.pop_back();
                }
            };
            visit_prefix(visit_prefix, 2, maximum_even_factors);
        }
        std::cout << "SU2_LEVEL_ODD_PREFIX_CARRIER result=PASS"
                  << " maximum_level=" << maximum_level
                  << " maximum_even_factors=" << maximum_even_factors
                  << " odd_count=" << odd_count
                  << " cuts=" << cuts
                  << " coefficients=" << coefficients << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
