#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

using Polynomial = std::unordered_map<std::uint64_t, cpp_int>;

std::uint64_t pair_key(int left, int right) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(left)) << 32U)
        | static_cast<std::uint32_t>(right);
}

std::pair<int, int> decode_pair(std::uint64_t key) {
    return {
        static_cast<int>(static_cast<std::uint32_t>(key >> 32U)),
        static_cast<int>(static_cast<std::uint32_t>(key))
    };
}

void add_scaled(Polynomial& target, const Polynomial& source, int sign) {
    for (const auto& [key, coefficient] : source) {
        target[key] += sign * coefficient;
        if (target[key] == 0) {
            target.erase(key);
        }
    }
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
    Polynomial output;
    for (const auto& [left_key, left_coefficient] : left) {
        const auto [a, b] = decode_pair(left_key);
        for (const auto& [right_key, right_coefficient] : right) {
            const auto [c, d] = decode_pair(right_key);
            for (int e = std::abs(a - c); e <= a + c; e += 2) {
                for (int f = std::abs(b - d); f <= b + d; f += 2) {
                    output[pair_key(e, f)] += left_coefficient * right_coefficient;
                }
            }
        }
    }
    return output;
}

Polynomial complete_character(int degree) {
    Polynomial output;
    if (degree < 0) {
        return output;
    }
    for (int left = 0; left <= degree; ++left) {
        output.emplace(pair_key(left, degree - left), 1);
    }
    return output;
}

int permutation_sign(const std::array<int, 4>& permutation, int length) {
    int inversions = 0;
    for (int i = 0; i < length; ++i) {
        for (int j = i + 1; j < length; ++j) {
            if (permutation[static_cast<std::size_t>(i)]
                > permutation[static_cast<std::size_t>(j)]) {
                ++inversions;
            }
        }
    }
    return (inversions & 1) == 0 ? 1 : -1;
}

Polynomial schur_character(const std::vector<int>& partition) {
    const int length = static_cast<int>(partition.size());
    if (length == 0) {
        return {{pair_key(0, 0), 1}};
    }
    if (length > 4) {
        throw std::runtime_error("only partitions of length at most four are supported");
    }

    std::array<int, 4> permutation{0, 1, 2, 3};
    Polynomial determinant;
    do {
        Polynomial term{{pair_key(0, 0), 1}};
        bool zero = false;
        for (int row = 0; row < length; ++row) {
            const int column = permutation[static_cast<std::size_t>(row)];
            const int degree = partition[static_cast<std::size_t>(row)] - row + column;
            if (degree < 0) {
                zero = true;
                break;
            }
            term = multiply(term, complete_character(degree));
        }
        if (!zero) {
            add_scaled(determinant, term, permutation_sign(permutation, length));
        }
    } while (std::next_permutation(
        permutation.begin(), permutation.begin() + length
    ));
    return determinant;
}

Polynomial difference_power(int exponent) {
    Polynomial states{{pair_key(0, 0), 1}};
    for (int step = 0; step < exponent; ++step) {
        Polynomial following;
        for (const auto& [key, coefficient] : states) {
            const auto [left, right] = decode_pair(key);
            for (int output = std::abs(left - 1); output <= left + 1; output += 2) {
                following[pair_key(output, right)] += coefficient;
            }
            for (int output = std::abs(right - 1); output <= right + 1; output += 2) {
                following[pair_key(left, output)] -= coefficient;
            }
        }
        states = std::move(following);
    }
    return states;
}

cpp_int constant_coefficient(
    const Polynomial& left,
    const Polynomial& right
) {
    cpp_int answer = 0;
    for (const auto& [key, coefficient] : left) {
        const auto found = right.find(key);
        if (found != right.end()) {
            answer += coefficient * found->second;
        }
    }
    return answer;
}

void partitions_rec(
    int remaining,
    int maximum_part,
    std::vector<int>& current,
    std::vector<std::vector<int>>& output
) {
    if (remaining == 0) {
        output.push_back(current);
        return;
    }
    if (current.size() == 4U) {
        return;
    }
    for (int part = std::min(remaining, maximum_part); part >= 1; --part) {
        current.push_back(part);
        partitions_rec(remaining - part, part, current, output);
        current.pop_back();
    }
}

std::vector<std::vector<int>> partitions_through(int maximum_degree) {
    std::vector<std::vector<int>> output{{}};
    for (int degree = 1; degree <= maximum_degree; ++degree) {
        std::vector<int> current;
        partitions_rec(degree, degree, current, output);
    }
    return output;
}

void print_partition(const std::vector<int>& partition) {
    std::cout << '[';
    for (std::size_t index = 0; index < partition.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << partition[index];
    }
    std::cout << ']';
}

int parse_positive(const char* text, const char* name) {
    const long value = std::strtol(text, nullptr, 10);
    if (value <= 0 || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int maximum_degree = 20;
        int maximum_half_minus_count = 8;
        if (argc != 1 && argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_schur_kernel [maximum_degree maximum_half_minus_count]"
            );
        }
        if (argc == 3) {
            maximum_degree = parse_positive(argv[1], "maximum degree");
            maximum_half_minus_count = parse_positive(
                argv[2], "maximum half-minus count"
            );
        }

        const auto partitions = partitions_through(maximum_degree);
        bool negative = false;
        #pragma omp parallel for schedule(dynamic)
        for (std::int64_t raw_index = 0;
             raw_index < static_cast<std::int64_t>(partitions.size());
             ++raw_index) {
            const auto& partition = partitions[static_cast<std::size_t>(raw_index)];
            const Polynomial schur = schur_character(partition);
            for (int half = 1; half <= maximum_half_minus_count; ++half) {
                const cpp_int value = constant_coefficient(
                    difference_power(2 * half), schur
                );
                if (value < 0) {
                    #pragma omp critical
                    {
                        if (!negative) {
                            negative = true;
                            std::cout << "SU2_SCHUR_KERNEL negative m=" << half
                                      << " partition=";
                            print_partition(partition);
                            std::cout << " value=" << value << '\n';
                        }
                    }
                    break;
                }
            }
        }
        std::cout << "SU2_SCHUR_KERNEL SEARCH: "
                  << (negative
                      ? "NEGATIVE WITNESS FOUND"
                      : "NO NEGATIVE IN TESTED DOMAIN")
                  << " partitions=" << partitions.size()
                  << " maximum_degree=" << maximum_degree
                  << " maximum_m=" << maximum_half_minus_count << '\n';
        return negative ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SCHUR_KERNEL FAILURE: " << error.what() << '\n';
        return 1;
    }
}
