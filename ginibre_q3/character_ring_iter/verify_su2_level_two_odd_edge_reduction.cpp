#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Character = std::vector<cpp_int>;
using Matrix = std::vector<Character>;

Matrix zero_matrix(int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        Character(static_cast<std::size_t>(dimension))
    );
}

template <class Function>
void for_each_output(int level, int left, int right, Function function) {
    const int maximum = std::min(left + right, 2 * level - left - right);
    for (int output = std::abs(left - right);
         output <= maximum;
         output += 2) {
        function(output);
    }
}

Matrix update(int level, const Matrix& current, int label) {
    Matrix result = zero_matrix(level + 1);
    for (int first = 0; first <= level; ++first) {
        for (int second = 0; second <= level; ++second) {
            const cpp_int value = current[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
            if (value == 0) {
                continue;
            }
            for_each_output(level, first, label, [&](const int output) {
                result[static_cast<std::size_t>(output)]
                    [static_cast<std::size_t>(second)] += value;
            });
            for_each_output(level, second, label, [&](const int output) {
                result[static_cast<std::size_t>(first)]
                    [static_cast<std::size_t>(output)] += value;
            });
        }
    }
    return result;
}

Character defect_column(int level, const Matrix& coefficients, int second) {
    Character result(static_cast<std::size_t>(level + 1));
    for (int first = 0; first <= level; ++first) {
        for_each_output(level, first, second, [&](const int output) {
            result[static_cast<std::size_t>(first)]
                += coefficients[static_cast<std::size_t>(output)][0];
        });
        result[static_cast<std::size_t>(first)]
            -= coefficients[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
    }
    return result;
}

Character multiply(int level, const Character& input, int label) {
    Character result(static_cast<std::size_t>(level + 1));
    for (int first = 0; first <= level; ++first) {
        const cpp_int value = input[static_cast<std::size_t>(first)];
        if (value == 0) {
            continue;
        }
        for_each_output(level, first, label, [&](const int output) {
            result[static_cast<std::size_t>(output)] += value;
        });
    }
    return result;
}

void add(Character& target, const Character& source, int sign) {
    for (std::size_t index = 0; index < target.size(); ++index) {
        target[index] += sign * source[index];
    }
}

void print_word(const std::vector<int>& word) {
    std::cout << '[';
    for (std::size_t index = 0; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_level_two_odd_edge_reduction "
                "MAXIMUM_LEVEL MAXIMUM_EVEN_FACTORS"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_level < 2 || maximum_factors < 0) {
            throw std::runtime_error("invalid bound");
        }
        long long words = 0;
        long long identities = 0;
        bool first_reservoir_suffices = true;
        bool second_reservoir_suffices = true;
        bool f2_nonnegative = true;
        long long f2_entries = 0;
        for (int level = 2; level <= maximum_level; ++level) {
            std::vector<int> word;
            const auto visit = [&](const auto& self, int next_even) -> void {
                ++words;
                Matrix coefficients = zero_matrix(level + 1);
                coefficients[0][0] = 1;
                Matrix exterior = zero_matrix(level + 1);
                exterior[1][0] = 1;
                exterior[0][1] = -1;
                for (const int label : word) {
                    coefficients = update(level, coefficients, label);
                    exterior = update(level, exterior, label);
                }
                std::vector<Character> defects(
                    static_cast<std::size_t>(level + 1)
                );
                for (int label = 0; label <= level; ++label) {
                    defects[static_cast<std::size_t>(label)]
                        = defect_column(level, coefficients, label);
                }
                if (level >= 2) {
                    const Character& f2 = defects[2];
                    for (int target = 0; target <= level; ++target) {
                        ++f2_entries;
                        const cpp_int value
                            = f2[static_cast<std::size_t>(target)];
                        if (value < 0 && f2_nonnegative) {
                            f2_nonnegative = false;
                            std::cout << "F2_NEGATIVE level=" << level
                                      << " even_word=";
                            print_word(word);
                            std::cout << " target=" << target
                                      << " value=" << value << '\n';
                        }
                    }
                }
                for (int q = 1; q <= level; q += 2) {
                    for (int r = q; r <= level; r += 2) {
                        Matrix extended = update(level, exterior, q);
                        extended = update(level, extended, r);
                        Character first(static_cast<std::size_t>(level + 1));
                        for_each_output(level, r, 1, [&](const int neighbor) {
                            add(first, multiply(
                                level,
                                defects[static_cast<std::size_t>(neighbor)],
                                q
                            ), 1);
                        });
                        Character second(static_cast<std::size_t>(level + 1));
                        for_each_output(level, q, 1, [&](const int neighbor) {
                            add(second, multiply(
                                level,
                                defects[static_cast<std::size_t>(neighbor)],
                                r
                            ), 1);
                        });
                        Character negative(static_cast<std::size_t>(level + 1));
                        for_each_output(level, q, r, [&](const int output) {
                            add(
                                negative,
                                defects[static_cast<std::size_t>(output)],
                                1
                            );
                        });
                        negative = multiply(level, negative, 1);
                        for (int target = 1; target <= level; ++target) {
                            ++identities;
                            const std::size_t index
                                = static_cast<std::size_t>(target);
                            const cpp_int right
                                = first[index] + second[index] - negative[index];
                            const cpp_int left = extended[index][0];
                            if (left != right) {
                                std::cout << "FAIL identity level=" << level
                                          << " even_word=";
                                print_word(word);
                                std::cout << " q=" << q << " r=" << r
                                          << " target=" << target
                                          << " left=" << left
                                          << " right=" << right << '\n';
                                std::exit(EXIT_FAILURE);
                            }
                            if (left < 0) {
                                std::cout << "FAIL positivity level=" << level
                                          << " even_word=";
                                print_word(word);
                                std::cout << " q=" << q << " r=" << r
                                          << " target=" << target
                                          << " value=" << left << '\n';
                                std::exit(EXIT_FAILURE);
                            }
                            if (first[index] - negative[index] < 0
                                && first_reservoir_suffices) {
                                first_reservoir_suffices = false;
                                std::cout << "FIRST_RESERVOIR_FAIL level="
                                          << level << " even_word=";
                                print_word(word);
                                std::cout << " q=" << q << " r=" << r
                                          << " target=" << target
                                          << " value="
                                          << first[index] - negative[index]
                                          << '\n';
                            }
                            if (second[index] - negative[index] < 0
                                && second_reservoir_suffices) {
                                second_reservoir_suffices = false;
                                std::cout << "SECOND_RESERVOIR_FAIL level="
                                          << level << " even_word=";
                                print_word(word);
                                std::cout << " q=" << q << " r=" << r
                                          << " target=" << target
                                          << " value="
                                          << second[index] - negative[index]
                                          << '\n';
                            }
                        }
                    }
                }
                if (word.size()
                    == static_cast<std::size_t>(maximum_factors)) {
                    return;
                }
                for (int label = next_even; label <= level; label += 2) {
                    word.push_back(label);
                    self(self, label);
                    word.pop_back();
                }
            };
            visit(visit, 2);
        }
        std::cout << "SU2_LEVEL_TWO_ODD_EDGE_REDUCTION PASS words=" << words
                  << " identities=" << identities
                  << " first_reservoir_suffices="
                  << (first_reservoir_suffices ? "true" : "false")
                  << " second_reservoir_suffices="
                  << (second_reservoir_suffices ? "true" : "false")
                  << " f2_nonnegative="
                  << (f2_nonnegative ? "true" : "false")
                  << " f2_entries=" << f2_entries
                  << " maximum_level=" << maximum_level
                  << " maximum_even_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
