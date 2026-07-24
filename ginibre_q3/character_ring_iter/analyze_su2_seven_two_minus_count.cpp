#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

cpp_int invariant_multiplicity(const std::vector<int>& labels) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    std::vector<cpp_int> current(static_cast<std::size_t>(total + 1));
    std::vector<cpp_int> next(static_cast<std::size_t>(total + 1));
    current[0U] = 1;
    int maximum = 0;
    for (int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int source = 0; source <= maximum; ++source) {
            const cpp_int& multiplicity =
                current[static_cast<std::size_t>(source)];
            if (multiplicity == 0) {
                continue;
            }
            for (int target = std::abs(source - label);
                 target <= source + label; target += 2) {
                next[static_cast<std::size_t>(target)] += multiplicity;
            }
        }
        maximum += label;
        current.swap(next);
    }
    return current[0U];
}

bool disjoint_support(
    const std::array<int, 2>& minus,
    const std::array<int, 5>& plus
) {
    for (int first : minus) {
        for (int second : plus) {
            if (first == second) {
                return false;
            }
        }
    }
    return true;
}

std::array<cpp_int, 2> negative_middle_by_orientation(
    const std::array<int, 2>& minus,
    const std::array<int, 5>& plus
) {
    std::array<cpp_int, 2> result{};
    for (std::size_t selected_minus = 0U;
         selected_minus < minus.size(); ++selected_minus) {
        const std::size_t other_minus = 1U - selected_minus;
        for (std::size_t first = 0U; first < plus.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < plus.size(); ++second) {
                const std::vector<int> triple{
                    minus[selected_minus], plus[first], plus[second]
                };
                const cpp_int triple_multiplicity =
                    invariant_multiplicity(triple);
                if (triple_multiplicity == 0) {
                    continue;
                }
                std::vector<int> complement{minus[other_minus]};
                for (std::size_t index = 0U;
                     index < plus.size(); ++index) {
                    if (index != first && index != second) {
                        complement.push_back(plus[index]);
                    }
                }
                result[selected_minus] += triple_multiplicity
                    * invariant_multiplicity(complement);
            }
        }
    }
    return result;
}

void print_case(
    const std::array<int, 2>& minus,
    const std::array<int, 5>& plus,
    const cpp_int& invariant,
    const cpp_int& negative
) {
    std::cout << " minus=[" << minus[0] << ',' << minus[1] << "] plus=[";
    for (std::size_t index = 0U; index < plus.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << plus[index];
    }
    std::cout << "] N7=" << invariant << " Tminus=" << negative
              << " margin=" << invariant - negative;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_seven_two_minus_count MAXIMUM_LABEL "
                "[--orientation]"
            );
        }
        const bool orientation_mode =
            argc == 3 && std::string(argv[2]) == "--orientation";
        if (argc == 3 && !orientation_mode) {
            throw std::runtime_error("unknown mode");
        }
        const long parsed = std::strtol(argv[1], nullptr, 10);
        if (parsed < 1 || parsed > std::numeric_limits<int>::max()) {
            throw std::runtime_error("invalid maximum label");
        }
        const int maximum_label = static_cast<int>(parsed);
        std::size_t cases = 0U;
        bool found = false;
        cpp_int minimum_margin = 0;
        bool have_margin = false;
        std::array<int, 2> minimum_minus{};
        std::array<int, 5> minimum_plus{};
        cpp_int minimum_invariant = 0;
        cpp_int minimum_negative = 0;

        for (int first_minus = 1;
             first_minus <= maximum_label && !found; ++first_minus) {
            for (int second_minus = first_minus;
                 second_minus <= maximum_label && !found; ++second_minus) {
                const std::array<int, 2> minus{
                    first_minus, second_minus
                };
                for (int a = 1; a <= maximum_label && !found; ++a) {
                    for (int b = a; b <= maximum_label && !found; ++b) {
                        for (int c = b; c <= maximum_label && !found; ++c) {
                            for (int d = c;
                                 d <= maximum_label && !found; ++d) {
                                for (int e = d;
                                     e <= maximum_label; ++e) {
                                    const std::array<int, 5> plus{
                                        a, b, c, d, e
                                    };
                                    if (!disjoint_support(minus, plus)) {
                                        continue;
                                    }
                                    std::vector<int> labels{
                                        first_minus, second_minus,
                                        a, b, c, d, e
                                    };
                                    const cpp_int invariant =
                                        invariant_multiplicity(labels);
                                    if (invariant == 0) {
                                        continue;
                                    }
                                    const auto by_orientation =
                                        negative_middle_by_orientation(
                                            minus, plus
                                        );
                                    const cpp_int negative =
                                        by_orientation[0U]
                                        + by_orientation[1U];
                                    const cpp_int margin =
                                        invariant - (
                                            orientation_mode
                                            ? std::max(
                                                by_orientation[0U],
                                                by_orientation[1U]
                                            )
                                            : negative
                                        );
                                    ++cases;
                                    if (!have_margin
                                        || margin < minimum_margin) {
                                        have_margin = true;
                                        minimum_margin = margin;
                                        minimum_minus = minus;
                                        minimum_plus = plus;
                                        minimum_invariant = invariant;
                                        minimum_negative = negative;
                                    }
                                    if (margin < 0) {
                                        std::cout
                                            << "SU2_SEVEN_TWO_MINUS_COUNT "
                                            << "result=FAIL";
                                        print_case(
                                            minus, plus, invariant, negative
                                        );
                                        std::cout << " cases=" << cases
                                                  << " orientation0="
                                                  << by_orientation[0U]
                                                  << " orientation1="
                                                  << by_orientation[1U]
                                                  << '\n';
                                        found = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (!found) {
            std::cout << "SU2_SEVEN_TWO_MINUS_COUNT result=PASS"
                      << " maximum_label=" << maximum_label
                      << " cases=" << cases
                      << " minimum";
            if (have_margin) {
                print_case(
                    minimum_minus, minimum_plus,
                    minimum_invariant, minimum_negative
                );
            } else {
                std::cout << " none";
            }
            std::cout << '\n';
        }
        return found ? EXIT_FAILURE : EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
