#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

std::vector<cpp_int> partial_character_coefficients(
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    const int width = total + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(width * width)
    );
    current[0] = 1;
    int support = 0;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int label = labels[index];
        const bool minus = ((minus_mask >> index) & 1U) != 0U;
        std::vector<cpp_int> next(current.size());
        for (int left = 0; left <= support; ++left) {
            for (int right = 0; right <= support; ++right) {
                const cpp_int& value = current[static_cast<std::size_t>(
                    left * width + right
                )];
                if (value == 0) {
                    continue;
                }
                for (int output = std::abs(left - label);
                     output <= left + label;
                     output += 2) {
                    next[static_cast<std::size_t>(
                        output * width + right
                    )] += value;
                }
                for (int output = std::abs(right - label);
                     output <= right + label;
                     output += 2) {
                    next[static_cast<std::size_t>(
                        left * width + output
                    )] += minus ? -value : value;
                }
            }
        }
        support += label;
        current = std::move(next);
    }
    std::vector<cpp_int> result(static_cast<std::size_t>(width));
    for (int target = 0; target <= total; ++target) {
        result[static_cast<std::size_t>(target)]
            = current[static_cast<std::size_t>(target * width)];
    }
    return result;
}

bool inspect_word(
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int& failed_target,
    cpp_int& failed_value
) {
    const std::vector<cpp_int> coefficients
        = partial_character_coefficients(labels, minus_mask);
    for (std::size_t target = 0U;
         target < coefficients.size(); ++target) {
        const cpp_int& value = coefficients[target];
        if (value < 0) {
            failed_target = static_cast<int>(target);
            failed_value = value;
            return false;
        }
    }
    return true;
}

cpp_int binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    k = std::min(k, n - k);
    cpp_int result = 1;
    for (int index = 1; index <= k; ++index) {
        result *= n - k + index;
        result /= index;
    }
    return result;
}

bool inspect_gamma(
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int& failed_index,
    cpp_int& failed_value
) {
    int total = 0;
    for (const int label : labels) {
        total += label;
    }
    const std::vector<cpp_int> coefficients
        = partial_character_coefficients(labels, minus_mask);
    const int middle = total / 2;
    std::vector<cpp_int> hilbert(
        static_cast<std::size_t>(total + 1), 0
    );
    cpp_int cumulative = 0;
    for (int degree = 0; degree <= middle; ++degree) {
        cumulative += coefficients[static_cast<std::size_t>(
            total - 2 * degree
        )];
        hilbert[static_cast<std::size_t>(degree)] = cumulative;
        hilbert[static_cast<std::size_t>(total - degree)]
            = cumulative;
    }
    std::vector<cpp_int> gamma(
        static_cast<std::size_t>(middle + 1), 0
    );
    for (int index = 0; index <= middle; ++index) {
        cpp_int value = hilbert[static_cast<std::size_t>(index)];
        for (int earlier = 0; earlier < index; ++earlier) {
            value -= gamma[static_cast<std::size_t>(earlier)]
                * binomial(
                    total - 2 * earlier, index - earlier
                );
        }
        gamma[static_cast<std::size_t>(index)] = value;
        if (value < 0) {
            failed_index = index;
            failed_value = value;
            return false;
        }
    }
    return true;
}

bool inspect_log_concavity(
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int& failed_index,
    cpp_int& failed_value
) {
    int total = 0;
    for (const int label : labels) {
        total += label;
    }
    const std::vector<cpp_int> coefficients
        = partial_character_coefficients(labels, minus_mask);
    const int middle = total / 2;
    std::vector<cpp_int> hilbert(
        static_cast<std::size_t>(total + 1), 0
    );
    cpp_int cumulative = 0;
    for (int degree = 0; degree <= middle; ++degree) {
        cumulative += coefficients[static_cast<std::size_t>(
            total - 2 * degree
        )];
        hilbert[static_cast<std::size_t>(degree)] = cumulative;
        hilbert[static_cast<std::size_t>(total - degree)]
            = cumulative;
    }
    for (int index = 1; index < total; ++index) {
        const cpp_int defect
            = hilbert[static_cast<std::size_t>(index)]
                * hilbert[static_cast<std::size_t>(index)]
              - hilbert[static_cast<std::size_t>(index - 1)]
                * hilbert[static_cast<std::size_t>(index + 1)];
        if (defect < 0) {
            failed_index = index;
            failed_value = defect;
            return false;
        }
    }
    return true;
}

void print_word(const std::vector<int>& labels, std::uint64_t minus_mask) {
    std::cout << '[';
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << (((minus_mask >> index) & 1U) != 0U ? '-' : '+')
                  << labels[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: search_su2_partial_character_cone MAXIMUM_LABEL "
                "MAXIMUM_FACTORS [gamma|logconcave]"
            );
        }
        const std::string mode = argc == 4
            ? std::string(argv[3]) : std::string("character");
        if (mode != "character" && mode != "gamma"
            && mode != "logconcave") {
            throw std::runtime_error(
                "optional mode must be gamma or logconcave"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 1
            || maximum_factors >= 63) {
            throw std::runtime_error("invalid search bound");
        }

        std::uint64_t tested = 0;
        bool failed = false;
        std::vector<int> labels;
        const auto visit = [&](const auto& self, int first) -> void {
            if (failed) {
                return;
            }
            if (!labels.empty()) {
                const std::uint64_t masks
                    = std::uint64_t{1} << labels.size();
                for (std::uint64_t mask = 0; mask < masks; ++mask) {
                    int failed_target = -1;
                    cpp_int failed_value = 0;
                    ++tested;
                    const bool passed = mode == "gamma"
                        ? inspect_gamma(
                            labels, mask, failed_target, failed_value
                        ) : mode == "logconcave"
                        ? inspect_log_concavity(
                            labels, mask, failed_target, failed_value
                        ) : inspect_word(
                            labels, mask, failed_target, failed_value
                        );
                    if (!passed) {
                        std::cout
                            << (mode == "gamma"
                                ? "GAMMA_FAIL factors="
                                : mode == "logconcave"
                                ? "LOGCONCAVITY_FAIL factors="
                                : "FAIL factors=");
                        print_word(labels, mask);
                        std::cout
                                  << (mode == "gamma"
                                      ? " gamma_index="
                                      : mode == "logconcave"
                                      ? " coefficient_index="
                                      : " target=")
                                  << failed_target
                                  << " coefficient=" << failed_value
                                  << " tested=" << tested << '\n';
                        failed = true;
                        return;
                    }
                }
            }
            if (labels.size() == static_cast<std::size_t>(maximum_factors)) {
                return;
            }
            for (int label = first; label <= maximum_label; ++label) {
                labels.push_back(label);
                self(self, label);
                labels.pop_back();
                if (failed) {
                    return;
                }
            }
        };
        visit(visit, 1);
        if (failed) {
            return EXIT_FAILURE;
        }
        std::cout
                  << (mode == "gamma"
                      ? "SU2_PARTIAL_CHARACTER_GAMMA PASS tested="
                      : mode == "logconcave"
                      ? "SU2_PARTIAL_CHARACTER_LOGCONCAVITY PASS tested="
                      : "SU2_PARTIAL_CHARACTER_CONE PASS tested=")
                  << tested
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
