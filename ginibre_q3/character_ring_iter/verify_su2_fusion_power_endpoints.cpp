#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("bound must be a positive integer");
    }
    return static_cast<int>(value);
}

Integer binomial_integer(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

Integer weight_multiplicity(int power, int label, int depth) {
    if (depth < 0) {
        return 0;
    }
    Integer result = 0;
    for (int image = 0; image <= power; ++image) {
        Integer term =
            binomial_integer(power, image)
            * binomial_integer(
                depth - image * (label + 1) + power - 1,
                power - 1
            );
        if (image % 2 == 0) {
            result += term;
        } else {
            result -= term;
        }
    }
    return result;
}

std::pair<int, int> affine_fold(int level, int label) {
    const int period_half = level + 2;
    const int period = 2 * period_half;
    const int residue = (label + 1) % period;
    if (residue == 0 || residue == period_half) {
        return {0, 0};
    }
    if (residue < period_half) {
        return {1, residue - 1};
    }
    return {-1, period - residue - 1};
}

bool fuses(
    int level,
    int left,
    int right,
    int target
) {
    return
        (left + right + target) % 2 == 0
        && std::abs(left - right) <= target
        && target <= left + right
        && left + right + target <= 2 * level;
}

Integer endpoint_formula(
    int level,
    int label,
    int power,
    int source,
    int target,
    int& runs,
    int& negative_runs
) {
    const int total = power * label;
    Integer result = 0;
    int active_sign = 0;
    int lower = 0;
    runs = 0;
    negative_runs = 0;
    const auto close_run =
        [
            &result,
            &runs,
            &negative_runs,
            power,
            label
        ](int sign, int first, int last) {
            if (sign == 0) {
                return;
            }
            const Integer interval =
                weight_multiplicity(power, label, last)
                - weight_multiplicity(power, label, first - 1);
            result += sign * interval;
            ++runs;
            if (sign < 0) {
                ++negative_runs;
            }
        };
    for (int depth = 0; depth <= total / 2; ++depth) {
        const int classical_label = total - 2 * depth;
        const auto [sign, folded_label] =
            affine_fold(level, classical_label);
        const int next_sign =
            sign != 0
                && fuses(
                    level,
                    source,
                    folded_label,
                    target
                )
                ? sign
                : 0;
        if (next_sign != active_sign) {
            close_run(active_sign, lower, depth - 1);
            active_sign = next_sign;
            lower = depth;
        }
    }
    close_run(active_sign, lower, total / 2);
    return result;
}

std::vector<Integer> multiply_fusion(
    const std::vector<Integer>& state,
    int level,
    int label
) {
    std::vector<Integer> next(
        static_cast<std::size_t>(level + 1)
    );
    for (int source = 0; source <= level; ++source) {
        if (state[static_cast<std::size_t>(source)] == 0) {
            continue;
        }
        const int lower = std::abs(source - label);
        const int upper = std::min(
            source + label,
            2 * level - source - label
        );
        for (int target = lower; target <= upper; target += 2) {
            next[static_cast<std::size_t>(target)] +=
                state[static_cast<std::size_t>(source)];
        }
    }
    return next;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_LEVEL");
        }
        const int maximum_level = parse_positive(argv[1]);
        std::uint64_t parameters = 0U;
        std::uint64_t coordinates = 0U;
        std::uint64_t signed_endpoint_runs = 0U;
        std::uint64_t negative_endpoint_runs = 0U;
        int maximum_runs = 0;
        Integer maximum_entry = 0;

        for (int level = 1; level <= maximum_level; ++level) {
            for (int label = 1; label <= level; ++label) {
                for (int source = 0; source <= level; ++source) {
                    ++parameters;
                    std::vector<Integer> state(
                        static_cast<std::size_t>(level + 1)
                    );
                    state[static_cast<std::size_t>(source)] = 1;
                    for (int power = 1; power <= 4; ++power) {
                        state = multiply_fusion(state, level, label);
                        for (int target = 0;
                             target <= level;
                             ++target) {
                            ++coordinates;
                            int runs = 0;
                            int negative_runs = 0;
                            const Integer expected = endpoint_formula(
                                level,
                                label,
                                power,
                                source,
                                target,
                                runs,
                                negative_runs
                            );
                            const Integer actual =
                                state[static_cast<std::size_t>(target)];
                            if (actual != expected || expected < 0) {
                                std::cerr
                                    << "FAILED_FUSION_POWER_ENDPOINT"
                                    << " ell=" << level
                                    << " d=" << label
                                    << " b=" << power
                                    << " s=" << source
                                    << " t=" << target
                                    << " actual=" << actual
                                    << " expected=" << expected << '\n';
                                throw std::runtime_error(
                                    "fusion-power endpoint mismatch"
                                );
                            }
                            maximum_entry = std::max(
                                maximum_entry,
                                expected
                            );
                            signed_endpoint_runs +=
                                static_cast<std::uint64_t>(runs);
                            negative_endpoint_runs +=
                                static_cast<std::uint64_t>(
                                    negative_runs
                                );
                            maximum_runs = std::max(
                                maximum_runs,
                                runs
                            );
                            if (runs > 4) {
                                throw std::runtime_error(
                                    "more than four affine endpoint runs"
                                );
                            }
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_FUSION_POWER_ENDPOINTS"
            << " maximum_level=" << maximum_level
            << " parameters=" << parameters
            << " coordinates=" << coordinates
            << " signed_endpoint_runs=" << signed_endpoint_runs
            << " negative_endpoint_runs=" << negative_endpoint_runs
            << " maximum_runs=" << maximum_runs
            << " maximum_entry=" << maximum_entry
            << " result=PASS_EXACT_ENDPOINT_FORMULA\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
