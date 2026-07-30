#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Vector = std::vector<long long>;

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0
        || parsed > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

bool selected_slope(int index, int level, int label) {
    if (index < 0) {
        return false;
    }
    const int period = 2 * level + 2;
    if (label == level) {
        return index >= period + level
            && index % period == level;
    }
    if (index < period + label) {
        return false;
    }
    const int residue = index % period;
    return residue == label || residue == period - label - 2;
}

long long reserve_coefficient(int index, int level, int label) {
    return static_cast<long long>(selected_slope(index, level, label))
        - static_cast<long long>(
            selected_slope(index - 1, level, label));
}

long long multiplied_coefficient(
    int output,
    int level,
    int factor,
    int label) {
    long long value = 0;
    const int maximum_character = 2 * factor;
    const int maximum_input = output + maximum_character;
    for (int input = 0; input <= maximum_input; ++input) {
        const long long coefficient
            = reserve_coefficient(input, level, label);
        if (coefficient == 0) {
            continue;
        }
        for (int character = 0;
             character <= maximum_character;
             ++character) {
            if (std::abs(input - character) <= output
                && output <= input + character) {
                value += coefficient;
            }
        }
    }
    return value;
}

std::string render(const Vector& values) {
    std::string result = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(values[index]);
    }
    return result + ']';
}

struct Failure {
    int level = -1;
    int factor = -1;
    int label = -1;
    int index = -1;
    long long value = 0;
    Vector tail_coefficients;
    Vector residual;
};

void record(
    Failure& failure,
    int level,
    int factor,
    int label,
    int index,
    long long value,
    const Vector& tail_coefficients,
    const Vector& residual) {
    if (failure.level >= 0) {
        return;
    }
    failure = {
        level,
        factor,
        label,
        index,
        value,
        tail_coefficients,
        residual};
}

std::string render_failure(const Failure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " factor=" + std::to_string(failure.factor)
        + " label=" + std::to_string(failure.label)
        + " index=" + std::to_string(failure.index)
        + " value=" + std::to_string(failure.value)
        + " tail_coefficients=" + render(failure.tail_coefficients)
        + " residual=" + render(failure.residual);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 30;
        if (argc > 2 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim1_generator "
                "[maximum_half_level]");
        }

        std::uint64_t decompositions = 0U;
        std::uint64_t tail_checks = 0U;
        std::uint64_t tail_failures = 0U;
        std::uint64_t residual_checks = 0U;
        std::uint64_t residual_failures = 0U;
        Failure first_tail;
        Failure first_residual;

        for (int level = 2; level <= maximum_level; ++level) {
            const int period = 2 * level + 2;
            const int horizon = 6 * period;
            const int reference_block = 4;
            for (int factor = 1; factor <= level / 2; ++factor) {
                for (int label = 0; label <= level; ++label) {
                    Vector product(static_cast<std::size_t>(horizon + 1), 0);
                    for (int index = 0; index <= horizon; ++index) {
                        product[static_cast<std::size_t>(index)]
                            = multiplied_coefficient(
                                index,
                                level,
                                factor,
                                label);
                    }

                    Vector tail_coefficients(
                        static_cast<std::size_t>(level + 1),
                        0);
                    long long cumulative = 0;
                    for (int residue = 0; residue <= level; ++residue) {
                        cumulative += product[static_cast<std::size_t>(
                            reference_block * period + residue)];
                        tail_coefficients[
                            static_cast<std::size_t>(residue)] = cumulative;
                        ++tail_checks;
                        if (cumulative < 0) {
                            ++tail_failures;
                        }
                    }

                    Vector residual(static_cast<std::size_t>(horizon + 1), 0);
                    for (int index = 0; index <= horizon; ++index) {
                        long long value
                            = product[static_cast<std::size_t>(index)];
                        for (int source = 0; source <= level; ++source) {
                            value -= tail_coefficients[
                                         static_cast<std::size_t>(source)]
                                * reserve_coefficient(
                                    index,
                                    level,
                                    source);
                        }
                        residual[static_cast<std::size_t>(index)] = value;
                        ++residual_checks;
                        if (value < 0) {
                            ++residual_failures;
                            record(
                                first_residual,
                                level,
                                factor,
                                label,
                                index,
                                value,
                                tail_coefficients,
                                residual);
                        }
                    }
                    if (first_tail.level < 0) {
                        for (int source = 0;
                             source <= level;
                             ++source) {
                            const long long value = tail_coefficients[
                                static_cast<std::size_t>(source)];
                            if (value < 0) {
                                record(
                                    first_tail,
                                    level,
                                    factor,
                                    label,
                                    source,
                                    value,
                                    tail_coefficients,
                                    residual);
                                break;
                            }
                        }
                    }
                    if (maximum_level <= 4) {
                        std::cout
                            << "DETAIL"
                            << " level=" << 2 * level
                            << " factor=" << factor
                            << " label=" << label
                            << " tail_coefficients="
                            << render(tail_coefficients)
                            << " residual=" << render(residual)
                            << '\n';
                    }
                    ++decompositions;
                }
            }
        }

        std::cout
            << "SU2_AIM1_GENERATOR"
            << " maximum_level=" << 2 * maximum_level
            << " decompositions=" << decompositions
            << " tail_checks=" << tail_checks
            << " tail_failures=" << tail_failures
            << " residual_checks=" << residual_checks
            << " residual_failures=" << residual_failures
            << '\n'
            << "FIRST_TAIL_FAILURE "
            << render_failure(first_tail) << '\n'
            << "FIRST_RESIDUAL_FAILURE "
            << render_failure(first_residual) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
