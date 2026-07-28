#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

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

Integer hinge_formula(int power, int label, int depth) {
    Integer result = 0;
    for (int image = 0; image <= power; ++image) {
        Integer term =
            binomial_integer(power, image)
            * binomial_integer(
                depth - image * (label + 1) + power - 2,
                power - 2
            );
        if (image % 2 == 0) {
            result += term;
        } else {
            result -= term;
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_LABEL");
        }
        const int maximum_label = parse_positive(argv[1]);
        std::uint64_t coordinates = 0U;
        std::uint64_t active_terms = 0U;
        Integer maximum_multiplicity = 0;

        for (int label = 1; label <= maximum_label; ++label) {
            std::map<int, Integer> decomposition{{0, 1}};
            for (int power = 1; power <= 4; ++power) {
                std::map<int, Integer> next;
                for (const auto& [source, multiplicity]
                     : decomposition) {
                    for (int target = std::abs(source - label);
                         target <= source + label;
                         target += 2) {
                        next[target] += multiplicity;
                    }
                }
                decomposition = std::move(next);
                if (power < 2) {
                    continue;
                }
                for (int depth = 0;
                     depth <= power * label / 2;
                     ++depth) {
                    ++coordinates;
                    const int target = power * label - 2 * depth;
                    const Integer expected =
                        hinge_formula(power, label, depth);
                    const auto found = decomposition.find(target);
                    const Integer actual =
                        found == decomposition.end()
                            ? Integer(0)
                            : found->second;
                    if (actual != expected || expected < 0) {
                        std::cerr
                            << "FAILED_TENSOR_POWER_HINGE"
                            << " d=" << label
                            << " b=" << power
                            << " r=" << depth
                            << " actual=" << actual
                            << " expected=" << expected << '\n';
                        throw std::runtime_error(
                            "tensor-power hinge formula mismatch"
                        );
                    }
                    maximum_multiplicity = std::max(
                        maximum_multiplicity,
                        expected
                    );
                    for (int image = 1; image <= power; ++image) {
                        if (depth - image * (label + 1) >= 0) {
                            ++active_terms;
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_TENSOR_POWER_HINGES"
            << " maximum_label=" << maximum_label
            << " coordinates=" << coordinates
            << " active_truncation_terms=" << active_terms
            << " maximum_multiplicity=" << maximum_multiplicity
            << " result=PASS_EXACT_HINGE_FORMULA\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
