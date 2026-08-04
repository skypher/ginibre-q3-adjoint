#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Profile = std::vector<Integer>;

int positive(const char* text, const char* name) {
    const std::string value{text};
    std::size_t used = 0U;
    const long parsed = std::stol(value, &used);
    if (used != value.size() || parsed <= 0) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

Integer at(const Profile& profile, int radius, int index) {
    if (index < -radius || index > radius) {
        return 0;
    }
    return profile[static_cast<std::size_t>(index + radius)];
}

Profile append_factor(const Profile& profile, int radius, int factor) {
    const int next_radius = radius + factor;
    Profile next(static_cast<std::size_t>(2 * next_radius + 1));
    for (int source = -radius; source <= radius; ++source) {
        const Integer& coefficient = at(profile, radius, source);
        if (coefficient == 0) {
            continue;
        }
        for (int shift = -factor; shift <= factor; ++shift) {
            next[static_cast<std::size_t>(source + shift + next_radius)]
                += coefficient;
        }
    }
    return next;
}

Profile square_profile(const Profile& profile, int radius) {
    const int squared_radius = 2 * radius;
    Profile square(static_cast<std::size_t>(2 * squared_radius + 1));
    for (int left = -radius; left <= radius; ++left) {
        const Integer& left_coefficient = at(profile, radius, left);
        if (left_coefficient == 0) {
            continue;
        }
        for (int right = -radius; right <= radius; ++right) {
            square[static_cast<std::size_t>(
                left + right + squared_radius
            )] += left_coefficient * at(profile, radius, right);
        }
    }
    return square;
}

Integer psi(const Profile& profile, int radius, int row, int index) {
    return at(profile, radius, index - row)
        - at(profile, radius, index + row + 1);
}

std::string render(const std::vector<int>& word) {
    std::string result{"["};
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(word[index]);
    }
    return result + ']';
}

struct Scan {
    std::uint64_t words = 0U;
    std::uint64_t tails = 0U;
    bool found = false;
};

bool inspect_word(const std::vector<int>& word, Scan& scan) {
    Profile profile{Integer{1}};
    int radius = 0;
    for (const int factor : word) {
        profile = append_factor(profile, radius, factor);
        radius += factor;
    }

    const int maximum_index = 2 * radius;
    const Profile square = square_profile(profile, radius);
    const auto multiplicity = [&square, radius](int label) {
        return at(square, 2 * radius, label)
            - at(square, 2 * radius, label + 1);
    };
    const int side = maximum_index + 1;
    std::vector<std::vector<Integer>> wedges(
        static_cast<std::size_t>(radius + 1),
        std::vector<Integer>(static_cast<std::size_t>(side * side))
    );
    for (int row = 1; row <= radius; ++row) {
        for (int left = 0; left <= maximum_index; ++left) {
            const Integer left_zero = psi(profile, radius, 0, left);
            const Integer left_row = psi(profile, radius, row, left);
            for (int right = left + 1; right <= maximum_index; ++right) {
                wedges[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(left * side + right)]
                    = left_zero * psi(profile, radius, row, right)
                    - psi(profile, radius, 0, right) * left_row;
            }
        }
    }

    for (int first = 1; first <= radius; ++first) {
        for (int second = first; second <= radius; ++second) {
            Integer tail = 0;
            for (int cutoff = maximum_index - 1; cutoff >= 0; --cutoff) {
                for (int right = cutoff + 1;
                     right <= maximum_index;
                     ++right) {
                    const std::size_t coordinate = static_cast<std::size_t>(
                        cutoff * side + right
                    );
                    tail += wedges[static_cast<std::size_t>(first)][coordinate]
                          * wedges[static_cast<std::size_t>(second)][coordinate];
                }
                ++scan.tails;
                if (tail < 0) {
                    std::cout
                        << "SU2_ORDINARY_FACTOR_WORD_WEDGE_TAIL"
                        << " result=COUNTEREXAMPLE"
                        << " word=" << render(word)
                        << " R=" << first
                        << " S=" << second
                        << " cutoff=" << cutoff
                        << " value=" << tail << '\n';
                    scan.found = true;
                    return false;
                }
            }
            Integer fusion_sum = 0;
            for (int label = std::abs(first - second);
                 label <= first + second;
                 ++label) {
                fusion_sum += multiplicity(label);
            }
            const Integer direct = multiplicity(0) * fusion_sum
                - multiplicity(first) * multiplicity(second);
            if (tail != direct) {
                throw std::runtime_error(
                    "wedge-tail identity disagrees with the direct current"
                );
            }
        }
    }
    ++scan.words;
    return true;
}

bool enumerate(
    int remaining,
    int minimum_factor,
    int maximum_factor,
    std::vector<int>& word,
    Scan& scan
) {
    if (remaining == 0) {
        return inspect_word(word, scan);
    }
    for (int factor = minimum_factor; factor <= maximum_factor; ++factor) {
        word.push_back(factor);
        if (!enumerate(
                remaining - 1,
                factor,
                maximum_factor,
                word,
                scan
            )) {
            return false;
        }
        word.pop_back();
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_factor = argc >= 2
            ? positive(argv[1], "maximum factor")
            : 3;
        const int maximum_length = argc >= 3
            ? positive(argv[2], "maximum length")
            : 4;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_ordinary_factor_word_wedge_tail "
                "[maximum_factor] [maximum_length]"
            );
        }
        Scan scan;
        std::vector<int> word;
        for (int length = 1; length <= maximum_length; ++length) {
            if (!enumerate(
                    length,
                    1,
                    maximum_factor,
                    word,
                    scan
                )) {
                return EXIT_SUCCESS;
            }
        }
        std::cout
            << "SU2_ORDINARY_FACTOR_WORD_WEDGE_TAIL"
            << " result=NO_COUNTEREXAMPLE"
            << " maximum_factor=" << maximum_factor
            << " maximum_length=" << maximum_length
            << " words=" << scan.words
            << " tails=" << scan.tails << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_ORDINARY_FACTOR_WORD_WEDGE_TAIL error="
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
