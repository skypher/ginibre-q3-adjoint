#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Profile = std::vector<Integer>;

std::uint64_t parse_positive(const char* text, const char* name) {
    const std::string value(text);
    std::size_t used = 0U;
    const unsigned long long parsed = std::stoull(value, &used);
    if (used != value.size() || parsed == 0U) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<std::uint64_t>(parsed);
}

std::uint64_t splitmix64(std::uint64_t& state) {
    state += UINT64_C(0x9e3779b97f4a7c15);
    std::uint64_t value = state;
    value = (value ^ (value >> 30U))
        * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U))
        * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

Profile fuse(const Profile& input, int level, int factor) {
    Profile output(static_cast<std::size_t>(level + 1), 0);
    for (int source = 0; source <= level; ++source) {
        const Integer& multiplicity =
            input[static_cast<std::size_t>(source)];
        if (multiplicity == 0) {
            continue;
        }
        const int lower = std::abs(source - factor);
        const int upper = std::min(
            source + factor,
            2 * level - source - factor
        );
        for (int target = lower; target <= upper; target += 2) {
            output[static_cast<std::size_t>(target)] += multiplicity;
        }
    }
    return output;
}

Profile square_profile(const Profile& root, int level) {
    Profile result(static_cast<std::size_t>(level + 1), 0);
    for (int label = 0; label <= level; ++label) {
        const Integer& multiplicity = root[static_cast<std::size_t>(label)];
        if (multiplicity == 0) {
            continue;
        }
        const Profile translated = fuse(root, level, label);
        for (int target = 0; target <= level; ++target) {
            result[static_cast<std::size_t>(target)] += multiplicity
                * translated[static_cast<std::size_t>(target)];
        }
    }
    return result;
}

Integer folded_character_coefficient(
    const Profile& profile,
    int degree,
    int level
) {
    if (degree < 0) {
        return 0;
    }
    const int period = 2 * (level + 2);
    const int residue = degree % period;
    if (residue <= level) {
        return profile[static_cast<std::size_t>(residue)];
    }
    if (residue == level + 1) {
        return 0;
    }
    const int reflected = period - residue - 2;
    return reflected < 0
        ? Integer(0)
        : -profile[static_cast<std::size_t>(reflected)];
}

Integer c2_packet_coordinate(
    const Profile& square,
    int first,
    int second,
    int level
) {
    if (first < second || second < 0) {
        throw std::runtime_error("invalid C2 packet coordinate");
    }
    const Integer d_first = folded_character_coefficient(
        square, first, level
    );
    const Integer d_first_plus_two = folded_character_coefficient(
        square, first + 2, level
    );
    const Integer d_first_plus_one = folded_character_coefficient(
        square, first + 1, level
    );
    const Integer d_second = folded_character_coefficient(
        square, second, level
    );
    const Integer d_second_minus_one = folded_character_coefficient(
        square, second - 1, level
    );
    const Integer d_second_plus_one = folded_character_coefficient(
        square, second + 1, level
    );
    return (d_first + d_first_plus_two) * d_second
        - d_first_plus_one * (d_second_minus_one + d_second_plus_one);
}

Integer anchored_current(
    const Profile& square,
    int left,
    int right,
    int level
) {
    const int lower = std::abs(left - right);
    const int upper = std::min(
        left + right,
        2 * level - left - right
    );
    Integer fusion_entry = 0;
    for (int label = lower; label <= upper; label += 2) {
        fusion_entry += square[static_cast<std::size_t>(label)];
    }
    return square[0] * fusion_entry
        - square[static_cast<std::size_t>(left)]
            * square[static_cast<std::size_t>(right)];
}

Integer complete_c2_packet(
    const Profile& square,
    int left,
    int right,
    int level
) {
    Integer total = 0;
    for (int index = 0; index < left; ++index) {
        for (int contraction = 0;
             contraction <= left - 1 - index;
             ++contraction) {
            total += c2_packet_coordinate(
                square,
                left + right - 2 - index - 2 * contraction,
                index,
                level
            );
        }
    }
    return total;
}

Integer radial_c2_packet(
    const Profile& square,
    int total_label,
    int contraction_depth,
    int level
) {
    Integer total = 0;
    for (int index = 0; index <= contraction_depth; ++index) {
        total += c2_packet_coordinate(
            square,
            total_label - index,
            index,
            level
        );
    }
    return total;
}

Integer radial_c2_closed_form(
    const Profile& square,
    int total_label,
    int contraction_depth,
    int level
) {
    return folded_character_coefficient(square, 0, level)
             * folded_character_coefficient(square, total_label + 2, level)
        + folded_character_coefficient(
              square, contraction_depth, level
          ) * folded_character_coefficient(
              square, total_label - contraction_depth, level
          )
        - folded_character_coefficient(
              square, contraction_depth + 1, level
          ) * folded_character_coefficient(
              square, total_label - contraction_depth + 1, level
          );
}

std::string render_word(const std::vector<int>& word) {
    std::string result = "[";
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(word[index]);
    }
    return result + ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error(
                "usage: verify_su2_finite_c2_packet "
                "SAMPLES MAXIMUM_LEVEL MAXIMUM_FACTORS"
            );
        }
        const std::uint64_t samples = parse_positive(argv[1], "samples");
        const int maximum_level = static_cast<int>(
            parse_positive(argv[2], "maximum level")
        );
        const int maximum_factors = static_cast<int>(
            parse_positive(argv[3], "maximum factors")
        );
        std::uint64_t current_checks = 0U;
        std::uint64_t packet_coordinate_checks = 0U;
        std::uint64_t radial_packet_checks = 0U;
        std::uint64_t radial_packet_negatives = 0U;
        bool printed_first_negative_radial_packet = false;
        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            std::uint64_t state = sample ^ UINT64_C(0x6a09e667f3bcc909);
            const int level = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level)
            );
            const int factors = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_factors)
            );
            Profile root(static_cast<std::size_t>(level + 1), 0);
            root[0] = 1;
            std::vector<int> word;
            word.reserve(static_cast<std::size_t>(factors));
            for (int index = 0; index < factors; ++index) {
                const int factor = 1 + static_cast<int>(
                    splitmix64(state)
                    % static_cast<std::uint64_t>(level)
                );
                word.push_back(factor);
                root = fuse(root, level, factor);
            }
            const Profile square = square_profile(root, level);
            for (int left = 1; left <= level; ++left) {
                for (int right = left; right <= level; ++right) {
                    const Integer direct = anchored_current(
                        square, left, right, level
                    );
                    const Integer packet = complete_c2_packet(
                        square, left, right, level
                    );
                    ++current_checks;
                    packet_coordinate_checks += static_cast<std::uint64_t>(
                        left * (left + 1) / 2
                    );
                    if (direct != packet) {
                        std::cerr
                            << "FCP_MISMATCH level=" << level
                            << " word=" << render_word(word)
                            << " R=" << left
                            << " S=" << right
                            << " direct=" << direct
                            << " packet=" << packet << '\n';
                        return EXIT_FAILURE;
                    }
                    for (int contraction = 0;
                         contraction < left;
                         ++contraction) {
                        const Integer radial = radial_c2_packet(
                            square,
                            left + right - 2 - 2 * contraction,
                            left - 1 - contraction,
                            level
                        );
                        const Integer radial_closed = radial_c2_closed_form(
                            square,
                            left + right - 2 - 2 * contraction,
                            left - 1 - contraction,
                            level
                        );
                        ++radial_packet_checks;
                        if (radial != radial_closed) {
                            std::cerr
                                << "FCP_RADIAL_MISMATCH level=" << level
                                << " word=" << render_word(word)
                                << " R=" << left
                                << " S=" << right
                                << " A="
                                << left + right - 2 - 2 * contraction
                                << " L=" << left - 1 - contraction
                                << " direct=" << radial
                                << " closed=" << radial_closed << '\n';
                            return EXIT_FAILURE;
                        }
                        if (radial < 0) {
                            ++radial_packet_negatives;
                            if (!printed_first_negative_radial_packet) {
                                printed_first_negative_radial_packet = true;
                                std::cout
                                    << "FIRST_NEGATIVE_RADIAL_PACKET"
                                    << " level=" << level
                                    << " word=" << render_word(word)
                                    << " R=" << left
                                    << " S=" << right
                                    << " A="
                                    << left + right - 2
                                        - 2 * contraction
                                    << " L=" << left - 1 - contraction
                                    << " value=" << radial
                                    << " full_current=" << direct << '\n';
                            }
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_FINITE_C2_PACKET"
            << " samples=" << samples
            << " maximum_level=" << maximum_level
            << " maximum_factors=" << maximum_factors
            << " current_checks=" << current_checks
            << " packet_coordinate_checks=" << packet_coordinate_checks
            << " radial_packet_checks=" << radial_packet_checks
            << " radial_packet_negatives=" << radial_packet_negatives
            << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
