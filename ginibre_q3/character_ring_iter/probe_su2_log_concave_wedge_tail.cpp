#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

Integer coefficient_at(const std::vector<Integer>& profile, int index) {
    if (index < 0 || index >= static_cast<int>(profile.size())) {
        return 0;
    }
    return profile[static_cast<std::size_t>(index)];
}

bool is_log_concave(const std::vector<Integer>& profile) {
    for (std::size_t index = 1; index + 1 < profile.size(); ++index) {
        if (
            profile[index] * profile[index]
            < profile[index - 1] * profile[index + 1]
        ) {
            return false;
        }
    }
    return true;
}

std::vector<Integer> multiply_character(
    const std::vector<Integer>& profile,
    int label
) {
    std::vector<Integer> result(
        profile.size() + static_cast<std::size_t>(label)
    );
    for (int target = 0; target < static_cast<int>(result.size());
         ++target) {
        const int lower = std::abs(target - label);
        const int upper = target + label;
        for (int source = lower; source <= upper; ++source) {
            result[static_cast<std::size_t>(target)] +=
                coefficient_at(profile, source);
        }
    }
    return result;
}

Integer inner_tail(
    const std::vector<Integer>& left,
    const std::vector<Integer>& right,
    int cutoff
) {
    const int size = std::max(
        static_cast<int>(left.size()),
        static_cast<int>(right.size())
    );
    Integer result = 0;
    for (int index = cutoff; index < size; ++index) {
        result +=
            coefficient_at(left, index)
            * coefficient_at(right, index);
    }
    return result;
}

bool check_profile(
    const std::vector<Integer>& profile,
    int max_label,
    std::size_t& tails,
    std::size_t& transforms,
    std::size_t& gap_tails,
    bool& gap_counterexample,
    std::size_t& complete_wall_gap_tails,
    bool& complete_wall_gap_counterexample,
    std::size_t& gap_suffix_tails,
    bool& gap_suffix_counterexample,
    std::size_t& one_wedge_tails,
    bool& one_wedge_counterexample,
    std::size_t& adjacent_gap_payments,
    bool& adjacent_gap_counterexample
) {
    const Integer norm = inner_tail(profile, profile, 0);
    if (norm == 0) {
        return true;
    }
    for (int q_half = 1; q_half <= max_label; ++q_half) {
        const std::vector<Integer> q_profile =
            multiply_character(profile, q_half);
        ++transforms;
        if (!is_log_concave(q_profile)) {
            std::cout
                << "SU2_LOG_CONCAVE_WEDGE_TAIL"
                << " transform_counterexample"
                << " profile={";
            for (std::size_t index = 0;
                 index < profile.size();
                 ++index) {
                if (index != 0) {
                    std::cout << ',';
                }
                std::cout << profile[index];
            }
            std::cout
                << "}"
                << " label=" << q_half
                << '\n';
            return false;
        }
        const int q_size = std::max(
            static_cast<int>(profile.size()),
            static_cast<int>(q_profile.size())
        );
        for (int cutoff = 0; cutoff <= q_size; ++cutoff) {
            for (int gap = 1; gap < q_size; ++gap) {
                Integer wedge_tail = 0;
                for (
                    int left = cutoff;
                    left + gap < q_size;
                    ++left
                ) {
                    const int right = left + gap;
                    wedge_tail +=
                        coefficient_at(profile, left)
                            * coefficient_at(q_profile, right)
                        - coefficient_at(profile, right)
                            * coefficient_at(q_profile, left);
                }
                ++one_wedge_tails;
                if (
                    wedge_tail < 0
                    && !one_wedge_counterexample
                ) {
                    one_wedge_counterexample = true;
                    std::cout
                        << "SU2_LOG_CONCAVE_WEDGE_TAIL"
                        << " one_wedge_counterexample"
                        << " profile={";
                    for (std::size_t index = 0;
                         index < profile.size();
                         ++index) {
                        if (index != 0) {
                            std::cout << ',';
                        }
                        std::cout << profile[index];
                    }
                    std::cout
                        << "}"
                        << " label=" << q_half
                        << " cutoff=" << cutoff
                        << " gap=" << gap
                        << " value=" << wedge_tail
                        << '\n';
                }
            }
        }
        for (int target = 1; target <= max_label; ++target) {
            const std::vector<Integer> target_profile =
                multiply_character(profile, target);
            const int max_cutoff = std::max(
                static_cast<int>(q_profile.size()),
                static_cast<int>(target_profile.size())
            );
            for (int cutoff = 0; cutoff <= max_cutoff; ++cutoff) {
                const Integer pp =
                    inner_tail(profile, profile, cutoff);
                const Integer qt =
                    inner_tail(q_profile, target_profile, cutoff);
                const Integer pq =
                    inner_tail(profile, q_profile, cutoff);
                const Integer pt =
                    inner_tail(profile, target_profile, cutoff);
                const Integer determinant = pp * qt - pq * pt;
                ++tails;
                if (determinant < 0) {
                    std::cout
                        << "SU2_LOG_CONCAVE_WEDGE_TAIL"
                        << " counterexample"
                        << " profile={";
                    for (std::size_t index = 0;
                         index < profile.size();
                         ++index) {
                        if (index != 0) {
                            std::cout << ',';
                        }
                        std::cout << profile[index];
                    }
                    std::cout
                        << "}"
                        << " q_half=" << q_half
                        << " target=" << target
                        << " cutoff=" << cutoff
                        << " value=" << determinant
                        << '\n';
                    return false;
                }
                std::vector<Integer> gap_determinants(
                    static_cast<std::size_t>(max_cutoff),
                    Integer{0}
                );
                for (int gap = 1; gap < max_cutoff; ++gap) {
                    Integer gap_determinant = 0;
                    for (
                        int left = cutoff;
                        left + gap < max_cutoff;
                        ++left
                    ) {
                        const int right = left + gap;
                        const Integer q_wedge =
                            coefficient_at(profile, left)
                                * coefficient_at(q_profile, right)
                            - coefficient_at(profile, right)
                                * coefficient_at(q_profile, left);
                        const Integer target_wedge =
                            coefficient_at(profile, left)
                                * coefficient_at(
                                    target_profile,
                                    right
                                )
                            - coefficient_at(profile, right)
                                * coefficient_at(
                                    target_profile,
                                    left
                                );
                        gap_determinant += q_wedge * target_wedge;
                    }
                    gap_determinants[static_cast<std::size_t>(gap)] =
                        gap_determinant;
                    ++gap_tails;
                    if (
                        gap_determinant < 0
                        && !gap_counterexample
                    ) {
                        gap_counterexample = true;
                        std::cout
                            << "SU2_LOG_CONCAVE_WEDGE_TAIL"
                            << " gap_counterexample"
                            << " profile={";
                        for (std::size_t index = 0;
                             index < profile.size();
                             ++index) {
                            if (index != 0) {
                                std::cout << ',';
                            }
                            std::cout << profile[index];
                        }
                        std::cout
                            << "}"
                            << " q_half=" << q_half
                            << " target=" << target
                            << " cutoff=" << cutoff
                            << " gap=" << gap
                            << " value=" << gap_determinant
                            << '\n';
                    }
                    if (cutoff == 0) {
                        ++complete_wall_gap_tails;
                        if (
                            gap_determinant < 0
                            && !complete_wall_gap_counterexample
                        ) {
                            complete_wall_gap_counterexample = true;
                            std::cout
                                << "SU2_LOG_CONCAVE_WEDGE_TAIL"
                                << " complete_wall_gap_counterexample"
                                << " profile={";
                            for (std::size_t index = 0;
                                 index < profile.size();
                                 ++index) {
                                if (index != 0) {
                                    std::cout << ',';
                                }
                                std::cout << profile[index];
                            }
                            std::cout
                                << "}"
                                << " q_half=" << q_half
                                << " target=" << target
                                << " gap=" << gap
                                << " value=" << gap_determinant
                                << '\n';
                        }
                    }
                }
                Integer gap_suffix = 0;
                for (int gap = max_cutoff - 1; gap >= 1; --gap) {
                    gap_suffix +=
                        gap_determinants[static_cast<std::size_t>(gap)];
                    ++gap_suffix_tails;
                    if (
                        gap_suffix < 0
                        && !gap_suffix_counterexample
                    ) {
                        gap_suffix_counterexample = true;
                        std::cout
                            << "SU2_LOG_CONCAVE_WEDGE_TAIL"
                            << " gap_suffix_counterexample"
                            << " profile={";
                        for (std::size_t index = 0;
                             index < profile.size();
                             ++index) {
                            if (index != 0) {
                                std::cout << ',';
                            }
                            std::cout << profile[index];
                        }
                        std::cout
                            << "}"
                            << " q_half=" << q_half
                            << " target=" << target
                            << " cutoff=" << cutoff
                            << " minimum_gap=" << gap
                            << " value=" << gap_suffix
                            << '\n';
                    }
                }
                for (int gap = 1; gap + 1 < max_cutoff; ++gap) {
                    const Integer adjacent_payment =
                        gap_determinants[static_cast<std::size_t>(gap)]
                        + gap_determinants[
                            static_cast<std::size_t>(gap + 1)
                        ];
                    ++adjacent_gap_payments;
                    if (
                        adjacent_payment < 0
                        && !adjacent_gap_counterexample
                    ) {
                        adjacent_gap_counterexample = true;
                        std::cout
                            << "SU2_LOG_CONCAVE_WEDGE_TAIL"
                            << " adjacent_gap_counterexample"
                            << " profile={";
                        for (std::size_t index = 0;
                             index < profile.size();
                             ++index) {
                            if (index != 0) {
                                std::cout << ',';
                            }
                            std::cout << profile[index];
                        }
                        std::cout
                            << "}"
                            << " q_half=" << q_half
                            << " target=" << target
                            << " cutoff=" << cutoff
                            << " first_gap=" << gap
                            << " value=" << adjacent_payment
                            << '\n';
                    }
                }
            }
        }
    }
    return true;
}

bool enumerate_profiles(
    int length,
    int maximum_entry,
    int max_label,
    std::size_t& profiles,
    std::size_t& tails,
    std::size_t& transforms,
    std::size_t& gap_tails,
    bool& gap_counterexample,
    std::size_t& complete_wall_gap_tails,
    bool& complete_wall_gap_counterexample,
    std::size_t& gap_suffix_tails,
    bool& gap_suffix_counterexample,
    std::size_t& one_wedge_tails,
    bool& one_wedge_counterexample,
    std::size_t& adjacent_gap_payments,
    bool& adjacent_gap_counterexample
) {
    std::vector<Integer> profile(
        static_cast<std::size_t>(length),
        Integer{1}
    );
    while (true) {
        if (is_log_concave(profile)) {
            for (int shift = 0; shift <= max_label; ++shift) {
                std::vector<Integer> shifted(
                    static_cast<std::size_t>(shift),
                    Integer{0}
                );
                shifted.insert(
                    shifted.end(),
                    profile.begin(),
                    profile.end()
                );
                ++profiles;
                if (!check_profile(
                        shifted,
                        max_label,
                        tails,
                        transforms,
                        gap_tails,
                        gap_counterexample,
                        complete_wall_gap_tails,
                        complete_wall_gap_counterexample,
                        gap_suffix_tails,
                        gap_suffix_counterexample,
                        one_wedge_tails,
                        one_wedge_counterexample,
                        adjacent_gap_payments,
                        adjacent_gap_counterexample
                    )) {
                    return false;
                }
            }
        }

        int position = length - 1;
        while (position >= 0) {
            Integer& entry = profile[static_cast<std::size_t>(position)];
            if (entry < maximum_entry) {
                ++entry;
                break;
            }
            entry = 1;
            --position;
        }
        if (position < 0) {
            break;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    int max_length = 6;
    int maximum_entry = 6;
    int max_label = 5;
    try {
        if (argc >= 2) {
            max_length = parse_positive(argv[1], "max_length");
        }
        if (argc >= 3) {
            maximum_entry = parse_positive(argv[2], "maximum_entry");
        }
        if (argc >= 4) {
            max_label = parse_positive(argv[3], "max_label");
        }
        if (argc > 4) {
            throw std::invalid_argument(
                "usage: probe_su2_log_concave_wedge_tail"
                " [max_length] [maximum_entry] [max_label]"
            );
        }
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }

    std::size_t profiles = 0;
    std::size_t tails = 0;
    std::size_t transforms = 0;
    std::size_t gap_tails = 0;
    bool gap_counterexample = false;
    std::size_t complete_wall_gap_tails = 0;
    bool complete_wall_gap_counterexample = false;
    std::size_t gap_suffix_tails = 0;
    bool gap_suffix_counterexample = false;
    std::size_t one_wedge_tails = 0;
    bool one_wedge_counterexample = false;
    std::size_t adjacent_gap_payments = 0;
    bool adjacent_gap_counterexample = false;
    for (int length = 1; length <= max_length; ++length) {
        if (!enumerate_profiles(
                length,
                maximum_entry,
                max_label,
                profiles,
                tails,
                transforms,
                gap_tails,
                gap_counterexample,
                complete_wall_gap_tails,
                complete_wall_gap_counterexample,
                gap_suffix_tails,
                gap_suffix_counterexample,
                one_wedge_tails,
                one_wedge_counterexample,
                adjacent_gap_payments,
                adjacent_gap_counterexample
            )) {
            return EXIT_SUCCESS;
        }
    }
    std::cout
        << "SU2_LOG_CONCAVE_WEDGE_TAIL"
        << " profiles=" << profiles
        << " transforms=" << transforms
        << " tails=" << tails
        << " gap_tails=" << gap_tails
        << " gap_counterexample=" << (gap_counterexample ? 1 : 0)
        << " complete_wall_gap_tails=" << complete_wall_gap_tails
        << " complete_wall_gap_counterexample="
        << (complete_wall_gap_counterexample ? 1 : 0)
        << " gap_suffix_tails=" << gap_suffix_tails
        << " gap_suffix_counterexample="
        << (gap_suffix_counterexample ? 1 : 0)
        << " one_wedge_tails=" << one_wedge_tails
        << " one_wedge_counterexample="
        << (one_wedge_counterexample ? 1 : 0)
        << " adjacent_gap_payments=" << adjacent_gap_payments
        << " adjacent_gap_counterexample="
        << (adjacent_gap_counterexample ? 1 : 0)
        << " max_length=" << max_length
        << " maximum_entry=" << maximum_entry
        << " max_label=" << max_label
        << " tail_result=NO_COUNTEREXAMPLE"
        << " adjacent_gap_result="
        << (adjacent_gap_counterexample ? "COUNTEREXAMPLE" : "NO_COUNTEREXAMPLE")
        << '\n';
    return EXIT_SUCCESS;
}
