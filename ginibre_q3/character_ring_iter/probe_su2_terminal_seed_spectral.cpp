#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

struct ModePair {
    long double squared_sum = 0.0L;
    long double sum = 0.0L;
    int first = 0;
    int second = 0;
};

bool same_eigenvalue(long double first, long double second) {
    const long double scale = 1.0L
        + std::max(std::abs(first), std::abs(second));
    return std::abs(first - second) <= 1.0e-13L * scale;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: probe_su2_terminal_seed_spectral MAXIMUM_LEVEL"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        if (maximum_level < 6) {
            throw std::runtime_error("require maximum level>=6");
        }

        const long double pi = std::acos(-1.0L);
        std::size_t parameter_rows = 0;
        std::size_t spectral_groups = 0;
        std::size_t negative_group_residues = 0;
        std::size_t upper_suffix_rows = 0;
        std::size_t negative_upper_suffixes = 0;
        long double minimum_group_residue =
            std::numeric_limits<long double>::infinity();
        long double minimum_upper_suffix =
            std::numeric_limits<long double>::infinity();
        bool residue_witness_present = false;
        int witness_level = 0;
        int witness_label = 0;
        long double witness_eigenvalue = 0.0L;
        long double witness_residue = 0.0L;

        for (int level = 6; level <= maximum_level; level += 2) {
            const int width = level + 1;
            const long double normalization =
                std::sqrt(2.0L / static_cast<long double>(level + 2));
            std::vector<std::vector<long double>> transform(
                static_cast<std::size_t>(width),
                std::vector<long double>(
                    static_cast<std::size_t>(width)
                )
            );
            for (int row = 0; row <= level; ++row) {
                for (int mode = 0; mode <= level; ++mode) {
                    transform[static_cast<std::size_t>(row)]
                             [static_cast<std::size_t>(mode)] =
                        normalization
                        * std::sin(
                            static_cast<long double>(
                                (row + 1) * (mode + 1)
                            ) * pi
                            / static_cast<long double>(level + 2)
                        );
                }
            }
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameter_rows;
                std::vector<long double> fusion_eigenvalue(
                    static_cast<std::size_t>(width)
                );
                for (int mode = 0; mode <= level; ++mode) {
                    fusion_eigenvalue[static_cast<std::size_t>(mode)] =
                        transform[static_cast<std::size_t>(label)]
                                 [static_cast<std::size_t>(mode)]
                        / transform[0][static_cast<std::size_t>(mode)];
                }
                std::vector<ModePair> pairs;
                for (int first = 0; first <= level; ++first) {
                    for (int second = first + 1;
                         second <= level; ++second) {
                        const long double sum =
                            fusion_eigenvalue[
                                static_cast<std::size_t>(first)
                            ]
                            + fusion_eigenvalue[
                                static_cast<std::size_t>(second)
                            ];
                        pairs.push_back({
                            sum * sum,
                            sum,
                            first,
                            second
                        });
                    }
                }
                std::sort(
                    pairs.begin(),
                    pairs.end(),
                    [](const ModePair& first, const ModePair& second) {
                        return first.squared_sum < second.squared_sum;
                    }
                );

                std::vector<long double> residues;
                std::vector<long double> eigenvalues;
                std::size_t begin = 0;
                while (begin < pairs.size()) {
                    std::size_t end = begin + 1;
                    while (
                        end < pairs.size()
                        && same_eigenvalue(
                            pairs[begin].squared_sum,
                            pairs[end].squared_sum
                        )
                    ) {
                        ++end;
                    }
                    if (
                        std::abs(pairs[begin].squared_sum) > 1.0e-12L
                    ) {
                        long double residue = 0.0L;
                        for (std::size_t pair = begin;
                             pair < end; ++pair) {
                            const int first = pairs[pair].first;
                            const int second = pairs[pair].second;
                            const long double seed_minor =
                                transform[
                                    static_cast<std::size_t>(label)
                                ][static_cast<std::size_t>(first)]
                                * transform[0][
                                    static_cast<std::size_t>(second)
                                ]
                                - transform[
                                    static_cast<std::size_t>(label)
                                ][static_cast<std::size_t>(second)]
                                * transform[0][
                                    static_cast<std::size_t>(first)
                                ];
                            const long double target_minor =
                                transform[
                                    static_cast<std::size_t>(level)
                                ][static_cast<std::size_t>(first)]
                                * transform[0][
                                    static_cast<std::size_t>(second)
                                ]
                                - transform[
                                    static_cast<std::size_t>(level)
                                ][static_cast<std::size_t>(second)]
                                * transform[0][
                                    static_cast<std::size_t>(first)
                                ];
                            residue += pairs[pair].sum
                                * seed_minor * target_minor;
                        }
                        ++spectral_groups;
                        minimum_group_residue =
                            std::min(minimum_group_residue, residue);
                        if (residue < -1.0e-11L) {
                            ++negative_group_residues;
                            if (!residue_witness_present) {
                                residue_witness_present = true;
                                witness_level = level;
                                witness_label = label;
                                witness_eigenvalue =
                                    pairs[begin].squared_sum;
                                witness_residue = residue;
                            }
                        }
                        eigenvalues.push_back(
                            pairs[begin].squared_sum
                        );
                        residues.push_back(residue);
                    }
                    begin = end;
                }

                long double suffix = 0.0L;
                for (std::size_t group = residues.size();
                     group > 0; --group) {
                    suffix += residues[group - 1];
                    ++upper_suffix_rows;
                    minimum_upper_suffix =
                        std::min(minimum_upper_suffix, suffix);
                    if (suffix < -1.0e-10L) {
                        ++negative_upper_suffixes;
                        std::cout
                            << "SU2_TERMINAL_SEED_SPECTRAL"
                            << " suffix_result=FAIL_NUMERICAL"
                            << " level=" << level
                            << " label=" << label
                            << " cutoff_eigenvalue="
                            << static_cast<double>(
                                eigenvalues[group - 1]
                            )
                            << " suffix="
                            << static_cast<double>(suffix)
                            << " negative_group_residues="
                            << negative_group_residues;
                        if (residue_witness_present) {
                            std::cout
                                << " first_negative_residue={level="
                                << witness_level
                                << ",label=" << witness_label
                                << ",eigenvalue="
                                << static_cast<double>(
                                    witness_eigenvalue
                                )
                                << ",residue="
                                << static_cast<double>(
                                    witness_residue
                                )
                                << '}';
                        }
                        std::cout
                            << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_SEED_SPECTRAL"
            << " maximum_level=" << maximum_level
            << " parameter_rows=" << parameter_rows
            << " spectral_groups=" << spectral_groups
            << " negative_group_residues="
            << negative_group_residues
            << " minimum_group_residue="
            << static_cast<double>(minimum_group_residue)
            << " upper_suffix_rows=" << upper_suffix_rows
            << " negative_upper_suffixes="
            << negative_upper_suffixes
            << " minimum_upper_suffix="
            << static_cast<double>(minimum_upper_suffix);
        if (residue_witness_present) {
            std::cout
                << " first_negative_residue={level=" << witness_level
                << ",label=" << witness_label
                << ",eigenvalue="
                << static_cast<double>(witness_eigenvalue)
                << ",residue="
                << static_cast<double>(witness_residue)
                << '}';
        }
        std::cout
            << " residue_cone="
            << (
                negative_group_residues == 0
                    ? "PASS_NUMERICAL_DISCOVERY" : "FAIL"
            )
            << " suffix_result=PASS_NUMERICAL_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_SEED_SPECTRAL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
