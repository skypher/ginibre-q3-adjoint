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
    long double eigenvalue = 0.0L;
    int first = 0;
    int second = 0;
};

struct SpectralGroup {
    long double eigenvalue = 0.0L;
    std::size_t multiplicity = 0;
    std::vector<long double> residue;
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
                "usage: probe_su2_terminal_boundary_spectral "
                "MAXIMUM_LEVEL"
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
        std::size_t residue_entries = 0;
        std::size_t negative_residue_entries = 0;
        std::size_t spectral_suffix_entries = 0;
        std::size_t negative_spectral_suffixes = 0;
        long double minimum_residue =
            std::numeric_limits<long double>::infinity();
        long double minimum_spectral_suffix =
            std::numeric_limits<long double>::infinity();
        bool residue_witness_present = false;
        int residue_witness_level = 0;
        int residue_witness_label = 0;
        int residue_witness_source = 0;
        int residue_witness_target = 0;
        long double residue_witness_eigenvalue = 0.0L;
        std::size_t residue_witness_multiplicity = 0;
        long double residue_witness_value = 0.0L;
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
            for (int label = 0; label <= level; ++label) {
                for (int mode = 0; mode <= level; ++mode) {
                    transform[static_cast<std::size_t>(label)]
                             [static_cast<std::size_t>(mode)] =
                        normalization
                        * std::sin(
                            static_cast<long double>(
                                (label + 1) * (mode + 1)
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
                            first,
                            second
                        });
                    }
                }
                std::sort(
                    pairs.begin(),
                    pairs.end(),
                    [](const ModePair& first, const ModePair& second) {
                        return first.eigenvalue < second.eigenvalue;
                    }
                );

                std::vector<SpectralGroup> groups;
                std::size_t begin = 0;
                while (begin < pairs.size()) {
                    std::size_t end = begin + 1;
                    while (
                        end < pairs.size()
                        && same_eigenvalue(
                            pairs[begin].eigenvalue,
                            pairs[end].eigenvalue
                        )
                    ) {
                        ++end;
                    }
                    ++spectral_groups;
                    if (
                        std::abs(pairs[begin].eigenvalue) <= 1.0e-12L
                    ) {
                        begin = end;
                        continue;
                    }
                    std::vector<long double> residue(
                        static_cast<std::size_t>(width * width)
                    );
                    for (std::size_t pair = begin;
                         pair < end; ++pair) {
                        std::vector<long double> boundary(
                            static_cast<std::size_t>(width)
                        );
                        const int first = pairs[pair].first;
                        const int second = pairs[pair].second;
                        for (int source = 1;
                             source <= level; ++source) {
                            boundary[
                                static_cast<std::size_t>(source)
                            ] =
                                transform[
                                    static_cast<std::size_t>(source)
                                ][static_cast<std::size_t>(first)]
                                * transform[0][
                                    static_cast<std::size_t>(second)
                                ]
                                - transform[
                                    static_cast<std::size_t>(source)
                                ][static_cast<std::size_t>(second)]
                                * transform[0][
                                    static_cast<std::size_t>(first)
                                ];
                        }
                        for (int source = 1;
                             source <= level; ++source) {
                            for (int target = 1;
                                 target <= level; ++target) {
                                residue[
                                    static_cast<std::size_t>(
                                        source * width + target
                                    )
                                ] += boundary[
                                    static_cast<std::size_t>(source)
                                ] * boundary[
                                    static_cast<std::size_t>(target)
                                ];
                            }
                        }
                    }
                    for (int source = 1;
                         source <= level; ++source) {
                        for (int target = 1;
                             target <= level; ++target) {
                            if ((source & 1) != (target & 1)) {
                                continue;
                            }
                            ++residue_entries;
                            const long double value = residue[
                                static_cast<std::size_t>(
                                    source * width + target
                                )
                            ];
                            minimum_residue =
                                std::min(minimum_residue, value);
                            if (value < -1.0e-11L) {
                                ++negative_residue_entries;
                                if (!residue_witness_present) {
                                    residue_witness_present = true;
                                    residue_witness_level = level;
                                    residue_witness_label = label;
                                    residue_witness_source = source;
                                    residue_witness_target = target;
                                    residue_witness_eigenvalue =
                                        pairs[begin].eigenvalue;
                                    residue_witness_multiplicity =
                                        end - begin;
                                    residue_witness_value = value;
                                }
                            }
                        }
                    }
                    groups.push_back({
                        pairs[begin].eigenvalue,
                        end - begin,
                        std::move(residue)
                    });
                    begin = end;
                }
                for (int source = 1;
                     source <= level; ++source) {
                    for (int target = 1;
                         target <= level; ++target) {
                        if ((source & 1) != (target & 1)) {
                            continue;
                        }
                        long double suffix = 0.0L;
                        for (std::size_t group = groups.size();
                             group > 0; --group) {
                            const SpectralGroup& current =
                                groups[group - 1];
                            suffix += current.residue[
                                static_cast<std::size_t>(
                                    source * width + target
                                )
                            ];
                            ++spectral_suffix_entries;
                            minimum_spectral_suffix = std::min(
                                minimum_spectral_suffix,
                                suffix
                            );
                            if (suffix < -1.0e-10L) {
                                ++negative_spectral_suffixes;
                                std::cout
                                    << "SU2_TERMINAL_BOUNDARY_SPECTRAL"
                                    << " suffix_result=FAIL_NUMERICAL"
                                    << " level=" << level
                                    << " label=" << label
                                    << " cutoff_eigenvalue="
                                    << static_cast<double>(
                                        current.eigenvalue
                                    )
                                    << " source=" << source
                                    << " target=" << target
                                    << " suffix="
                                    << static_cast<double>(suffix)
                                    << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_BOUNDARY_SPECTRAL"
            << " maximum_level=" << maximum_level
            << " parameter_rows=" << parameter_rows
            << " spectral_groups=" << spectral_groups
            << " residue_entries=" << residue_entries
            << " negative_residue_entries="
            << negative_residue_entries
            << " minimum_residue="
            << static_cast<double>(minimum_residue)
            << " spectral_suffix_entries="
            << spectral_suffix_entries
            << " negative_spectral_suffixes="
            << negative_spectral_suffixes
            << " minimum_spectral_suffix="
            << static_cast<double>(minimum_spectral_suffix);
        if (residue_witness_present) {
            std::cout
                << " first_negative_residue={level="
                << residue_witness_level
                << ",label=" << residue_witness_label
                << ",eigenvalue="
                << static_cast<double>(residue_witness_eigenvalue)
                << ",multiplicity="
                << residue_witness_multiplicity
                << ",source=" << residue_witness_source
                << ",target=" << residue_witness_target
                << ",residue="
                << static_cast<double>(residue_witness_value)
                << '}';
        }
        std::cout
            << " residue_cone="
            << (
                negative_residue_entries == 0
                    ? "PASS_NUMERICAL_DISCOVERY" : "FAIL"
            )
            << " suffix_result=PASS_NUMERICAL_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_BOUNDARY_SPECTRAL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
