#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

struct SignedBlock {
    int label = 0;
    int sign = 1;
    int multiplicity = 0;
};

struct Witness {
    bool initialized = false;
    cpp_int value = 0;
    int target = 0;
    std::vector<SignedBlock> blocks;
};

struct Statistics {
    std::uint64_t words = 0U;
    std::uint64_t corner_entries = 0U;
    std::uint64_t boundary_entries = 0U;
    std::uint64_t zero_corners = 0U;
    std::uint64_t zero_boundary_entries = 0U;
    std::uint64_t cyclic_entries = 0U;
    Witness minimum_corner;
    Witness minimum_boundary;
    Witness minimum_positive_corner;
    Witness minimum_positive_boundary;
};

using LinearForm = std::map<int, long long>;
using QuadraticForm = std::map<std::pair<int, int>, long long>;

void add_cyclic_charge(
    LinearForm& form,
    int rank,
    int parity,
    int charge,
    long long coefficient
) {
    const int order = 2 * rank + 1;
    int residue = charge % order;
    if (residue < 0) {
        residue += order;
    }
    if (residue > rank) {
        residue = order - residue;
        coefficient *= static_cast<long long>(parity);
    }
    if (residue == 0 && parity == -1) {
        return;
    }
    form[residue] += coefficient;
}

LinearForm multiplied_cyclic_coefficient(
    int rank,
    int parity,
    int output,
    int feature,
    int local_sign
) {
    LinearForm result;
    add_cyclic_charge(result, rank, parity, output - feature, 1);
    add_cyclic_charge(
        result,
        rank,
        parity,
        output + feature,
        static_cast<long long>(local_sign)
    );
    return result;
}

LinearForm cyclic_unit(int rank, int parity, int charge) {
    LinearForm result;
    add_cyclic_charge(result, rank, parity, charge, 1);
    return result;
}

void add_form_product(
    QuadraticForm& destination,
    const LinearForm& left,
    const LinearForm& right,
    long long scale
) {
    for (const auto& [left_charge, left_coefficient] : left) {
        for (const auto& [right_charge, right_coefficient] : right) {
            destination[std::minmax(left_charge, right_charge)]
                += scale * left_coefficient * right_coefficient;
        }
    }
}

void add_three_charge_q(
    QuadraticForm& destination,
    const LinearForm& zero,
    const LinearForm& one,
    const LinearForm& two,
    long long scale
) {
    add_form_product(destination, zero, zero, scale);
    add_form_product(destination, one, one, -2 * scale);
    add_form_product(destination, zero, two, scale);
}

void erase_zero_terms(QuadraticForm& form) {
    for (auto iterator = form.begin(); iterator != form.end();) {
        if (iterator->second == 0) {
            iterator = form.erase(iterator);
        } else {
            ++iterator;
        }
    }
}

QuadraticForm local_cyclic_q(int rank, int radius, int sign) {
    QuadraticForm result;
    if (sign == 1) {
        add_three_charge_q(
            result,
            cyclic_unit(rank, sign, 0),
            cyclic_unit(rank, sign, 1),
            cyclic_unit(rank, sign, 2),
            2
        );
    }
    for (int feature = 1; feature <= radius; ++feature) {
        add_three_charge_q(
            result,
            multiplied_cyclic_coefficient(
                rank, sign, 0, feature, sign
            ),
            multiplied_cyclic_coefficient(
                rank, sign, 1, feature, sign
            ),
            multiplied_cyclic_coefficient(
                rank, sign, 2, feature, sign
            ),
            1
        );
    }
    erase_zero_terms(result);
    return result;
}

QuadraticForm expected_cyclic_turan(int rank, int radius, int sign) {
    QuadraticForm result;
    const auto add_h = [&](int charge, long long scale) {
        const LinearForm center = cyclic_unit(rank, sign, charge);
        const LinearForm lower = cyclic_unit(rank, sign, charge - 1);
        const LinearForm upper = cyclic_unit(rank, sign, charge + 1);
        add_form_product(result, center, center, scale);
        add_form_product(result, lower, upper, -scale);
    };
    add_h(radius, 2);
    add_h(radius + 1, -2);
    erase_zero_terms(result);
    return result;
}

bool run_turan_symbolic_check(int maximum_rank) {
    std::uint64_t identities = 0U;
    for (int rank = 2; rank <= maximum_rank; ++rank) {
        for (int radius = 1; radius < rank; ++radius) {
            for (int sign : {1, -1}) {
                const QuadraticForm actual
                    = local_cyclic_q(rank, radius, sign);
                const QuadraticForm expected
                    = expected_cyclic_turan(rank, radius, sign);
                ++identities;
                if (actual != expected) {
                    std::cout
                        << "SU2_ODD_ORBIT_TURAN result=FAIL rank=" << rank
                        << " radius=" << radius << " sign=" << sign
                        << " actual_terms=" << actual.size()
                        << " expected_terms=" << expected.size() << '\n';
                    return false;
                }
            }
        }
    }
    std::cout << "SU2_ODD_ORBIT_TURAN ranks=2.." << maximum_rank
              << " identities=" << identities << " result=PASS\n";
    return true;
}

int even_lift_index(int rank, int orbit_label) {
    const int level = 2 * rank - 1;
    const int even_label = (orbit_label & 1) == 0
        ? orbit_label
        : level - orbit_label;
    return even_label / 2;
}

template <class Function>
void for_each_orbit_output(int rank, int left, int right, Function function) {
    const int level = 2 * rank - 1;
    const int lower = std::abs(left - right);
    const int upper = std::min(left + right, 2 * level - left - right);
    for (int output = lower; output <= upper; output += 2) {
        function(std::min(output, level - output));
    }
}

std::vector<cpp_int> apply_factor(
    int rank,
    const std::vector<cpp_int>& current,
    int label,
    int sign
) {
    if (rank < 1 || label < 1 || label >= rank
        || (sign != -1 && sign != 1)) {
        throw std::runtime_error("invalid orbit factor");
    }
    const std::size_t expected
        = static_cast<std::size_t>(rank * rank);
    if (current.size() != expected) {
        throw std::runtime_error("orbit state has wrong size");
    }

    std::vector<cpp_int> next(expected, 0);
    for (int left = 0; left < rank; ++left) {
        for (int right = 0; right < rank; ++right) {
            const cpp_int& value = current[static_cast<std::size_t>(
                left * rank + right
            )];
            if (value == 0) {
                continue;
            }
            for_each_orbit_output(rank, left, label, [&](int output) {
                next[static_cast<std::size_t>(
                    output * rank + right
                )] += value;
            });
            for_each_orbit_output(rank, right, label, [&](int output) {
                cpp_int& destination = next[static_cast<std::size_t>(
                    left * rank + output
                )];
                destination += sign == 1 ? value : -value;
            });
        }
    }
    return next;
}

bool run_rank_three_depth(int maximum_power) {
    if (maximum_power < 0) {
        throw std::runtime_error("invalid rank-three depth");
    }
    constexpr int rank = 3;
    constexpr int first_even_lift_orbit_label = 2;
    constexpr int second_even_lift_orbit_label = 1;
    std::uint64_t words = 0U;
    std::uint64_t boundary_entries = 0U;
    cpp_int minimum_corner = 0;
    cpp_int minimum_boundary = 0;
    bool initialized = false;
    for (int first_sign : {1, -1}) {
        for (int second_sign : {1, -1}) {
            std::vector<cpp_int> first_state(9U, 0);
            first_state[0] = 1;
            for (int first_power = 0;
                 first_power <= maximum_power;
                 ++first_power) {
                std::vector<cpp_int> state = first_state;
                for (int second_power = 0;
                     second_power <= maximum_power;
                     ++second_power) {
                    const int exchange_sign
                        = (first_sign == -1 && (first_power & 1) != 0)
                            != (second_sign == -1
                                && (second_power & 1) != 0)
                        ? -1
                        : 1;
                    const cpp_int& corner = state[0];
                    if (!initialized || corner < minimum_corner) {
                        minimum_corner = corner;
                    }
                    if (corner < 0
                        || (exchange_sign == -1 && corner != 0)) {
                        std::cout
                            << "SU2_ODD_ORBIT_RANK3_DEPTH result=FAIL_CORNER"
                            << " first_sign=" << first_sign
                            << " first_power=" << first_power
                            << " second_sign=" << second_sign
                            << " second_power=" << second_power
                            << " corner=" << corner << '\n';
                        return false;
                    }
                    for (int target = 0; target < rank; ++target) {
                        const cpp_int& value = state[
                            static_cast<std::size_t>(target * rank)
                        ];
                        if (!initialized || value < minimum_boundary) {
                            minimum_boundary = value;
                        }
                        ++boundary_entries;
                        if (value < 0) {
                            std::cout
                                << "SU2_ODD_ORBIT_RANK3_DEPTH "
                                << "result=FAIL_BOUNDARY"
                                << " first_sign=" << first_sign
                                << " first_power=" << first_power
                                << " second_sign=" << second_sign
                                << " second_power=" << second_power
                                << " target=" << target
                                << " value=" << value << '\n';
                            return false;
                        }
                    }
                    initialized = true;
                    ++words;
                    if (second_power != maximum_power) {
                        state = apply_factor(
                            rank,
                            state,
                            second_even_lift_orbit_label,
                            second_sign
                        );
                    }
                }
                if (first_power != maximum_power) {
                    first_state = apply_factor(
                        rank,
                        first_state,
                        first_even_lift_orbit_label,
                        first_sign
                    );
                }
            }
        }
    }
    std::cout << "SU2_ODD_ORBIT_RANK3_DEPTH maximum_power="
              << maximum_power << " words=" << words
              << " boundary_entries=" << boundary_entries
              << " minimum_corner=" << minimum_corner
              << " minimum_boundary=" << minimum_boundary
              << " result=PASS\n";
    return true;
}

void print_word(const std::vector<SignedBlock>& blocks) {
    std::cout << '[';
    bool first = true;
    for (const SignedBlock& block : blocks) {
        for (int copy = 0; copy < block.multiplicity; ++copy) {
            if (!first) {
                std::cout << ',';
            }
            std::cout << (block.sign == 1 ? '+' : '-') << block.label;
            first = false;
        }
    }
    std::cout << ']';
}

void consider(
    Witness& witness,
    const cpp_int& value,
    int target,
    const std::vector<SignedBlock>& blocks
) {
    if (!witness.initialized || value < witness.value) {
        witness.initialized = true;
        witness.value = value;
        witness.target = target;
        witness.blocks = blocks;
    }
}

std::vector<cpp_int> cyclic_product(
    int rank,
    const std::vector<SignedBlock>& blocks
) {
    const int order = 2 * rank + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(order * order), 0
    );
    current[0] = 1;
    for (const SignedBlock& block : blocks) {
        const int radius = even_lift_index(rank, block.label);
        for (int copy = 0; copy < block.multiplicity; ++copy) {
            std::vector<cpp_int> next(current.size(), 0);
            for (int left = 0; left < order; ++left) {
                for (int right = 0; right < order; ++right) {
                    const cpp_int& value = current[static_cast<std::size_t>(
                        left * order + right
                    )];
                    if (value == 0) {
                        continue;
                    }
                    for (int step = -radius; step <= radius; ++step) {
                        const int left_target = (left + step + order) % order;
                        next[static_cast<std::size_t>(
                            left_target * order + right
                        )] += value;
                        const int right_target
                            = (right + step + order) % order;
                        cpp_int& destination = next[static_cast<std::size_t>(
                            left * order + right_target
                        )];
                        destination += block.sign == 1 ? value : -value;
                    }
                }
            }
            current = std::move(next);
        }
    }
    return current;
}

cpp_int cyclic_monge_coefficient(
    int rank,
    int orbit_target,
    const std::vector<cpp_int>& cyclic
) {
    const int order = 2 * rank + 1;
    const int target = even_lift_index(rank, orbit_target);
    const int next_target = target + 1;
    const auto entry = [&](int left, int right) -> const cpp_int& {
        return cyclic[static_cast<std::size_t>(left * order + right)];
    };
    return entry(target, 0) - entry(target, 1)
        - entry(next_target, 0) + entry(next_target, 1);
}

bool inspect_state(
    int rank,
    const std::vector<cpp_int>& state,
    const std::vector<SignedBlock>& blocks,
    bool check_cyclic_model,
    Statistics& statistics
) {
    int exchange_sign = 1;
    for (const SignedBlock& block : blocks) {
        if (block.sign == -1 && (block.multiplicity & 1) != 0) {
            exchange_sign = -exchange_sign;
        }
    }
    for (int left = 0; left < rank; ++left) {
        for (int right = 0; right < rank; ++right) {
            const cpp_int& value = state[static_cast<std::size_t>(
                left * rank + right
            )];
            const cpp_int& transpose = state[static_cast<std::size_t>(
                right * rank + left
            )];
            if (value != exchange_sign * transpose) {
                std::cout << "SU2_ODD_ORBIT_GKS result=FAIL_EXCHANGE"
                          << " rank=" << rank << " word=";
                print_word(blocks);
                std::cout << " state=(" << left << ',' << right << ')'
                          << " value=" << value
                          << " transpose=" << transpose
                          << " exchange_sign=" << exchange_sign << '\n';
                return false;
            }
        }
    }

    const cpp_int& corner = state[0];
    ++statistics.words;
    ++statistics.corner_entries;
    consider(statistics.minimum_corner, corner, 0, blocks);
    if (corner == 0) {
        ++statistics.zero_corners;
    } else if (corner > 0) {
        consider(statistics.minimum_positive_corner, corner, 0, blocks);
    }
    if (corner < 0) {
        std::cout << "SU2_ODD_ORBIT_GKS result=FAIL_GKS2_STAR"
                  << " rank=" << rank << " level=" << 2 * rank - 1
                  << " word=";
        print_word(blocks);
        std::cout << " corner=" << corner << '\n';
        return false;
    }
    if (exchange_sign == -1 && corner != 0) {
        std::cout << "SU2_ODD_ORBIT_GKS result=FAIL_ODD_CORNER"
                  << " rank=" << rank << " word=";
        print_word(blocks);
        std::cout << " corner=" << corner << '\n';
        return false;
    }

    std::vector<cpp_int> cyclic;
    if (check_cyclic_model) {
        cyclic = cyclic_product(rank, blocks);
    }
    for (int target = 0; target < rank; ++target) {
        const cpp_int& value = state[static_cast<std::size_t>(target * rank)];
        if (check_cyclic_model) {
            const cpp_int cyclic_value = cyclic_monge_coefficient(
                rank, target, cyclic
            );
            ++statistics.cyclic_entries;
            if (cyclic_value != value) {
                std::cout << "SU2_ODD_ORBIT_GKS result=FAIL_CYCLIC_MONGE"
                          << " rank=" << rank << " word=";
                print_word(blocks);
                std::cout << " target=" << target
                          << " orbit=" << value
                          << " cyclic=" << cyclic_value << '\n';
                return false;
            }
        }
        ++statistics.boundary_entries;
        consider(statistics.minimum_boundary, value, target, blocks);
        if (value == 0) {
            ++statistics.zero_boundary_entries;
        } else if (value > 0) {
            consider(
                statistics.minimum_positive_boundary,
                value,
                target,
                blocks
            );
        }
        if (value < 0) {
            std::cout << "SU2_ODD_ORBIT_GKS result=FAIL_PARTIAL_CHARACTER"
                      << " rank=" << rank << " level=" << 2 * rank - 1
                      << " word=";
            print_word(blocks);
            std::cout << " target=" << target << " value=" << value
                      << '\n';
            return false;
        }
    }
    return true;
}

bool enumerate_exact_degree(
    int rank,
    int label,
    int remaining,
    const std::vector<cpp_int>& state,
    std::vector<SignedBlock>& blocks,
    bool check_cyclic_model,
    Statistics& statistics
) {
    if (label == rank) {
        return remaining != 0
            || inspect_state(
                rank,
                state,
                blocks,
                check_cyclic_model,
                statistics
            );
    }

    if (!enumerate_exact_degree(
            rank,
            label + 1,
            remaining,
            state,
            blocks,
            check_cyclic_model,
            statistics
        )) {
        return false;
    }

    for (int sign : {1, -1}) {
        std::vector<cpp_int> updated = state;
        for (int multiplicity = 1;
             multiplicity <= remaining;
             ++multiplicity) {
            updated = apply_factor(rank, updated, label, sign);
            blocks.push_back(SignedBlock{label, sign, multiplicity});
            const bool passed = enumerate_exact_degree(
                rank,
                label + 1,
                remaining - multiplicity,
                updated,
                blocks,
                check_cyclic_model,
                statistics
            );
            blocks.pop_back();
            if (!passed) {
                return false;
            }
        }
    }
    return true;
}

void print_witness(const char* name, const Witness& witness) {
    std::cout << ' ' << name << '=';
    if (!witness.initialized) {
        std::cout << "none";
        return;
    }
    std::cout << witness.value << "@target" << witness.target << ':';
    print_word(witness.blocks);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "rank3-depth") {
            const int maximum_power = std::stoi(argv[2]);
            return run_rank_three_depth(maximum_power)
                ? EXIT_SUCCESS
                : EXIT_FAILURE;
        }
        if (argc == 3 && std::string(argv[1]) == "turan-check") {
            const int maximum_rank = std::stoi(argv[2]);
            if (maximum_rank < 2) {
                throw std::runtime_error("invalid Turan-check rank");
            }
            return run_turan_symbolic_check(maximum_rank)
                ? EXIT_SUCCESS
                : EXIT_FAILURE;
        }
        const bool check_cyclic_model
            = argc == 4 && std::string(argv[1]) == "cyclic-check";
        if (argc != 3 && !check_cyclic_model) {
            std::cerr << "usage: " << argv[0]
                      << " [cyclic-check] MAXIMUM_RANK MAXIMUM_FACTORS\n"
                      << "       " << argv[0]
                      << " turan-check MAXIMUM_RANK\n"
                      << "       " << argv[0]
                      << " rank3-depth MAXIMUM_POWER\n";
            return EXIT_FAILURE;
        }
        const int argument_offset = check_cyclic_model ? 1 : 0;
        const int maximum_rank = std::stoi(argv[1 + argument_offset]);
        const int maximum_factors = std::stoi(argv[2 + argument_offset]);
        if (maximum_rank < 2 || maximum_factors < 0) {
            throw std::runtime_error("invalid search bounds");
        }

        Statistics total;
        for (int rank = 2; rank <= maximum_rank; ++rank) {
            Statistics rank_statistics;
            for (int factors = 0; factors <= maximum_factors; ++factors) {
                std::vector<cpp_int> state(
                    static_cast<std::size_t>(rank * rank), 0
                );
                state[0] = 1;
                std::vector<SignedBlock> blocks;
                const std::uint64_t words_before = rank_statistics.words;
                if (!enumerate_exact_degree(
                        rank,
                        1,
                        factors,
                        state,
                        blocks,
                        check_cyclic_model,
                        rank_statistics
                    )) {
                    return EXIT_FAILURE;
                }
                std::cout << "progress rank=" << rank
                          << " level=" << 2 * rank - 1
                          << " factors=" << factors
                          << " degree_words="
                          << rank_statistics.words - words_before
                          << " cumulative_words=" << rank_statistics.words
                          << " boundary_entries="
                          << rank_statistics.boundary_entries
                          << " cyclic_entries="
                          << rank_statistics.cyclic_entries << '\n';
            }
            total.words += rank_statistics.words;
            total.corner_entries += rank_statistics.corner_entries;
            total.boundary_entries += rank_statistics.boundary_entries;
            total.zero_corners += rank_statistics.zero_corners;
            total.zero_boundary_entries
                += rank_statistics.zero_boundary_entries;
            total.cyclic_entries += rank_statistics.cyclic_entries;
            if (rank_statistics.minimum_corner.initialized) {
                consider(
                    total.minimum_corner,
                    rank_statistics.minimum_corner.value,
                    rank_statistics.minimum_corner.target,
                    rank_statistics.minimum_corner.blocks
                );
            }
            if (rank_statistics.minimum_boundary.initialized) {
                consider(
                    total.minimum_boundary,
                    rank_statistics.minimum_boundary.value,
                    rank_statistics.minimum_boundary.target,
                    rank_statistics.minimum_boundary.blocks
                );
            }
            if (rank_statistics.minimum_positive_corner.initialized) {
                consider(
                    total.minimum_positive_corner,
                    rank_statistics.minimum_positive_corner.value,
                    rank_statistics.minimum_positive_corner.target,
                    rank_statistics.minimum_positive_corner.blocks
                );
            }
            if (rank_statistics.minimum_positive_boundary.initialized) {
                consider(
                    total.minimum_positive_boundary,
                    rank_statistics.minimum_positive_boundary.value,
                    rank_statistics.minimum_positive_boundary.target,
                    rank_statistics.minimum_positive_boundary.blocks
                );
            }
        }

        std::cout << "SU2_ODD_ORBIT_GKS ranks=2.." << maximum_rank
                  << " maximum_factors=" << maximum_factors
                  << " words=" << total.words
                  << " corner_entries=" << total.corner_entries
                  << " boundary_entries=" << total.boundary_entries
                  << " zero_corners=" << total.zero_corners
                  << " zero_boundary_entries="
                  << total.zero_boundary_entries
                  << " cyclic_entries=" << total.cyclic_entries;
        print_witness("minimum_corner", total.minimum_corner);
        print_witness(
            "minimum_positive_corner",
            total.minimum_positive_corner
        );
        print_witness("minimum_boundary", total.minimum_boundary);
        print_witness(
            "minimum_positive_boundary",
            total.minimum_positive_boundary
        );
        std::cout << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
