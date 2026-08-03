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

using Vector = std::vector<int>;

int parse_nonnegative(const char* text, const char* name) {
    const std::string input(text);
    std::size_t used = 0U;
    const long long value = std::stoll(input, &used, 10);
    if (used != input.size() || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

void print_vector(const Vector& values) {
    std::cout << '[';
    for (std::size_t i = 0U; i < values.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << values[i];
    }
    std::cout << ']';
}

bool finite_valid_top(const Vector& degrees, const Vector& top, int level) {
    if (degrees.size() != top.size()) {
        return false;
    }
    int height = 0;
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        const int alpha = top[i];
        const int beta = degrees[i] - alpha;
        if (degrees[i] < 0 || degrees[i] > level || alpha < 0
            || beta < 0 || beta > height || alpha > level - height) {
            return false;
        }
        height -= beta;
        height += alpha;
    }
    return height == 0;
}

void enumerate_tops(const Vector& degrees, int level, std::size_t index,
                    int height, Vector& top, std::vector<Vector>& output) {
    if (index == degrees.size()) {
        if (height == 0) {
            output.push_back(top);
        }
        return;
    }
    const int degree = degrees[index];
    const int lower = std::max(0, degree - height);
    const int upper = std::min(degree, level - height);
    for (int alpha = lower; alpha <= upper; ++alpha) {
        const int beta = degree - alpha;
        top[index] = alpha;
        enumerate_tops(
            degrees, level, index + 1U, height - beta + alpha, top, output
        );
    }
}

std::vector<Vector> finite_tops(const Vector& degrees, int level) {
    Vector top(degrees.size(), 0);
    std::vector<Vector> output;
    enumerate_tops(degrees, level, 0U, 0, top, output);
    return output;
}

Vector bender_knuth_swap(const Vector& degrees, const Vector& top,
                         std::size_t position) {
    if (position + 1U >= degrees.size()) {
        throw std::runtime_error("swap position out of range");
    }
    int height = 0;
    for (std::size_t i = 0U; i < position; ++i) {
        height += 2 * top[i] - degrees[i];
    }
    const int beta_left = degrees[position] - top[position];
    const int beta_right = degrees[position + 1U] - top[position + 1U];
    const int transfer = std::max(0, beta_left + beta_right - height);
    Vector result = top;
    result[position] = top[position + 1U] + transfer;
    result[position + 1U] = top[position] - transfer;
    return result;
}

Vector fusion_outputs(int left, int label, int level) {
    const int lower = std::abs(left - label);
    const int upper = std::min(left + label, 2 * level - left - label);
    Vector output;
    for (int value = lower; value <= upper; value += 2) {
        output.push_back(value);
    }
    return output;
}

Vector common_intermediates(int left, int first, int right, int second,
                            int level) {
    const Vector first_outputs = fusion_outputs(left, first, level);
    const Vector reverse_outputs = fusion_outputs(right, second, level);
    Vector output;
    std::set_intersection(
        first_outputs.begin(), first_outputs.end(), reverse_outputs.begin(),
        reverse_outputs.end(), std::back_inserter(output)
    );
    return output;
}

Vector finite_recoupling_swap(const Vector& degrees, const Vector& top,
                              std::size_t position, int level) {
    if (position + 1U >= degrees.size()) {
        throw std::runtime_error("swap position out of range");
    }
    Vector heights(degrees.size() + 1U, 0);
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        heights[i + 1U] = heights[i] + 2 * top[i] - degrees[i];
    }
    const int left = heights[position];
    const int middle = heights[position + 1U];
    const int right = heights[position + 2U];
    const int first = degrees[position];
    const int second = degrees[position + 1U];
    const Vector old_middle = common_intermediates(
        left, first, right, second, level
    );
    const Vector new_middle = common_intermediates(
        left, second, right, first, level
    );
    if (old_middle.size() != new_middle.size()) {
        throw std::runtime_error("finite associativity multiplicity mismatch");
    }
    const auto found = std::find(old_middle.begin(), old_middle.end(), middle);
    if (found == old_middle.end()) {
        throw std::runtime_error("finite path has invalid local intermediate");
    }
    const std::size_t rank = static_cast<std::size_t>(
        std::distance(old_middle.begin(), found)
    );
    heights[position + 1U] = new_middle[rank];
    Vector swapped_degrees = degrees;
    std::swap(swapped_degrees[position], swapped_degrees[position + 1U]);
    Vector result(degrees.size(), 0);
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        const int numerator
            = swapped_degrees[i] + heights[i + 1U] - heights[i];
        if ((numerator & 1) != 0) {
            throw std::runtime_error("recoupled top is nonintegral");
        }
        result[i] = numerator / 2;
    }
    if (!finite_valid_top(swapped_degrees, result, level)) {
        throw std::runtime_error("recoupling did not preserve finite path");
    }
    return result;
}

struct CrystalToken {
    int label = 0;
    int down = 0;
};

void crystal_r_swap(CrystalToken& left, CrystalToken& right) {
    const int channel = std::min(right.down, left.label - left.down);
    const int distance = left.down + right.down - channel;
    const int new_left_down = std::min(distance, right.label - channel);
    const int new_right_down = channel
        + std::max(0, distance - (right.label - channel));
    const CrystalToken new_left{right.label, new_left_down};
    const CrystalToken new_right{left.label, new_right_down};
    left = new_left;
    right = new_right;
}

Vector finite_crystal_swap(const Vector& degrees, const Vector& top,
                           std::size_t position) {
    if (position + 1U >= degrees.size()) {
        throw std::runtime_error("swap position out of range");
    }
    CrystalToken first{degrees[position], degrees[position] - top[position]};
    CrystalToken second{
        degrees[position + 1U], degrees[position + 1U] - top[position + 1U]
    };
    crystal_r_swap(first, second);
    Vector swapped_degrees = degrees;
    std::swap(swapped_degrees[position], swapped_degrees[position + 1U]);
    Vector result = top;
    result[position] = swapped_degrees[position] - first.down;
    result[position + 1U]
        = swapped_degrees[position + 1U] - second.down;
    return result;
}

struct FusionPath {
    Vector degrees;
    Vector top;
};

FusionPath apply_recoupling_swap(
    const FusionPath& path, std::size_t position, int level
) {
    FusionPath result;
    result.degrees = path.degrees;
    std::swap(result.degrees[position], result.degrees[position + 1U]);
    result.top = finite_recoupling_swap(path.degrees, path.top, position, level);
    return result;
}

bool same_path(const FusionPath& left, const FusionPath& right) {
    return left.degrees == right.degrees && left.top == right.top;
}

bool same_tokens(
    const std::vector<CrystalToken>& left,
    const std::vector<CrystalToken>& right
) {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t index = 0U; index < left.size(); ++index) {
        if (left[index].label != right[index].label
            || left[index].down != right[index].down) {
            return false;
        }
    }
    return true;
}

void crystal_r_check(int maximum_label) {
    std::uint64_t involutions = 0U;
    std::uint64_t braid_tests = 0U;
    for (int first_label = 0; first_label <= maximum_label; ++first_label) {
        for (int second_label = 0; second_label <= maximum_label;
             ++second_label) {
            for (int first_down = 0; first_down <= first_label; ++first_down) {
                for (int second_down = 0; second_down <= second_label;
                     ++second_down) {
                    CrystalToken left{first_label, first_down};
                    CrystalToken right{second_label, second_down};
                    const CrystalToken original_left = left;
                    const CrystalToken original_right = right;
                    crystal_r_swap(left, right);
                    crystal_r_swap(left, right);
                    ++involutions;
                    if (left.label != original_left.label
                        || left.down != original_left.down
                        || right.label != original_right.label
                        || right.down != original_right.down) {
                        throw std::runtime_error("crystal R is not involutive");
                    }
                }
            }
        }
    }
    for (int first_label = 0; first_label <= maximum_label; ++first_label) {
        for (int second_label = 0; second_label <= maximum_label;
             ++second_label) {
            for (int third_label = 0; third_label <= maximum_label;
                 ++third_label) {
                for (int first_down = 0; first_down <= first_label;
                     ++first_down) {
                    for (int second_down = 0; second_down <= second_label;
                         ++second_down) {
                        for (int third_down = 0; third_down <= third_label;
                             ++third_down) {
                            const std::vector<CrystalToken> initial{
                                {first_label, first_down},
                                {second_label, second_down},
                                {third_label, third_down}
                            };
                            std::vector<CrystalToken> left = initial;
                            crystal_r_swap(left[0U], left[1U]);
                            crystal_r_swap(left[1U], left[2U]);
                            crystal_r_swap(left[0U], left[1U]);
                            std::vector<CrystalToken> right = initial;
                            crystal_r_swap(right[1U], right[2U]);
                            crystal_r_swap(right[0U], right[1U]);
                            crystal_r_swap(right[1U], right[2U]);
                            ++braid_tests;
                            if (!same_tokens(left, right)) {
                                throw std::runtime_error(
                                    "crystal R violates Yang--Baxter"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "SU2_A1_CRYSTAL_R result=PASS maximum_label="
              << maximum_label << " involutions=" << involutions
              << " braid_tests=" << braid_tests << '\n';
}

bool increment_degrees(Vector& degrees, int maximum_label) {
    for (std::size_t reverse = 0U; reverse < degrees.size(); ++reverse) {
        const std::size_t i = degrees.size() - 1U - reverse;
        if (degrees[i] < maximum_label) {
            ++degrees[i];
            std::fill(
                degrees.begin() + static_cast<std::ptrdiff_t>(i + 1U),
                degrees.end(), 1
            );
            return true;
        }
    }
    return false;
}

struct ColouredTicket {
    Vector degrees;
    Vector top;
    Vector colours;
};

bool same_coloured_ticket(
    const ColouredTicket& left, const ColouredTicket& right
) {
    return left.degrees == right.degrees && left.top == right.top
        && left.colours == right.colours;
}

std::vector<std::size_t> positions_of_colour(
    const Vector& colours, int colour
) {
    std::vector<std::size_t> positions;
    for (std::size_t index = 0U; index < colours.size(); ++index) {
        if (colours[index] == colour) {
            positions.push_back(index);
        }
    }
    return positions;
}

Vector restrict_to_positions(
    const Vector& values, const std::vector<std::size_t>& positions
) {
    Vector result;
    result.reserve(positions.size());
    for (const std::size_t position : positions) {
        result.push_back(values[position]);
    }
    return result;
}

bool valid_coloured_ticket(const ColouredTicket& ticket, int level) {
    if (ticket.degrees.size() != ticket.top.size()
        || ticket.degrees.size() != ticket.colours.size()) {
        return false;
    }
    for (int colour = 0; colour <= 1; ++colour) {
        const std::vector<std::size_t> positions = positions_of_colour(
            ticket.colours, colour
        );
        if (!finite_valid_top(
                restrict_to_positions(ticket.degrees, positions),
                restrict_to_positions(ticket.top, positions), level
            )) {
            return false;
        }
    }
    return true;
}

std::vector<ColouredTicket> coloured_tickets(
    const Vector& degrees, std::uint64_t colour_mask, int level
) {
    ColouredTicket base;
    base.degrees = degrees;
    base.colours.resize(degrees.size(), 0);
    for (std::size_t index = 0U; index < degrees.size(); ++index) {
        base.colours[index] = ((colour_mask >> index) & UINT64_C(1)) != 0U
            ? 1
            : 0;
    }
    const std::vector<std::size_t> zero_positions = positions_of_colour(
        base.colours, 0
    );
    const std::vector<std::size_t> one_positions = positions_of_colour(
        base.colours, 1
    );
    const std::vector<Vector> zero_tops = finite_tops(
        restrict_to_positions(degrees, zero_positions), level
    );
    const std::vector<Vector> one_tops = finite_tops(
        restrict_to_positions(degrees, one_positions), level
    );
    std::vector<ColouredTicket> result;
    for (const Vector& zero_top : zero_tops) {
        for (const Vector& one_top : one_tops) {
            ColouredTicket ticket = base;
            ticket.top.assign(degrees.size(), 0);
            for (std::size_t index = 0U; index < zero_positions.size();
                 ++index) {
                ticket.top[zero_positions[index]] = zero_top[index];
            }
            for (std::size_t index = 0U; index < one_positions.size();
                 ++index) {
                ticket.top[one_positions[index]] = one_top[index];
            }
            if (!valid_coloured_ticket(ticket, level)) {
                throw std::runtime_error("invalid enumerated coloured ticket");
            }
            result.push_back(std::move(ticket));
        }
    }
    return result;
}

ColouredTicket apply_coloured_ticket_swap(
    const ColouredTicket& ticket, std::size_t position, int level
) {
    if (position + 1U >= ticket.degrees.size()) {
        throw std::runtime_error("coloured swap position out of range");
    }
    ColouredTicket result = ticket;
    const int left_colour = ticket.colours[position];
    const int right_colour = ticket.colours[position + 1U];
    std::swap(result.degrees[position], result.degrees[position + 1U]);
    std::swap(result.colours[position], result.colours[position + 1U]);
    if (left_colour != right_colour) {
        std::swap(result.top[position], result.top[position + 1U]);
    } else {
        const std::vector<std::size_t> positions = positions_of_colour(
            ticket.colours, left_colour
        );
        const auto found = std::find(
            positions.begin(), positions.end(), position
        );
        if (found == positions.end() || found + 1 == positions.end()
            || *(found + 1) != position + 1U) {
            throw std::runtime_error("same-colour entries are not adjacent");
        }
        const std::size_t local_position = static_cast<std::size_t>(
            std::distance(positions.begin(), found)
        );
        const Vector local_degrees = restrict_to_positions(
            ticket.degrees, positions
        );
        const Vector local_top = restrict_to_positions(ticket.top, positions);
        const Vector swapped_top = finite_recoupling_swap(
            local_degrees, local_top, local_position, level
        );
        for (std::size_t index = 0U; index < positions.size(); ++index) {
            result.top[positions[index]] = swapped_top[index];
        }
    }
    if (!valid_coloured_ticket(result, level)) {
        throw std::runtime_error("coloured transport left finite path space");
    }
    return result;
}

void coloured_ticket_check(int level, int vertices, int maximum_label) {
    if (vertices >= 63) {
        throw std::runtime_error("coloured-ticket vertices must be below 63");
    }
    std::uint64_t degree_words = 0U;
    std::uint64_t tickets = 0U;
    std::uint64_t swaps = 0U;
    std::uint64_t braid_tests = 0U;
    std::uint64_t distant_tests = 0U;
    Vector degrees(static_cast<std::size_t>(vertices), 1);
    bool more = true;
    while (more) {
        ++degree_words;
        const std::uint64_t colour_count = UINT64_C(1) << vertices;
        for (std::uint64_t colour_mask = 0U; colour_mask < colour_count;
             ++colour_mask) {
            const std::vector<ColouredTicket> current = coloured_tickets(
                degrees, colour_mask, level
            );
            tickets += static_cast<std::uint64_t>(current.size());
            for (const ColouredTicket& ticket : current) {
                for (std::size_t position = 0U;
                     position + 1U < ticket.degrees.size(); ++position) {
                    const ColouredTicket swapped = apply_coloured_ticket_swap(
                        ticket, position, level
                    );
                    const ColouredTicket round_trip =
                        apply_coloured_ticket_swap(swapped, position, level);
                    ++swaps;
                    if (!same_coloured_ticket(ticket, round_trip)) {
                        throw std::runtime_error(
                            "coloured ticket transport is not involutive"
                        );
                    }
                }
                for (std::size_t position = 0U;
                     position + 2U < ticket.degrees.size(); ++position) {
                    const ColouredTicket left = apply_coloured_ticket_swap(
                        apply_coloured_ticket_swap(
                            apply_coloured_ticket_swap(ticket, position, level),
                            position + 1U, level
                        ),
                        position, level
                    );
                    const ColouredTicket right = apply_coloured_ticket_swap(
                        apply_coloured_ticket_swap(
                            apply_coloured_ticket_swap(
                                ticket, position + 1U, level
                            ),
                            position, level
                        ),
                        position + 1U, level
                    );
                    ++braid_tests;
                    if (!same_coloured_ticket(left, right)) {
                        throw std::runtime_error(
                            "coloured ticket transport violates braid"
                        );
                    }
                }
                for (std::size_t left_position = 0U;
                     left_position + 1U < ticket.degrees.size();
                     ++left_position) {
                    for (std::size_t right_position = left_position + 2U;
                         right_position + 1U < ticket.degrees.size();
                         ++right_position) {
                        const ColouredTicket left = apply_coloured_ticket_swap(
                            apply_coloured_ticket_swap(
                                ticket, left_position, level
                            ),
                            right_position, level
                        );
                        const ColouredTicket right = apply_coloured_ticket_swap(
                            apply_coloured_ticket_swap(
                                ticket, right_position, level
                            ),
                            left_position, level
                        );
                        ++distant_tests;
                        if (!same_coloured_ticket(left, right)) {
                            throw std::runtime_error(
                                "coloured ticket transport violates distant commutation"
                            );
                        }
                    }
                }
            }
        }
        more = increment_degrees(degrees, maximum_label);
    }
    std::cout << "SU2_FINITE_COLOURED_TICKETS result=PASS level=" << level
              << " vertices=" << vertices
              << " maximum_label=" << maximum_label
              << " degree_words=" << degree_words
              << " tickets=" << tickets << " swaps=" << swaps
              << " braid_tests=" << braid_tests
              << " distant_tests=" << distant_tests << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "crystal-r") {
            crystal_r_check(parse_nonnegative(argv[2], "maximum label"));
            return 0;
        }
        const bool coloured_ticket_mode = argc == 5
            && std::string(argv[1]) == "coloured-ticket";
        const bool recoupling_mode = argc == 5
            && std::string(argv[1]) == "recoupling";
        const bool coherence_mode = argc == 5
            && std::string(argv[1]) == "coherence";
        if (argc != 4 && !coloured_ticket_mode && !recoupling_mode
            && !coherence_mode) {
            throw std::runtime_error(
                "usage: probe_su2_finite_bk_exchange LEVEL VERTICES MAX_LABEL "
                "| coloured-ticket LEVEL VERTICES MAX_LABEL "
                "| recoupling LEVEL VERTICES MAX_LABEL"
                "| coherence LEVEL VERTICES MAX_LABEL"
                "| crystal-r MAX_LABEL"
            );
        }
        const int offset = (coloured_ticket_mode || recoupling_mode
                            || coherence_mode)
            ? 1
            : 0;
        const int level = parse_nonnegative(argv[1 + offset], "level");
        const int vertices = parse_nonnegative(argv[2 + offset], "vertices");
        const int maximum_label
            = parse_nonnegative(argv[3 + offset], "maximum label");
        if (level <= 0 || vertices <= 1 || maximum_label <= 0
            || maximum_label > level) {
            throw std::runtime_error("bounds must satisfy 1<=MAX_LABEL<=LEVEL and VERTICES>=2");
        }
        if (coloured_ticket_mode) {
            coloured_ticket_check(level, vertices, maximum_label);
            return 0;
        }

        std::uint64_t degree_words = 0U;
        std::uint64_t finite_paths = 0U;
        std::uint64_t tested_swaps = 0U;
        std::uint64_t braid_tests = 0U;
        std::uint64_t distant_tests = 0U;
        Vector degrees(static_cast<std::size_t>(vertices), 1);
        bool more = true;
        while (more) {
            ++degree_words;
            const std::vector<Vector> tops = finite_tops(degrees, level);
            finite_paths += static_cast<std::uint64_t>(tops.size());
            for (const Vector& top : tops) {
                for (std::size_t position = 0U;
                     position + 1U < degrees.size(); ++position) {
                    ++tested_swaps;
                    Vector swapped_degrees = degrees;
                    std::swap(
                        swapped_degrees[position],
                        swapped_degrees[position + 1U]
                    );
                    const Vector swapped_top = (recoupling_mode || coherence_mode)
                        ? finite_recoupling_swap(degrees, top, position, level)
                        : bender_knuth_swap(degrees, top, position);
                    if (!finite_valid_top(swapped_degrees, swapped_top, level)) {
                        std::cout << "SU2_FINITE_BK_ORDINARY_TRANSPORT "
                                  << "result=FAIL level=" << level
                                  << " degrees=";
                        print_vector(degrees);
                        std::cout << " top=";
                        print_vector(top);
                        std::cout << " position=" << position
                                  << " swapped_degrees=";
                        print_vector(swapped_degrees);
                        std::cout << " swapped_top=";
                        print_vector(swapped_top);
                        std::cout << '\n';
                        return 1;
                    }
                    if (recoupling_mode || coherence_mode) {
                        const Vector crystal_top = finite_crystal_swap(
                            degrees, top, position
                        );
                        if (crystal_top != swapped_top) {
                            std::cout << "SU2_FINITE_BK_RECOUPLING "
                                      << "result=CRYSTAL_MISMATCH level="
                                      << level << " degrees=";
                            print_vector(degrees);
                            std::cout << " top=";
                            print_vector(top);
                            std::cout << " position=" << position
                                      << " rank_top=";
                            print_vector(swapped_top);
                            std::cout << " crystal_top=";
                            print_vector(crystal_top);
                            std::cout << '\n';
                            return 1;
                        }
                    }
                    const Vector round_trip = (recoupling_mode || coherence_mode)
                        ? finite_recoupling_swap(
                            swapped_degrees, swapped_top, position, level
                        )
                        : bender_knuth_swap(
                            swapped_degrees, swapped_top, position
                        );
                    if (round_trip != top) {
                        throw std::runtime_error("BK round trip is not involutive");
                    }
                }
                if (coherence_mode) {
                    const FusionPath initial{degrees, top};
                    for (std::size_t position = 0U;
                         position + 2U < degrees.size(); ++position) {
                        const FusionPath left = apply_recoupling_swap(
                            apply_recoupling_swap(
                                apply_recoupling_swap(initial, position, level),
                                position + 1U, level
                            ),
                            position, level
                        );
                        const FusionPath right = apply_recoupling_swap(
                            apply_recoupling_swap(
                                apply_recoupling_swap(initial, position + 1U, level),
                                position, level
                            ),
                            position + 1U, level
                        );
                        ++braid_tests;
                        if (!same_path(left, right)) {
                            std::cout << "SU2_FINITE_BK_RECOUPLING "
                                      << "result=BRAID_FAIL level=" << level
                                      << " degrees=";
                            print_vector(degrees);
                            std::cout << " top=";
                            print_vector(top);
                            std::cout << " position=" << position
                                      << " left_top=";
                            print_vector(left.top);
                            std::cout << " right_top=";
                            print_vector(right.top);
                            std::cout << '\n';
                            return 1;
                        }
                    }
                    for (std::size_t left_position = 0U;
                         left_position + 1U < degrees.size(); ++left_position) {
                        for (std::size_t right_position = left_position + 2U;
                             right_position + 1U < degrees.size();
                             ++right_position) {
                            const FusionPath left = apply_recoupling_swap(
                                apply_recoupling_swap(
                                    initial, left_position, level
                                ),
                                right_position, level
                            );
                            const FusionPath right = apply_recoupling_swap(
                                apply_recoupling_swap(
                                    initial, right_position, level
                                ),
                                left_position, level
                            );
                            ++distant_tests;
                            if (!same_path(left, right)) {
                                std::cout << "SU2_FINITE_BK_RECOUPLING "
                                          << "result=DISTANT_FAIL level=" << level
                                          << " degrees=";
                                print_vector(degrees);
                                std::cout << " top=";
                                print_vector(top);
                                std::cout << " positions=(" << left_position
                                          << ',' << right_position
                                          << ") left_top=";
                                print_vector(left.top);
                                std::cout << " right_top=";
                                print_vector(right.top);
                                std::cout << '\n';
                                return 1;
                            }
                        }
                    }
                }
            }
            more = increment_degrees(degrees, maximum_label);
        }
        std::cout << ((recoupling_mode || coherence_mode)
                          ? "SU2_FINITE_BK_RECOUPLING"
                          : "SU2_FINITE_BK_ORDINARY_TRANSPORT")
                  << " result=PASS level="
                  << level << " vertices=" << vertices
                  << " maximum_label=" << maximum_label
                  << " degree_words=" << degree_words
                  << " finite_paths=" << finite_paths
                  << " tested_swaps=" << tested_swaps;
        if (coherence_mode) {
            std::cout << " braid_tests=" << braid_tests
                      << " distant_tests=" << distant_tests;
        }
        std::cout << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
