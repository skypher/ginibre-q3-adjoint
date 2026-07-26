#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

namespace {

template <class F>
void outputs(const int k, const int a, const int b, F f) {
    const int upper = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= upper; c += 2) {
        f(c);
    }
}

std::vector<int> fusion_interval(const int k, const int a, const int b) {
    std::vector<int> result;
    outputs(k, a, b, [&](const int c) {
        result.push_back(c);
    });
    return result;
}

bool contains_output(const int k, const int a, const int b, const int target) {
    return std::abs(a - b) <= target
        && target <= std::min(a + b, 2 * k - a - b)
        && ((a + b + target) % 2 == 0);
}

int profile(
    const int k,
    const std::vector<int>& lhs,
    const std::vector<int>& rhs,
    const int target
) {
    int count = 0;
    for (const int x : lhs) {
        for (const int y : rhs) {
            count += contains_output(k, x, y, target) ? 1 : 0;
        }
    }
    return count;
}

int overlap_size(const std::vector<int>& lhs, const std::vector<int>& rhs) {
    int count = 0;
    for (const int x : lhs) {
        count += std::binary_search(rhs.begin(), rhs.end(), x) ? 1 : 0;
    }
    return count;
}

int cap(const int k, const int label) {
    return std::min(label, k - label) + 1;
}

}  // namespace

int main() {
    constexpr int min_k = 6;
    constexpr int max_k = 32;
    std::array<std::int64_t, 3> checked{};
    std::int64_t two_by_three_checked = 0;
    std::int64_t rank_two_six_checked = 0;
    std::int64_t rank_two_nonzero_six_checked = 0;
    std::int64_t rank_two_one_nonzero_six_checked = 0;
    std::int64_t rank_two_distinct_six_checked = 0;
    std::int64_t rank_four_short_checked = 0;
    std::int64_t one_by_three_checked = 0;
    std::int64_t negative_controls = 0;
    std::array<int, 5> distinct_payment_minimum{};
    std::array<std::array<int, 5>, 5> distinct_payment_witness{};
    std::array<std::array<int, 4>, 5> repeated_payment_minimum{};
    distinct_payment_minimum.fill(std::numeric_limits<int>::max());
    for (auto& row : repeated_payment_minimum) {
        row.fill(std::numeric_limits<int>::max());
    }

    for (int k = min_k; k <= max_k; ++k) {
        for (int a = 0; a <= k; ++a) {
            for (int b = a; b <= k; ++b) {
                const auto lhs = fusion_interval(k, a, b);
                if (lhs.size() < 2U) {
                    continue;
                }
                for (int c = 0; c <= k; ++c) {
                    for (int d = c; d <= k; ++d) {
                        const auto rhs = fusion_interval(k, c, d);
                        if (rhs.size() < 2U) {
                            continue;
                        }
                        const int overlap = overlap_size(lhs, rhs);
                        if (overlap < 1 || overlap > 4) {
                            continue;
                        }
                        const int b0 = profile(k, lhs, rhs, 0);
                        const int b2 = profile(k, lhs, rhs, 2);
                        const int b4 = profile(k, lhs, rhs, 4);
                        const int b6 = profile(k, lhs, rhs, 6);
                        if (lhs.size() >= 3U && rhs.size() >= 3U) {
                            if (overlap >= 3 && overlap <= 4
                                && k >= 7
                                && a != b
                                && ((a != c && a != d)
                                    || (b != c && b != d))
                                && !(std::max(lhs.front(), rhs.front()) == 1
                                    && std::min(lhs.back(), rhs.back()) == k)) {
                                const std::size_t rank =
                                    static_cast<std::size_t>(overlap);
                                const int distinct_payment =
                                    b0 + 2 * b2 + 3 * b4 + 2 * b6;
                                if (distinct_payment
                                    < distinct_payment_minimum[rank]) {
                                    distinct_payment_minimum[rank] =
                                        distinct_payment;
                                    distinct_payment_witness[rank] =
                                        {k, a, b, c, d};
                                }
                                const std::array<int, 4> repeated{{
                                    b0 + 4 * b2 + 2 * b4 + b6,
                                    b0 + 2 * b2 + 4 * b4 + b6,
                                    b0 + 2 * b2 + 3 * b4 + 4 * b6,
                                    b0 + 2 * b2 + 3 * b4 + 3 * b6
                                }};
                                if (overlap == 4 && repeated[0] < 80) {
                                    const bool first_geometry =
                                        k == 8
                                        && lhs.front() == 1
                                        && lhs.back() == 7
                                        && rhs.front() == 1
                                        && rhs.back() == 7;
                                    const bool second_geometry =
                                        k == 8
                                        && ((lhs.front() == 2
                                             && lhs.back() == 8
                                             && rhs.front() == 0
                                             && rhs.back() == 8)
                                            || (rhs.front() == 2
                                                && rhs.back() == 8
                                                && lhs.front() == 0
                                                && lhs.back() == 8));
                                    if (!first_geometry && !second_geometry) {
                                        std::cerr
                                            << "rank-four short-profile"
                                            << " classification failure k="
                                            << k << " pairs=[" << a << ','
                                            << b << "]|[" << c << ',' << d
                                            << "] profile=(" << b0 << ','
                                            << b2 << ',' << b4 << ',' << b6
                                            << ")\n";
                                        return EXIT_FAILURE;
                                    }
                                    ++rank_four_short_checked;
                                }
                                for (std::size_t profile_index = 0U;
                                     profile_index < repeated.size();
                                     ++profile_index) {
                                    repeated_payment_minimum[rank]
                                        [profile_index] = std::min(
                                            repeated_payment_minimum[rank]
                                                [profile_index],
                                            repeated[profile_index]
                                        );
                                }
                            }
                            if (overlap >= 3) {
                                continue;
                            }
                            const bool pass = b0 == overlap
                                && (overlap == 1
                                    ? b2 >= 3 && b4 >= 6
                                    : b2 >= 6 && b4 >= 6);
                            if (!pass) {
                                std::cerr
                                    << "band failure k=" << k
                                    << " pairs=[" << a << ',' << b << "]|["
                                    << c << ',' << d << "] d=" << overlap
                                    << " profile=(" << b0 << ',' << b2 << ','
                                    << b4 << ")\n";
                                return EXIT_FAILURE;
                            }
                            if (overlap == 2) {
                                if (b6 < 3) {
                                    std::cerr
                                        << "rank-two output-six failure k="
                                        << k << " pairs=[" << a << ',' << b
                                        << "]|[" << c << ',' << d
                                        << "] profile=(" << b0 << ',' << b2
                                        << ',' << b4 << ',' << b6 << ")\n";
                                    return EXIT_FAILURE;
                                }
                                ++rank_two_six_checked;
                                if (lhs.front() != 0
                                    && rhs.front() != 0) {
                                    if (b6 < 5) {
                                        std::cerr
                                            << "nonzero rank-two output-six"
                                            << " failure k=" << k
                                            << " pairs=[" << a << ',' << b
                                            << "]|[" << c << ',' << d
                                            << "] profile=(" << b0 << ','
                                            << b2 << ',' << b4 << ',' << b6
                                            << ")\n";
                                        return EXIT_FAILURE;
                                    }
                                    ++rank_two_nonzero_six_checked;
                                }
                                if (k >= 7
                                    && (lhs.front() != 0
                                        || rhs.front() != 0)) {
                                    if (b6 < 5) {
                                        std::cerr
                                            << "one-nonzero rank-two"
                                            << " output-six failure k=" << k
                                            << " pairs=[" << a << ',' << b
                                            << "]|[" << c << ',' << d
                                            << "] profile=(" << b0 << ','
                                            << b2 << ',' << b4 << ',' << b6
                                            << ")\n";
                                        return EXIT_FAILURE;
                                    }
                                    ++rank_two_one_nonzero_six_checked;
                                }
                                if (k >= 10
                                    && lhs.front() != 0
                                    && rhs.front() != 0) {
                                    if (b6 < 6) {
                                        std::cerr
                                            << "distinct rank-two output-six"
                                            << " failure k=" << k
                                            << " pairs=[" << a << ',' << b
                                            << "]|[" << c << ',' << d
                                            << "] profile=(" << b0 << ','
                                            << b2 << ',' << b4 << ',' << b6
                                            << ")\n";
                                        return EXIT_FAILURE;
                                    }
                                    ++rank_two_distinct_six_checked;
                                }
                            }
                            ++checked[static_cast<std::size_t>(overlap)];
                        } else {
                            if (overlap == 1
                                && ((lhs.size() == 2U && rhs.size() >= 3U
                                     && std::max(cap(k, a), cap(k, b)) >= 3)
                                    || (rhs.size() == 2U
                                        && lhs.size() >= 3U
                                        && std::max(cap(k, c), cap(k, d)) >= 3))
                                && !(std::binary_search(lhs.begin(), lhs.end(), 0)
                                    && std::binary_search(
                                        rhs.begin(), rhs.end(), 0
                                    ))) {
                                if (b2 < 3 || b4 < 3) {
                                    std::cerr
                                        << "1x3 band failure k=" << k
                                        << " pairs=[" << a << ',' << b
                                        << "]|[" << c << ',' << d
                                        << "] profile=(" << b0 << ',' << b2
                                        << ',' << b4 << ")\n";
                                    return EXIT_FAILURE;
                                }
                                ++one_by_three_checked;
                            }
                            if (overlap == 2
                                && ((lhs.size() == 2U && rhs.size() >= 3U
                                     && std::max(cap(k, a), cap(k, b)) >= 3)
                                    || (rhs.size() == 2U
                                        && lhs.size() >= 3U
                                        && std::max(cap(k, c), cap(k, d)) >= 3))
                                && !(std::binary_search(lhs.begin(), lhs.end(), 0)
                                    && std::binary_search(
                                        rhs.begin(), rhs.end(), 0
                                    ))) {
                                if (b2 < 5 || b4 < 5) {
                                    std::cerr
                                        << "2x3 band failure k=" << k
                                        << " pairs=[" << a << ',' << b
                                        << "]|[" << c << ',' << d
                                        << "] profile=(" << b0 << ',' << b2
                                        << ',' << b4 << ")\n";
                                    return EXIT_FAILURE;
                                }
                                ++two_by_three_checked;
                            }
                            const bool would_pass = overlap == 1
                                ? b2 >= 3 && b4 >= 6
                                : b2 >= 6 && b4 >= 6;
                            negative_controls += would_pass ? 0 : 1;
                        }
                    }
                }
            }
        }
    }

    if (checked[1] == 0 || checked[2] == 0
        || rank_two_six_checked == 0
        || rank_two_nonzero_six_checked == 0
        || rank_two_one_nonzero_six_checked == 0
        || rank_two_distinct_six_checked == 0
        || rank_four_short_checked == 0
        || one_by_three_checked == 0 || two_by_three_checked == 0
        || negative_controls == 0) {
        std::cerr << "coverage failure d1=" << checked[1]
                  << " d2=" << checked[2]
                  << " d2_b6=" << rank_two_six_checked
                  << " d2_nonzero_b6="
                  << rank_two_nonzero_six_checked
                  << " d2_one_nonzero_b6="
                  << rank_two_one_nonzero_six_checked
                  << " d2_distinct_b6="
                  << rank_two_distinct_six_checked
                  << " d4_short=" << rank_four_short_checked
                  << " one_by_three=" << one_by_three_checked
                  << " two_by_three=" << two_by_three_checked
                  << " negative_controls=" << negative_controls << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "SU2_INTERVAL_BAND3 levels=" << (max_k - min_k + 1)
              << " d1=" << checked[1]
              << " d2=" << checked[2]
              << " d2_b6=" << rank_two_six_checked
              << " d2_nonzero_b6="
              << rank_two_nonzero_six_checked
              << " d2_one_nonzero_b6="
              << rank_two_one_nonzero_six_checked
              << " d2_distinct_b6="
              << rank_two_distinct_six_checked
              << " d4_short=" << rank_four_short_checked
              << " one_by_three=" << one_by_three_checked
              << " two_by_three=" << two_by_three_checked
              << " shallow_rejected=" << negative_controls
              << " d3_distinct_payment="
              << distinct_payment_minimum[3]
              << " d3_witness=("
              << distinct_payment_witness[3][0] << ';'
              << distinct_payment_witness[3][1] << ','
              << distinct_payment_witness[3][2] << '|'
              << distinct_payment_witness[3][3] << ','
              << distinct_payment_witness[3][4] << ')'
              << " d3_repeated_payments=("
              << repeated_payment_minimum[3][0] << ','
              << repeated_payment_minimum[3][1] << ','
              << repeated_payment_minimum[3][2] << ','
              << repeated_payment_minimum[3][3] << ')'
              << " d4_distinct_payment="
              << distinct_payment_minimum[4]
              << " d4_witness=("
              << distinct_payment_witness[4][0] << ';'
              << distinct_payment_witness[4][1] << ','
              << distinct_payment_witness[4][2] << '|'
              << distinct_payment_witness[4][3] << ','
              << distinct_payment_witness[4][4] << ')'
              << " d4_repeated_payments=("
              << repeated_payment_minimum[4][0] << ','
              << repeated_payment_minimum[4][1] << ','
              << repeated_payment_minimum[4][2] << ','
              << repeated_payment_minimum[4][3] << ')'
              << " result=PASS\n";
    return EXIT_SUCCESS;
}
