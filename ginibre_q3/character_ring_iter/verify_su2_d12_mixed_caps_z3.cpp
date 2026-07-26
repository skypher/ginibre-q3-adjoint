#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <z3++.h>

namespace {

z3::expr zabs(const z3::expr& value) {
    return z3::ite(value >= 0, value, -value);
}

z3::expr zmin(const z3::expr& lhs, const z3::expr& rhs) {
    return z3::ite(lhs <= rhs, lhs, rhs);
}

z3::expr zmax(const z3::expr& lhs, const z3::expr& rhs) {
    return z3::ite(lhs >= rhs, lhs, rhs);
}

z3::expr fusion(
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& target
) {
    return zabs(a - b) <= target
        && target <= a + b
        && target <= 2 * k - a - b
        && (a + b + target) % 2 == 0;
}

z3::expr interval_overlap_at_least(
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& c,
    const z3::expr& d,
    const int depth
) {
    const z3::expr lo1 = zabs(a - b);
    const z3::expr hi1 = zmin(a + b, 2 * k - a - b);
    const z3::expr lo2 = zabs(c - d);
    const z3::expr hi2 = zmin(c + d, 2 * k - c - d);
    const z3::expr lo = zmax(lo1, lo2);
    const z3::expr hi = zmin(hi1, hi2);
    return lo1 % 2 == lo2 % 2 && lo + 2 * (depth - 1) <= hi;
}

z3::expr labels_equal(
    z3::context& ctx,
    const std::array<z3::expr, 5>& labels,
    const std::array<z3::expr, 5>& values
) {
    z3::expr result = ctx.bool_val(true);
    for (std::size_t i = 0; i < labels.size(); ++i) {
        result = result && labels[i] == values[i];
    }
    return result;
}

bool fusion_int(
    const int k,
    const int a,
    const int b,
    const int target
) {
    return std::abs(a - b) <= target
        && target <= std::min(a + b, 2 * k - a - b)
        && (a + b + target) % 2 == 0;
}

int fourfold_rank_int(
    const int k,
    const int a,
    const int b,
    const int c,
    const int d
) {
    const int lo1 = std::abs(a - b);
    const int hi1 = std::min(a + b, 2 * k - a - b);
    const int lo2 = std::abs(c - d);
    const int hi2 = std::min(c + d, 2 * k - c - d);
    const int lo = std::max(lo1, lo2);
    const int hi = std::min(hi1, hi2);
    if (lo1 % 2 != lo2 % 2 || lo > hi) {
        return 0;
    }
    return (hi - lo) / 2 + 1;
}

bool check_base_witness() {
    constexpr int k = 13;
    constexpr std::array<int, 2> minus{4, 6};
    constexpr std::array<int, 5> plus{10, 10, 10, 10, 10};
    std::array<int, 2> colors{};
    int demand = 0;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            std::array<int, 3> rest{};
            int next = 0;
            for (int h = 0; h < 5; ++h) {
                if (h != i && h != j) {
                    rest[static_cast<std::size_t>(next)] = h;
                    ++next;
                }
            }
            for (std::size_t orientation = 0;
                 orientation < minus.size(); ++orientation) {
                const int r = minus[orientation];
                const int s = minus[1U - orientation];
                if (!fusion_int(k, plus[static_cast<std::size_t>(i)],
                                plus[static_cast<std::size_t>(j)], r)) {
                    continue;
                }
                const int rank = fourfold_rank_int(
                    k, s, plus[static_cast<std::size_t>(rest[0])],
                    plus[static_cast<std::size_t>(rest[1])],
                    plus[static_cast<std::size_t>(rest[2])]
                );
                if (rank <= 0 || rank > 2) {
                    return false;
                }
                ++colors[orientation];
                demand += rank;
            }
        }
    }
    return colors[0] == 10 && colors[1] == 10 && demand == 30;
}

bool check_n2_incidence_witness() {
    constexpr int k = 9;
    constexpr std::array<int, 2> minus{3, 5};
    constexpr std::array<int, 5> plus{1, 1, 4, 4, 4};
    std::array<int, 2> colors{};
    int demand = 0;
    bool rank_two_incidence_one = false;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            std::array<int, 3> rest{};
            int next = 0;
            for (int h = 0; h < 5; ++h) {
                if (h != i && h != j) {
                    rest[static_cast<std::size_t>(next)] = h;
                    ++next;
                }
            }
            for (std::size_t orientation = 0;
                 orientation < minus.size(); ++orientation) {
                if (!fusion_int(
                        k, plus[static_cast<std::size_t>(i)],
                        plus[static_cast<std::size_t>(j)],
                        minus[orientation])) {
                    continue;
                }
                const int rank = fourfold_rank_int(
                    k, minus[1U - orientation],
                    plus[static_cast<std::size_t>(rest[0])],
                    plus[static_cast<std::size_t>(rest[1])],
                    plus[static_cast<std::size_t>(rest[2])]
                );
                if (rank <= 0) {
                    continue;
                }
                if (rank > 2) {
                    return false;
                }
                ++colors[orientation];
                demand += rank;
                const int incidence =
                    (plus[static_cast<std::size_t>(i)] == 1
                     || plus[static_cast<std::size_t>(i)] == k - 1)
                    + (plus[static_cast<std::size_t>(j)] == 1
                       || plus[static_cast<std::size_t>(j)] == k - 1);
                rank_two_incidence_one =
                    rank_two_incidence_one
                    || (rank == 2 && incidence == 1);
            }
        }
    }
    return colors[0] == 6 && colors[1] == 6
        && demand == 24 && rank_two_incidence_one;
}

struct MaskResult {
    bool passed = false;
    std::string diagnostic;
};

MaskResult verify_mask(const int parity_mask) {
    z3::context ctx;
    z3::solver solver(ctx);

    const auto label = [&](const char* name, const int shift) {
        const std::string half_name = std::string(name) + "_half";
        return 2 * ctx.int_const(half_name.c_str())
            + ((parity_mask >> shift) & 1);
    };
    const z3::expr k = label("k", 0);
    const z3::expr q = label("q", 1);
    const z3::expr a = label("a", 2);
    const std::array<z3::expr, 5> plus{
        label("p0", 3),
        label("p1", 4),
        label("p2", 5),
        label("p3", 6),
        label("p4", 7)
    };

    solver.add(k >= 9);
    solver.add(q >= 2 && q <= k - 2);
    solver.add(a >= 2 && a <= k - 2);
    solver.add(q < a);
    for (std::size_t i = 0; i < plus.size(); ++i) {
        solver.add(plus[i] >= 1 && plus[i] <= k - 1);
        solver.add(plus[i] != q && plus[i] != a);
        if (i + 1U < plus.size()) {
            solver.add(plus[i] <= plus[i + 1U]);
        }
    }

    std::vector<z3::expr> active;
    std::vector<z3::expr> rank;
    std::vector<int> orientation;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            std::array<int, 3> rest{};
            int next = 0;
            for (int h = 0; h < 5; ++h) {
                if (h != i && h != j) {
                    rest[static_cast<std::size_t>(next)] = h;
                    ++next;
                }
            }
            for (int o = 0; o < 2; ++o) {
                const z3::expr r = o == 0 ? q : a;
                const z3::expr s = o == 0 ? a : q;
                const z3::expr complement_nonzero =
                    interval_overlap_at_least(
                    k, s,
                    plus[static_cast<std::size_t>(rest[0])],
                    plus[static_cast<std::size_t>(rest[1])],
                    plus[static_cast<std::size_t>(rest[2])], 1
                );
                const z3::expr complement_rank_two =
                    interval_overlap_at_least(
                    k, s,
                    plus[static_cast<std::size_t>(rest[0])],
                    plus[static_cast<std::size_t>(rest[1])],
                    plus[static_cast<std::size_t>(rest[2])], 2
                );
                const z3::expr complement_rank_three =
                    interval_overlap_at_least(
                    k, s,
                    plus[static_cast<std::size_t>(rest[0])],
                    plus[static_cast<std::size_t>(rest[1])],
                    plus[static_cast<std::size_t>(rest[2])], 3
                );
                const z3::expr is_active = fusion(
                    k, plus[static_cast<std::size_t>(i)],
                    plus[static_cast<std::size_t>(j)], r
                ) && complement_nonzero;
                const z3::expr complement_rank = z3::ite(
                    complement_rank_two, ctx.int_val(2),
                    z3::ite(
                        complement_nonzero, ctx.int_val(1),
                        ctx.int_val(0)
                    )
                );
                solver.add(z3::implies(
                    is_active, !complement_rank_three
                ));
                active.push_back(is_active);
                rank.push_back(complement_rank);
                orientation.push_back(o);
            }
        }
    }

    const z3::expr zero = ctx.int_val(0);
    const z3::expr one = ctx.int_val(1);
    z3::expr c = zero;
    z3::expr h = zero;
    z3::expr demand = zero;
    std::array<z3::expr, 2> color_count{zero, zero};
    for (std::size_t edge = 0; edge < active.size(); ++edge) {
        c = c + z3::ite(active[edge], one, zero);
        h = h + z3::ite(
            active[edge] && rank[edge] == 2, one, zero
        );
        demand = demand + z3::ite(
            active[edge], rank[edge], zero
        );
        const std::size_t color =
            static_cast<std::size_t>(orientation[edge]);
        color_count[color] = color_count[color] + z3::ite(
            active[edge], one, zero
        );
    }
    solver.add(color_count[0] > 0 && color_count[1] > 0);

    z3::expr cap_two_count = zero;
    for (const z3::expr& label_value : plus) {
        cap_two_count = cap_two_count + z3::ite(
            label_value == 1 || label_value == k - 1, one, zero
        );
    }
    z3::expr rank_two_incidence_one = ctx.bool_val(false);
    std::size_t edge = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            const z3::expr cap_i =
                plus[i] == 1 || plus[i] == k - 1;
            const z3::expr cap_j =
                plus[j] == 1 || plus[j] == k - 1;
            const z3::expr incidence_one =
                (cap_i && !cap_j) || (!cap_i && cap_j);
            for (int orientation_index = 0;
                 orientation_index < 2; ++orientation_index) {
                rank_two_incidence_one =
                    rank_two_incidence_one
                    || (incidence_one
                        && active[edge]
                        && rank[edge] == 2);
                ++edge;
            }
        }
    }
    const z3::expr dense =
        (h == 0 && c >= 8) || (h > 0 && demand >= 12);
    const z3::expr one_label = ctx.int_val(1);
    const z3::expr two_label = ctx.int_val(2);
    const z3::expr four_label = ctx.int_val(4);
    const z3::expr j_one = k - 1;
    const std::array<z3::expr, 5> cap_three_low{
        one_label, one_label, one_label, four_label, four_label
    };
    const std::array<z3::expr, 5> cap_three_high{
        four_label, four_label, j_one, j_one, j_one
    };
    const std::array<z3::expr, 5> cap_three_bridge_low{
        one_label, one_label, two_label, two_label, j_one
    };
    const std::array<z3::expr, 5> cap_three_bridge_high{
        one_label, two_label, two_label, j_one, j_one
    };
    const std::array<z3::expr, 5> cap_three_middle_low{
        one_label, one_label, q + 1, q + 1, j_one
    };
    const std::array<z3::expr, 5> cap_three_middle_high{
        one_label, one_label, a - 1, a - 1, j_one
    };
    const z3::expr family_three_0 =
        q == 2 && a == 3
        && labels_equal(ctx, plus, cap_three_low);
    const z3::expr family_three_1 =
        q == 2 && a == 5
        && labels_equal(ctx, plus, cap_three_low);
    const z3::expr family_three_2 =
        q == 2 && a == k - 5
        && labels_equal(ctx, plus, cap_three_high);
    const z3::expr family_three_3 =
        q == 2 && a == k - 3
        && labels_equal(ctx, plus, cap_three_high);
    const z3::expr family_three_4 =
        q == 3 && a == k - 2
        && labels_equal(ctx, plus, cap_three_bridge_low);
    const z3::expr family_three_5 =
        q == k - 3 && a == k - 2
        && labels_equal(ctx, plus, cap_three_bridge_high);
    const z3::expr family_three_6 =
        2 * q == k - 3 && 2 * a == k + 1
        && labels_equal(ctx, plus, cap_three_middle_low);
    const z3::expr family_three_7 =
        2 * q == k - 1 && 2 * a == k + 3
        && labels_equal(ctx, plus, cap_three_middle_high);
    const z3::expr classified_three =
        family_three_0 || family_three_1 || family_three_2
        || family_three_3 || family_three_4 || family_three_5
        || family_three_6 || family_three_7;
    const z3::expr three_demand_table =
        (family_three_0 && demand == 13)
        || (family_three_1 && demand == 12)
        || (family_three_2 && demand == 12)
        || (family_three_3 && demand == 13)
        || (family_three_4 && demand == 12)
        || (family_three_5 && demand == 12)
        || (family_three_6 && demand == 12)
        || (family_three_7
            && ((k == 9 && demand == 13)
                || (k >= 11 && demand == 12)));
    solver.add(
        (cap_two_count == 0 && demand > 32)
        || (cap_two_count == 1 && demand > 24)
        || (cap_two_count == 2 && demand > 27)
        || (cap_two_count == 3 && demand > 13)
        || (cap_two_count >= 4 && dense)
        || (cap_two_count == 0 && h == 0 && c > 10)
        || (cap_two_count == 1 && h == 0 && c > 10)
        || (cap_two_count == 2 && h == 0 && c > 9)
        || (cap_two_count == 3 && h == 0 && c >= 8)
        || (cap_two_count == 2 && h > 0 && dense
            && !rank_two_incidence_one)
        || (cap_two_count == 3 && dense && !classified_three)
        || (cap_two_count == 3 && dense && !three_demand_table)
    );
    const z3::check_result result = solver.check();
    if (result != z3::unsat) {
        std::ostringstream message;
        message << "mask=" << parity_mask
                << " result=" << result;
        if (result == z3::sat) {
            message << " model=" << solver.get_model();
        }
        return {false, message.str()};
    }
    return {true, {}};
}

unsigned worker_limit() {
    const unsigned hardware = std::max(
        1U, std::thread::hardware_concurrency()
    );
    double load_average = 0.0;
    const int load_status = getloadavg(&load_average, 1);
    const unsigned loaded = load_status == 1
        ? static_cast<unsigned>(std::max(0.0, std::ceil(load_average)))
        : 0U;
    const unsigned cpu_limit =
        hardware > loaded ? hardware - loaded : 1U;

    std::ifstream memory("/proc/meminfo");
    std::string key;
    std::uint64_t value = 0;
    std::string unit;
    std::uint64_t available_kib = 0;
    while (memory >> key >> value >> unit) {
        if (key == "MemAvailable:") {
            available_kib = value;
            break;
        }
    }
    constexpr std::uint64_t kib_per_worker = 128U * 1024U;
    const unsigned memory_limit = available_kib == 0
        ? hardware
        : static_cast<unsigned>(std::max<std::uint64_t>(
            1U, available_kib / kib_per_worker
        ));
    return std::max(1U, std::min(cpu_limit, memory_limit));
}

}  // namespace

int main(int argc, char** argv) {
    if (!check_base_witness() || !check_n2_incidence_witness()) {
        std::cerr << "SU2_D12_MIXED_CAPS_Z3 FAIL base_witness\n";
        return EXIT_FAILURE;
    }
    if (argc == 2) {
        const int parity_mask = std::atoi(argv[1]);
        if (parity_mask < 0 || parity_mask >= 256) {
            std::cerr << "parity mask must lie in [0,255]\n";
            return EXIT_FAILURE;
        }
        const MaskResult result = verify_mask(parity_mask);
        if (!result.passed) {
            std::cerr << "SU2_D12_MIXED_CAPS_Z3 FAIL "
                      << result.diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout
            << "SU2_D12_MIXED_CAPS_Z3 mask=" << parity_mask
            << " k>=9 witness_T=30"
            << " placement_counterexample=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    }
    if (argc != 1) {
        std::cerr << "usage: verify_su2_d12_mixed_caps_z3 [PARITY_MASK]\n";
        return EXIT_FAILURE;
    }

    const unsigned workers = worker_limit();
    std::atomic<int> next_mask{0};
    std::atomic<int> completed{0};
    std::atomic<bool> failed{false};
    std::mutex diagnostic_mutex;
    std::string diagnostic;
    std::vector<std::jthread> pool;
    pool.reserve(workers);
    for (unsigned worker = 0; worker < workers; ++worker) {
        pool.emplace_back([&]() {
            while (!failed.load()) {
                const int mask = next_mask.fetch_add(1);
                if (mask >= 256) {
                    return;
                }
                const MaskResult result = verify_mask(mask);
                if (!result.passed) {
                    {
                        std::lock_guard<std::mutex> lock(diagnostic_mutex);
                        if (diagnostic.empty()) {
                            diagnostic = result.diagnostic;
                        }
                    }
                    failed.store(true);
                    return;
                }
                completed.fetch_add(1);
            }
        });
    }
    pool.clear();

    if (failed.load() || completed.load() != 256) {
        std::cerr << "SU2_D12_MIXED_CAPS_Z3 FAIL"
                  << " completed=" << completed.load()
                  << " diagnostic=" << diagnostic << '\n';
        return EXIT_FAILURE;
    }
    std::cout
        << "SU2_D12_MIXED_CAPS_Z3 k>=9 masks=256"
        << " workers=" << workers
        << " witness_T=30"
        << " bounds=n0:32,n1:24,n2:27,n3:13,n4+:no-dense"
        << " rank1=n0:10,n1:10,n2:9,n3:no-dense"
        << " n2_rank2=has-incidence-one"
        << " n3_families=8"
        << " result=PASS\n";
    return EXIT_SUCCESS;
}
