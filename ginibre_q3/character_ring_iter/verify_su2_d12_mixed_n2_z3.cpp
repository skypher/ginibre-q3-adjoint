#include <array>
#include <algorithm>
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

template <class F>
void outputs(const int k, const int a, const int b, F f) {
    const int upper = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= upper; c += 2) {
        f(c);
    }
}

bool contains_output(const int k, const int a, const int b, const int target) {
    return std::abs(a - b) <= target
        && target <= std::min(a + b, 2 * k - a - b)
        && ((a + b + target) % 2 == 0);
}

std::vector<std::int64_t> decomposition(
    const int k,
    const std::vector<int>& labels
) {
    std::vector<std::int64_t> cur(static_cast<std::size_t>(k + 1));
    std::vector<std::int64_t> next(static_cast<std::size_t>(k + 1));
    cur[0] = 1;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int x = 0; x <= k; ++x) {
            const auto value = cur[static_cast<std::size_t>(x)];
            if (value == 0) {
                continue;
            }
            outputs(k, x, label, [&](const int y) {
                next[static_cast<std::size_t>(y)] += value;
            });
        }
        cur.swap(next);
    }
    return cur;
}

std::int64_t check_small_levels() {
    std::int64_t checked = 0;
    for (int k = 6; k <= 8; ++k) {
        for (int q = 2; q <= k - 2; ++q) {
            for (int a = q + 1; a <= k - 2; ++a) {
                for (int p0 = 1; p0 <= k - 1; ++p0)
                for (int p1 = p0; p1 <= k - 1; ++p1)
                for (int p2 = p1; p2 <= k - 1; ++p2)
                for (int p3 = p2; p3 <= k - 1; ++p3)
                for (int p4 = p3; p4 <= k - 1; ++p4) {
                    const std::array<int, 5> plus{p0, p1, p2, p3, p4};
                    bool disjoint = true;
                    for (const int p : plus) {
                        disjoint = disjoint && p != q && p != a;
                    }
                    if (!disjoint) {
                        continue;
                    }
                    int c = 0;
                    int h = 0;
                    std::array<int, 2> colors{};
                    std::int64_t best_local = 0;
                    bool valid = true;
                    for (std::size_t i = 0; i < plus.size(); ++i) {
                        for (std::size_t j = i + 1; j < plus.size(); ++j) {
                            std::vector<int> rest;
                            for (std::size_t u = 0; u < plus.size(); ++u) {
                                if (u != i && u != j) {
                                    rest.push_back(plus[u]);
                                }
                            }
                            for (std::size_t o = 0; o < 2U; ++o) {
                                const int r = o == 0U ? q : a;
                                const int s = o == 0U ? a : q;
                                if (!contains_output(k, plus[i], plus[j], r)) {
                                    continue;
                                }
                                auto block = rest;
                                block.push_back(s);
                                const auto complement =
                                    decomposition(k, block);
                                const auto rank = complement[0];
                                if (rank == 0) {
                                    continue;
                                }
                                if (rank > 2) {
                                    valid = false;
                                    continue;
                                }
                                ++c;
                                ++colors[o];
                                h += rank == 2 ? 1 : 0;
                                const auto active_profile = decomposition(
                                    k, {r, plus[i], plus[j]}
                                );
                                std::int64_t local = 0;
                                for (int output = 0;
                                     output <= std::min(k, 4); ++output) {
                                    local += active_profile[
                                        static_cast<std::size_t>(output)
                                    ] * complement[
                                        static_cast<std::size_t>(output)
                                    ];
                                }
                                best_local = std::max(best_local, local);
                            }
                        }
                    }
                    if (!valid || colors[0] == 0 || colors[1] == 0) {
                        continue;
                    }
                    const int demand = c + h;
                    const bool dense =
                        (h == 0 && c >= 8) || (h > 0 && demand >= 12);
                    if (!dense) {
                        continue;
                    }
                    ++checked;
                    if (best_local < demand) {
                        std::cerr << "small-level failure k=" << k
                                  << " minus=[" << q << ',' << a
                                  << "] plus=[" << p0 << ',' << p1 << ','
                                  << p2 << ',' << p3 << ',' << p4
                                  << "] local=" << best_local
                                  << " demand=" << demand << '\n';
                        return -1;
                    }
                }
            }
        }
    }
    return checked;
}

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

z3::expr interval_rank(
    z3::context& ctx,
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& c,
    const z3::expr& d
) {
    const z3::expr lo1 = zabs(a - b);
    const z3::expr hi1 = zmin(a + b, 2 * k - a - b);
    const z3::expr lo2 = zabs(c - d);
    const z3::expr hi2 = zmin(c + d, 2 * k - c - d);
    const z3::expr lo = zmax(lo1, lo2);
    const z3::expr hi = zmin(hi1, hi2);
    return z3::ite(
        lo1 % 2 == lo2 % 2 && lo <= hi,
        (hi - lo) / 2 + 1,
        ctx.int_val(0)
    );
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

z3::expr step_two_count(
    z3::context& ctx,
    const z3::expr& base,
    const z3::expr& raw_lo,
    const z3::expr& raw_hi,
    const z3::expr& compatible
) {
    const z3::expr lo = z3::ite(
        (raw_lo - base) % 2 == 0, raw_lo, raw_lo + 1
    );
    const z3::expr hi = z3::ite(
        (raw_hi - base) % 2 == 0, raw_hi, raw_hi - 1
    );
    return z3::ite(
        compatible && lo <= hi,
        (hi - lo) / 2 + 1,
        ctx.int_val(0)
    );
}

z3::expr four_tensor_output(
    z3::context& ctx,
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& c,
    const z3::expr& d,
    const int target
) {
    const z3::expr lo1 = zabs(a - b);
    const z3::expr hi1 = zmin(a + b, 2 * k - a - b);
    const z3::expr lo2 = zabs(c - d);
    const z3::expr hi2 = zmin(c + d, 2 * k - c - d);
    z3::expr result = ctx.int_val(0);
    for (int delta = -target; delta <= target; delta += 2) {
        const int triangle_lo = (target - delta) / 2;
        const z3::expr triangle_hi =
            k - (target + delta) / 2;
        const z3::expr raw_lo = zmax(
            zmax(lo1, lo2 - delta), ctx.int_val(triangle_lo)
        );
        const z3::expr raw_hi = zmin(
            zmin(hi1, hi2 - delta), triangle_hi
        );
        const z3::expr compatible =
            lo1 % 2 == (lo2 - delta) % 2;
        result = result + step_two_count(
            ctx, lo1, raw_lo, raw_hi, compatible
        );
    }
    return result;
}

z3::expr bounded_product(
    z3::context& ctx,
    const z3::expr& small,
    const z3::expr& value,
    const int maximum
) {
    z3::expr result = ctx.int_val(0);
    for (int n = 1; n <= maximum; ++n) {
        result = result + z3::ite(
            small >= n, value, ctx.int_val(0)
        );
    }
    return result;
}

struct QueryResult {
    bool passed = false;
    std::string diagnostic;
};

constexpr std::array<const char*, 4> schema_names{
    "equal-low", "equal-high", "opposite-low", "opposite-high"
};

bool check_priority_witness() {
    constexpr int k = 13;
    constexpr std::array<int, 2> minus{4, 6};
    constexpr std::array<int, 5> plus{8, 8, 8, 12, 12};
    int demand = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                if (!contains_output(k, plus[i], plus[j], minus[o])) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto complement = decomposition(k, block);
                if (complement[0] <= 0 || complement[0] > 2) {
                    continue;
                }
                demand += static_cast<int>(complement[0]);
            }
        }
    }
    const auto selected_local = [&](const std::size_t orientation) {
        const auto active_profile =
            decomposition(k, {minus[orientation], plus[0], plus[3]});
        const auto complement_profile = decomposition(
            k, {minus[1U - orientation], plus[1], plus[2], plus[4]}
        );
        std::int64_t local = 0;
        for (int output : {0, 2, 4}) {
            local += active_profile[static_cast<std::size_t>(output)]
                * complement_profile[static_cast<std::size_t>(output)];
        }
        return local;
    };
    return demand == 27
        && selected_local(0) == 32
        && selected_local(1) == 26;
}

QueryResult verify_query(
    const int schema,
    const int selected_orientation,
    const int parity_mask
) {
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
    for (const z3::expr& plus_label : plus) {
        solver.add(plus_label >= 1 && plus_label <= k - 1);
        solver.add(plus_label != q && plus_label != a);
    }

    const bool selected_high = schema == 1 || schema == 3;
    const bool other_high = schema == 1 || schema == 2;
    solver.add(plus[0] == (selected_high ? k - 1 : ctx.int_val(1)));
    solver.add(plus[4] == (other_high ? k - 1 : ctx.int_val(1)));
    for (std::size_t i = 1; i <= 3; ++i) {
        solver.add(plus[i] >= 2 && plus[i] <= k - 2);
    }
    solver.add(plus[2] <= plus[3]);

    std::vector<z3::expr> active;
    std::vector<z3::expr> rank;
    std::vector<z3::expr> local_bound;
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
                const z3::expr active_triple = fusion(
                    k, plus[static_cast<std::size_t>(i)],
                    plus[static_cast<std::size_t>(j)], r
                );
                const z3::expr is_active =
                    active_triple && complement_nonzero;
                const z3::expr complement_rank = z3::ite(
                    complement_rank_two, ctx.int_val(2),
                    z3::ite(
                        complement_nonzero, ctx.int_val(1), ctx.int_val(0)
                    )
                );
                solver.add(z3::implies(
                    is_active, !complement_rank_three
                ));

                const z3::expr active_0 = z3::ite(
                    active_triple, ctx.int_val(1), ctx.int_val(0)
                );
                const z3::expr active_2 = interval_rank(
                    ctx, k, r, plus[static_cast<std::size_t>(i)],
                    plus[static_cast<std::size_t>(j)], ctx.int_val(2)
                ).simplify();
                const z3::expr active_4 = interval_rank(
                    ctx, k, r, plus[static_cast<std::size_t>(i)],
                    plus[static_cast<std::size_t>(j)], ctx.int_val(4)
                ).simplify();
                solver.add(active_0 >= 0 && active_0 <= 1);
                solver.add(active_2 >= 0 && active_2 <= 3);
                solver.add(active_4 >= 0 && active_4 <= 5);
                const z3::expr complement_2 = four_tensor_output(
                    ctx, k, s,
                    plus[static_cast<std::size_t>(rest[0])],
                    plus[static_cast<std::size_t>(rest[1])],
                    plus[static_cast<std::size_t>(rest[2])], 2
                ).simplify();
                const z3::expr complement_4 = four_tensor_output(
                    ctx, k, s,
                    plus[static_cast<std::size_t>(rest[0])],
                    plus[static_cast<std::size_t>(rest[1])],
                    plus[static_cast<std::size_t>(rest[2])], 4
                ).simplify();
                const z3::expr bound = (
                    bounded_product(ctx, active_0, complement_rank, 1)
                    + bounded_product(ctx, active_2, complement_2, 3)
                    + bounded_product(ctx, active_4, complement_4, 5)
                ).simplify();

                active.push_back(is_active);
                rank.push_back(complement_rank);
                local_bound.push_back(bound);
                orientation.push_back(o);
            }
        }
    }

    const z3::expr zero = ctx.int_val(0);
    const z3::expr one = ctx.int_val(1);
    z3::expr c = zero;
    z3::expr h = zero;
    std::array<z3::expr, 2> color_count{zero, zero};
    for (std::size_t edge = 0; edge < active.size(); ++edge) {
        c = c + z3::ite(active[edge], one, zero);
        h = h + z3::ite(
            active[edge] && rank[edge] == 2, one, zero
        );
        const auto color =
            static_cast<std::size_t>(orientation[edge]);
        color_count[color] = color_count[color]
            + z3::ite(active[edge], one, zero);
    }
    const z3::expr demand = c + h;
    solver.add(color_count[0] > 0 && color_count[1] > 0);
    solver.add(h > 0 && demand >= 12);

    const std::size_t selected_edge =
        static_cast<std::size_t>(selected_orientation);
    solver.add(active[selected_edge] && rank[selected_edge] == 2);

    if (selected_orientation == 1) {
        std::size_t edge = 0;
        for (std::size_t i = 0; i < plus.size(); ++i) {
            for (std::size_t j = i + 1; j < plus.size(); ++j) {
                const z3::expr cap_i =
                    plus[i] == 1 || plus[i] == k - 1;
                const z3::expr cap_j =
                    plus[j] == 1 || plus[j] == k - 1;
                const z3::expr incidence_one =
                    (cap_i && !cap_j) || (!cap_i && cap_j);
                solver.add(z3::implies(
                    incidence_one,
                    !(active[edge] && rank[edge] == 2)
                ));
                edge += 2U;
            }
        }
    }

    solver.add(local_bound[selected_edge] < demand);
    const z3::check_result result = solver.check();
    if (result == z3::unsat) {
        return {true, {}};
    }
    std::ostringstream diagnostic;
    diagnostic << "schema=" << schema_names[
        static_cast<std::size_t>(schema)
    ] << " orientation=" << (selected_orientation == 0 ? 'q' : 'a')
      << " mask=" << parity_mask << " result=" << result;
    if (result == z3::sat) {
        diagnostic << " model=" << solver.get_model();
    } else {
        diagnostic << " reason=" << solver.reason_unknown();
    }
    return {false, diagnostic.str()};
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
    const std::int64_t small_cases = check_small_levels();
    if (small_cases < 0 || !check_priority_witness()) {
        std::cerr << "SU2_D12_MIXED_N2_Z3 FAIL base_witness\n";
        return EXIT_FAILURE;
    }

    if (argc == 4) {
        const int schema = std::atoi(argv[1]);
        const int orientation = std::atoi(argv[2]);
        const int parity_mask = std::atoi(argv[3]);
        if (schema < 0 || schema >= 4
            || orientation < 0 || orientation >= 2
            || parity_mask < 0 || parity_mask >= 256) {
            std::cerr
                << "usage: verify_su2_d12_mixed_n2_z3"
                << " [SCHEMA(0..3) ORIENTATION(0..1) MASK(0..255)]\n";
            return EXIT_FAILURE;
        }
        const QueryResult result =
            verify_query(schema, orientation, parity_mask);
        if (!result.passed) {
            std::cerr << "SU2_D12_MIXED_N2_Z3 FAIL "
                      << result.diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout
            << "SU2_D12_MIXED_N2_Z3 schema="
            << schema_names[static_cast<std::size_t>(schema)]
            << " orientation=" << (orientation == 0 ? 'q' : 'a')
            << " mask=" << parity_mask
            << " local_counterexample=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    }
    if (argc != 1) {
        std::cerr
            << "usage: verify_su2_d12_mixed_n2_z3"
            << " [SCHEMA(0..3) ORIENTATION(0..1) MASK(0..255)]\n";
        return EXIT_FAILURE;
    }

    constexpr int task_count = 4 * 2 * 256;
    const unsigned workers = worker_limit();
    std::atomic<int> next_task{0};
    std::atomic<int> completed{0};
    std::atomic<bool> failed{false};
    std::mutex diagnostic_mutex;
    std::string diagnostic;
    std::vector<std::jthread> pool;
    pool.reserve(workers);
    for (unsigned worker = 0; worker < workers; ++worker) {
        pool.emplace_back([&]() {
            while (!failed.load()) {
                const int task = next_task.fetch_add(1);
                if (task >= task_count) {
                    return;
                }
                const int parity_mask = task % 256;
                const int orientation = (task / 256) % 2;
                const int schema = task / (2 * 256);
                const QueryResult result =
                    verify_query(schema, orientation, parity_mask);
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

    if (failed.load() || completed.load() != task_count) {
        std::cerr << "SU2_D12_MIXED_N2_Z3 FAIL"
                  << " completed=" << completed.load()
                  << " diagnostic=" << diagnostic << '\n';
        return EXIT_FAILURE;
    }
    std::cout
        << "SU2_D12_MIXED_N2_Z3 k>=9"
        << " schemas=4 orientations=q-first"
        << " masks=256 queries=" << task_count
        << " workers=" << workers
        << " transition_cases=" << small_cases
        << " priority_witness=T27,q32,a26"
        << " local_counterexamples=UNSAT result=PASS\n";
    return EXIT_SUCCESS;
}
