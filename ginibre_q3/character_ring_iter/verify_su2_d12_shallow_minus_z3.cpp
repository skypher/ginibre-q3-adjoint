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

using IndexTriple = std::array<int, 3>;

int bit_count(int value) {
    int result = 0;
    while (value != 0) {
        result += value & 1;
        value >>= 1;
    }
    return result;
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
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& output
) {
    return zabs(first - second) <= output
        && output <= first + second
        && output <= 2 * k - first - second
        && (first + second + output) % 2 == 0;
}

z3::expr interval_rank(
    z3::context& ctx,
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& c,
    const z3::expr& d
) {
    const z3::expr lower_first = zabs(a - b);
    const z3::expr upper_first = zmin(a + b, 2 * k - a - b);
    const z3::expr lower_second = zabs(c - d);
    const z3::expr upper_second = zmin(c + d, 2 * k - c - d);
    const z3::expr lower = zmax(lower_first, lower_second);
    const z3::expr upper = zmin(upper_first, upper_second);
    return z3::ite(
        lower_first % 2 == lower_second % 2 && lower <= upper,
        (upper - lower) / 2 + 1,
        ctx.int_val(0)
    );
}

z3::expr step_two_count(
    z3::context& ctx,
    const z3::expr& base,
    const z3::expr& raw_lower,
    const z3::expr& raw_upper,
    const z3::expr& compatible
) {
    const z3::expr lower = z3::ite(
        (raw_lower - base) % 2 == 0, raw_lower, raw_lower + 1
    );
    const z3::expr upper = z3::ite(
        (raw_upper - base) % 2 == 0, raw_upper, raw_upper - 1
    );
    return z3::ite(
        compatible && lower <= upper,
        (upper - lower) / 2 + 1,
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
    const int output
) {
    const z3::expr lower_first = zabs(a - b);
    const z3::expr upper_first = zmin(a + b, 2 * k - a - b);
    const z3::expr lower_second = zabs(c - d);
    const z3::expr upper_second = zmin(c + d, 2 * k - c - d);
    z3::expr result = ctx.int_val(0);
    for (int delta = -output; delta <= output; delta += 2) {
        const z3::expr raw_lower = zmax(
            zmax(lower_first, lower_second - delta),
            ctx.int_val((output - delta) / 2)
        );
        const z3::expr raw_upper = zmin(
            zmin(upper_first, upper_second - delta),
            k - (output + delta) / 2
        );
        result = result + step_two_count(
            ctx,
            lower_first,
            raw_lower,
            raw_upper,
            lower_first % 2 == (lower_second - delta) % 2
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

IndexTriple pair_complement(const int first, const int second) {
    IndexTriple result{};
    int next = 0;
    for (int index = 0; index < 5; ++index) {
        if (index != first && index != second) {
            result[static_cast<std::size_t>(next)] = index;
            ++next;
        }
    }
    return result;
}

z3::expr triple_output(
    z3::context& ctx,
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& c,
    const z3::expr& output
) {
    return interval_rank(ctx, k, a, b, c, output);
}

z3::expr negative_rank(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 2>& minus,
    const std::array<z3::expr, 5>& plus,
    const int orientation,
    const int first,
    const int second
) {
    const IndexTriple rest = pair_complement(first, second);
    const z3::expr active = fusion(
        k,
        plus[static_cast<std::size_t>(first)],
        plus[static_cast<std::size_t>(second)],
        minus[static_cast<std::size_t>(orientation)]
    );
    const z3::expr complement = interval_rank(
        ctx,
        k,
        minus[static_cast<std::size_t>(1 - orientation)],
        plus[static_cast<std::size_t>(rest[0])],
        plus[static_cast<std::size_t>(rest[1])],
        plus[static_cast<std::size_t>(rest[2])]
    );
    return z3::ite(active, complement, ctx.int_val(0));
}

z3::expr local_cut(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 2>& minus,
    const std::array<z3::expr, 5>& plus,
    const int orientation,
    const int first,
    const int second
) {
    const IndexTriple rest = pair_complement(first, second);
    z3::expr result = ctx.int_val(0);
    for (int output = 0; output <= 4; output += 2) {
        const z3::expr active_profile = triple_output(
            ctx,
            k,
            minus[static_cast<std::size_t>(orientation)],
            plus[static_cast<std::size_t>(first)],
            plus[static_cast<std::size_t>(second)],
            ctx.int_val(output)
        );
        const z3::expr complement_profile = four_tensor_output(
            ctx,
            k,
            minus[static_cast<std::size_t>(1 - orientation)],
            plus[static_cast<std::size_t>(rest[0])],
            plus[static_cast<std::size_t>(rest[1])],
            plus[static_cast<std::size_t>(rest[2])],
            output
        );
        result = result + bounded_product(
            ctx, active_profile, complement_profile, output + 1
        );
    }
    return result;
}

z3::expr shallow_fivefold(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 2>& minus,
    const std::array<z3::expr, 5>& plus,
    const IndexTriple& rest
) {
    const z3::expr lower = zabs(minus[0] - minus[1]);
    const z3::expr upper = zmin(
        minus[0] + minus[1], 2 * k - minus[0] - minus[1]
    );
    const z3::expr first = triple_output(
        ctx,
        k,
        plus[static_cast<std::size_t>(rest[0])],
        plus[static_cast<std::size_t>(rest[1])],
        plus[static_cast<std::size_t>(rest[2])],
        lower
    );
    const z3::expr second = triple_output(
        ctx,
        k,
        plus[static_cast<std::size_t>(rest[0])],
        plus[static_cast<std::size_t>(rest[1])],
        plus[static_cast<std::size_t>(rest[2])],
        lower + 2
    );
    return first + z3::ite(
        upper >= lower + 2, second, ctx.int_val(0)
    );
}

z3::expr equal_plus_reservoir(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 2>& minus,
    const std::array<z3::expr, 5>& plus
) {
    z3::expr result = ctx.int_val(0);
    for (int first = 0; first < 5; ++first) {
        for (int second = first + 1; second < 5; ++second) {
            const IndexTriple rest = pair_complement(first, second);
            const z3::expr contribution = shallow_fivefold(
                ctx, k, minus, plus, rest
            );
            result = result + z3::ite(
                plus[static_cast<std::size_t>(first)]
                    == plus[static_cast<std::size_t>(second)],
                contribution,
                ctx.int_val(0)
            );
        }
    }
    return result;
}

z3::expr positive_reservoir(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 2>& minus,
    const std::array<z3::expr, 5>& plus
) {
    z3::expr result = ctx.int_val(0);
    for (int selected = 0; selected < 5; ++selected) {
        std::array<int, 4> rest{};
        int next = 0;
        for (int index = 0; index < 5; ++index) {
            if (index != selected) {
                rest[static_cast<std::size_t>(next)] = index;
                ++next;
            }
        }
        const z3::expr active = fusion(
            k,
            minus[0],
            minus[1],
            plus[static_cast<std::size_t>(selected)]
        );
        const z3::expr complement = interval_rank(
            ctx,
            k,
            plus[static_cast<std::size_t>(rest[0])],
            plus[static_cast<std::size_t>(rest[1])],
            plus[static_cast<std::size_t>(rest[2])],
            plus[static_cast<std::size_t>(rest[3])]
        );
        result = result + z3::ite(
            active, complement, ctx.int_val(0)
        );
    }

    for (int first = 0; first < 5; ++first) {
        for (int second = first + 1; second < 5; ++second) {
            const IndexTriple rest = pair_complement(first, second);
            const z3::expr selected = interval_rank(
                ctx,
                k,
                minus[0],
                minus[1],
                plus[static_cast<std::size_t>(first)],
                plus[static_cast<std::size_t>(second)]
            );
            const z3::expr complement = fusion(
                k,
                plus[static_cast<std::size_t>(rest[0])],
                plus[static_cast<std::size_t>(rest[1])],
                plus[static_cast<std::size_t>(rest[2])]
            );
            result = result + z3::ite(
                complement, selected, ctx.int_val(0)
            );
        }
    }
    return result;
}

struct QueryResult {
    bool passed = false;
    std::string diagnostic;
};

QueryResult verify_query(
    const int shallow_type,
    const int level_parity,
    const int plus_parity_mask,
    const int selected_cut
) {
    z3::context ctx;
    z3::solver solver(ctx);

    const z3::expr k =
        2 * ctx.int_const("k_half") + level_parity;
    const int plus_parity = plus_parity_mask >= 0
        ? bit_count(plus_parity_mask) % 2
        : -1;
    if (shallow_type == 3 && plus_parity >= 0
        && plus_parity != 1) {
        return {true, {}};
    }
    if (shallow_type >= 4 && plus_parity >= 0
        && plus_parity != 0) {
        return {true, {}};
    }
    const auto parity_variable = [&ctx](
        const char* name, const int parity
    ) {
        if (parity < 0) {
            return z3::expr(ctx.int_const(name));
        }
        const std::string half_name = std::string(name) + "_half";
        return z3::expr(
            2 * ctx.int_const(half_name.c_str()) + parity
        );
    };

    z3::expr q = ctx.int_const("q");
    z3::expr a = ctx.int_const("a");
    if (shallow_type == 0) {
        q = ctx.int_val(1);
        a = parity_variable(
            "a", plus_parity < 0 ? -1 : (plus_parity + 1) % 2
        );
    } else if (shallow_type == 1) {
        q = parity_variable(
            "q",
            plus_parity < 0
                ? -1
                : (plus_parity + 1 - level_parity + 2) % 2
        );
        a = k - 1;
    } else if (shallow_type == 2) {
        q = parity_variable(
            "q",
            plus_parity < 0
                ? -1
                : (plus_parity + level_parity) % 2
        );
        a = k;
    } else if (shallow_type == 3) {
        q = k - 1;
        a = k;
    } else if (shallow_type == 4) {
        q = ctx.int_val(1);
        a = ctx.int_val(1);
    } else if (shallow_type == 5) {
        q = k - 1;
        a = k - 1;
    } else {
        q = k;
        a = k;
    }
    const std::array<z3::expr, 2> minus{q, a};
    const auto plus_variable = [&](
        const char* name, const int index
    ) {
        const int parity = plus_parity_mask < 0
            ? -1
            : (plus_parity_mask >> index) & 1;
        return parity_variable(name, parity);
    };
    const std::array<z3::expr, 5> plus{
        plus_variable("p0", 0),
        plus_variable("p1", 1),
        plus_variable("p2", 2),
        plus_variable("p3", 3),
        plus_variable("p4", 4)
    };

    solver.add(k >= 3);
    solver.add(minus[0] >= 1 && minus[1] <= k);
    if (shallow_type < 4) {
        solver.add(minus[0] < minus[1]);
    } else {
        solver.add(minus[0] == minus[1]);
    }
    if (shallow_type == 1) {
        solver.add(minus[0] >= 2);
    } else if (shallow_type == 2) {
        solver.add(minus[0] >= 2 && minus[0] <= k - 2);
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
        solver.add(plus[i] >= 1 && plus[i] <= k);
        solver.add(plus[i] != minus[0] && plus[i] != minus[1]);
        if (selected_cut < 0 && i + 1U < plus.size()) {
            solver.add(plus[i] <= plus[i + 1U]);
        }
    }
    if (selected_cut >= 0) {
        solver.add(plus[0] <= plus[1]);
        solver.add(plus[2] <= plus[3]);
        solver.add(plus[3] <= plus[4]);
    }

    std::vector<z3::expr> ranks;
    std::vector<z3::expr> locals;
    z3::expr demand = ctx.int_val(0);
    for (int orientation = 0; orientation < 2; ++orientation) {
        for (int first = 0; first < 5; ++first) {
            for (int second = first + 1; second < 5; ++second) {
                const z3::expr rank = negative_rank(
                    ctx,
                    k,
                    minus,
                    plus,
                    orientation,
                    first,
                    second
                );
                solver.add(rank <= 2);
                const int cut_index = static_cast<int>(ranks.size());
                ranks.push_back(rank);
                locals.push_back(
                    selected_cut >= 0 && cut_index != selected_cut
                        ? ctx.int_val(0)
                        : local_cut(
                            ctx,
                            k,
                            minus,
                            plus,
                            orientation,
                            first,
                            second
                        )
                );
                demand = demand + rank;
            }
        }
    }
    solver.add(demand > 0);

    z3::expr supply = equal_plus_reservoir(
        ctx, k, minus, plus
    ) + positive_reservoir(ctx, k, minus, plus);
    if (shallow_type >= 4) {
        supply = supply + 1;
    }
    if (selected_cut >= 0) {
        const std::size_t selected =
            static_cast<std::size_t>(selected_cut);
        solver.add(ranks[selected] > 0);
        for (std::size_t i = 0; i < ranks.size(); ++i) {
            solver.add(ranks[i] <= ranks[selected]);
        }
        solver.add(locals[selected] + supply < demand);
    } else {
        for (std::size_t i = 0; i < ranks.size(); ++i) {
            solver.add(z3::implies(
                ranks[i] > 0, locals[i] + supply < demand
            ));
        }
    }

    const z3::check_result result = solver.check();
    if (result == z3::unsat) {
        return {true, {}};
    }
    std::ostringstream diagnostic;
    diagnostic << "type=" << shallow_type
               << " parity=" << level_parity
               << " plus_mask=" << plus_parity_mask
               << " selected_cut=" << selected_cut
               << " result=" << result;
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
    constexpr std::uint64_t kib_per_worker = 1024U * 1024U;
    const unsigned memory_limit = available_kib == 0
        ? hardware
        : static_cast<unsigned>(std::max<std::uint64_t>(
            1U, available_kib / kib_per_worker
        ));
    return std::max(1U, std::min(cpu_limit, memory_limit));
}

}  // namespace

int main(int argc, char** argv) {
    constexpr int shallow_types = 7;
    constexpr int parity_count = 2;
    constexpr int plus_masks = 32;
    constexpr int task_count =
        shallow_types * parity_count * plus_masks * 2;

    if (argc == 3) {
        const int type = std::atoi(argv[1]);
        const int parity = std::atoi(argv[2]);
        if (type < 0 || type >= shallow_types
            || parity < 0 || parity >= parity_count) {
            std::cerr
                << "usage: verify_su2_d12_shallow_minus_z3 "
                << "[TYPE(0..6) LEVEL_PARITY(0..1)]\n";
            return EXIT_FAILURE;
        }
        const QueryResult result = verify_query(type, parity, -1, -1);
        if (!result.passed) {
            std::cerr << "SU2_D12_SHALLOW_MINUS_Z3 FAIL "
                      << result.diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_D12_SHALLOW_MINUS_Z3"
                  << " type=" << type
                  << " parity=" << parity
                  << " counterexample=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    }
    if (argc == 4) {
        const int type = std::atoi(argv[1]);
        const int parity = std::atoi(argv[2]);
        const int plus_mask = std::atoi(argv[3]);
        if (type < 0 || type >= shallow_types
            || parity < 0 || parity >= parity_count
            || plus_mask < 0 || plus_mask >= plus_masks) {
            std::cerr
                << "usage: verify_su2_d12_shallow_minus_z3 "
                << "[TYPE(0..6) LEVEL_PARITY(0..1) [PLUS_MASK(0..31)]]\n";
            return EXIT_FAILURE;
        }
        const QueryResult result = verify_query(
            type, parity, plus_mask, -1
        );
        if (!result.passed) {
            std::cerr << "SU2_D12_SHALLOW_MINUS_Z3 FAIL "
                      << result.diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_D12_SHALLOW_MINUS_Z3"
                  << " type=" << type
                  << " parity=" << parity
                  << " plus_mask=" << plus_mask
                  << " counterexample=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    }
    if (argc == 5) {
        const int type = std::atoi(argv[1]);
        const int parity = std::atoi(argv[2]);
        const int plus_mask = std::atoi(argv[3]);
        const int selected_cut = std::atoi(argv[4]);
        if (type < 0 || type >= shallow_types
            || parity < 0 || parity >= parity_count
            || plus_mask < 0 || plus_mask >= plus_masks
            || (selected_cut != 0 && selected_cut != 10)) {
            std::cerr
                << "usage: verify_su2_d12_shallow_minus_z3 "
                << "[TYPE(0..6) LEVEL_PARITY(0..1) "
                << "PLUS_MASK(0..31) SELECTED_CUT(0|10)]\n";
            return EXIT_FAILURE;
        }
        const QueryResult result = verify_query(
            type, parity, plus_mask, selected_cut
        );
        if (!result.passed) {
            std::cerr << "SU2_D12_SHALLOW_MINUS_Z3 FAIL "
                      << result.diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_D12_SHALLOW_MINUS_Z3"
                  << " type=" << type
                  << " parity=" << parity
                  << " plus_mask=" << plus_mask
                  << " selected_cut=" << selected_cut
                  << " counterexample=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    }
    if (argc != 1) {
        std::cerr
            << "usage: verify_su2_d12_shallow_minus_z3 "
            << "[TYPE(0..6) LEVEL_PARITY(0..1) "
            << "[PLUS_MASK(0..31) [SELECTED_CUT(0|10)]]]\n";
        return EXIT_FAILURE;
    }

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
                const int selected_cut =
                    task % 2 == 0 ? 0 : 10;
                const int plus_mask =
                    (task / 2) % plus_masks;
                const int parity =
                    (task / (2 * plus_masks)) % parity_count;
                const int type =
                    task / (2 * plus_masks * parity_count);
                const QueryResult result = verify_query(
                    type, parity, plus_mask, selected_cut
                );
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
        std::cerr << "SU2_D12_SHALLOW_MINUS_Z3 FAIL"
                  << " completed=" << completed.load()
                  << " diagnostic=" << diagnostic << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "SU2_D12_SHALLOW_MINUS_Z3"
              << " scope=all types=7 parities=2 plus_masks=32"
              << " selected_orbits=2 queries="
              << task_count
              << " workers=" << workers
              << " counterexamples=UNSAT result=PASS\n";
    return EXIT_SUCCESS;
}
