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

using Cut = std::array<int, 3>;

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

std::vector<Cut> negative_cuts() {
    std::vector<Cut> result;
    for (int minus = 0; minus < 4; ++minus) {
        for (int first_plus = 4; first_plus < 7; ++first_plus) {
            for (int second_plus = first_plus + 1;
                 second_plus < 7; ++second_plus) {
                result.push_back(Cut{minus, first_plus, second_plus});
            }
        }
    }
    for (int omitted_minus = 0; omitted_minus < 4; ++omitted_minus) {
        Cut cut{};
        int next = 0;
        for (int minus = 0; minus < 4; ++minus) {
            if (minus != omitted_minus) {
                cut[static_cast<std::size_t>(next)] = minus;
                ++next;
            }
        }
        result.push_back(cut);
    }
    return result;
}

std::array<int, 4> complement(const Cut& cut) {
    std::array<int, 4> result{};
    int next = 0;
    for (int index = 0; index < 7; ++index) {
        if (std::find(cut.begin(), cut.end(), index) == cut.end()) {
            result[static_cast<std::size_t>(next)] = index;
            ++next;
        }
    }
    return result;
}

z3::expr triple_output(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const Cut& cut,
    const int output,
    const bool reflect
) {
    const z3::expr first = reflect
        ? k - labels[static_cast<std::size_t>(cut[0])]
        : labels[static_cast<std::size_t>(cut[0])];
    return interval_rank(
        ctx,
        k,
        first,
        labels[static_cast<std::size_t>(cut[1])],
        labels[static_cast<std::size_t>(cut[2])],
        ctx.int_val(output)
    );
}

z3::expr fourfold_output(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const std::array<int, 4>& rest,
    const int output,
    const bool reflect
) {
    const z3::expr first = reflect
        ? k - labels[static_cast<std::size_t>(rest[0])]
        : labels[static_cast<std::size_t>(rest[0])];
    return four_tensor_output(
        ctx,
        k,
        first,
        labels[static_cast<std::size_t>(rest[1])],
        labels[static_cast<std::size_t>(rest[2])],
        labels[static_cast<std::size_t>(rest[3])],
        output
    );
}

z3::expr selected_local(
    z3::context& ctx,
    const z3::expr& k,
    const z3::expr& shallow,
    const std::array<z3::expr, 7>& labels,
    const Cut& cut,
    const int level_parity
) {
    const std::array<int, 4> rest = complement(cut);
    z3::expr result = ctx.int_val(0);
    for (int j = 0; j <= 4; ++j) {
        const int low = 2 * j;
        result = result + bounded_product(
            ctx,
            triple_output(ctx, k, labels, cut, low, false),
            fourfold_output(ctx, k, labels, rest, low, false),
            low + 1
        );

        const int reflected_low = level_parity + 2 * j;
        const z3::expr high = shallow - 2 * j;
        const z3::expr reflected_term = bounded_product(
            ctx,
            triple_output(
                ctx, k, labels, cut, reflected_low, true
            ),
            fourfold_output(
                ctx, k, labels, rest, reflected_low, true
            ),
            reflected_low + 1
        );
        result = result + z3::ite(
            high > 8, reflected_term, ctx.int_val(0)
        );
    }
    return result;
}

z3::expr pair_reservoir(
    z3::context& ctx,
    const z3::expr& k,
    const z3::expr& shallow,
    const std::array<z3::expr, 7>& labels,
    const int level_parity
) {
    z3::expr result = ctx.int_val(0);
    const auto add_pair = [&](const int first, const int second) {
        std::array<int, 5> rest{};
        int next = 0;
        for (int index = 0; index < 7; ++index) {
            if (index != first && index != second) {
                rest[static_cast<std::size_t>(next)] = index;
                ++next;
            }
        }
        z3::expr contribution = ctx.int_val(0);
        for (int j = 0; j <= 4; ++j) {
            const int low = 2 * j;
            const z3::expr low_pair = fusion(
                k,
                labels[static_cast<std::size_t>(rest[0])],
                labels[static_cast<std::size_t>(rest[1])],
                ctx.int_val(low)
            );
            const z3::expr low_triple = interval_rank(
                ctx,
                k,
                labels[static_cast<std::size_t>(rest[2])],
                labels[static_cast<std::size_t>(rest[3])],
                labels[static_cast<std::size_t>(rest[4])],
                ctx.int_val(low)
            );
            contribution = contribution + z3::ite(
                low_pair, low_triple, ctx.int_val(0)
            );

            const int reflected_low = level_parity + 2 * j;
            const z3::expr high = shallow - 2 * j;
            const z3::expr high_pair = fusion(
                k,
                k - labels[static_cast<std::size_t>(rest[0])],
                labels[static_cast<std::size_t>(rest[1])],
                ctx.int_val(reflected_low)
            );
            const z3::expr high_triple = interval_rank(
                ctx,
                k,
                k - labels[static_cast<std::size_t>(rest[2])],
                labels[static_cast<std::size_t>(rest[3])],
                labels[static_cast<std::size_t>(rest[4])],
                ctx.int_val(reflected_low)
            );
            contribution = contribution + z3::ite(
                high > 8 && high_pair,
                high_triple,
                ctx.int_val(0)
            );
        }
        result = result + z3::ite(
            labels[static_cast<std::size_t>(first)]
                == labels[static_cast<std::size_t>(second)],
            contribution,
            ctx.int_val(0)
        );
    };

    for (int first = 0; first < 4; ++first) {
        for (int second = first + 1; second < 4; ++second) {
            add_pair(first, second);
        }
    }
    for (int first = 4; first < 7; ++first) {
        for (int second = first + 1; second < 7; ++second) {
            add_pair(first, second);
        }
    }
    return result;
}

struct QueryResult {
    bool passed = false;
    std::string diagnostic;
};

QueryResult verify_query(
    const int selected_cut,
    const int level_parity,
    const int shallow_orbit
) {
    z3::context ctx;
    z3::solver solver(ctx);

    const z3::expr k =
        2 * ctx.int_const("k_half") + level_parity;
    const z3::expr shallow = k - level_parity;
    const std::array<z3::expr, 7> labels{
        2 * ctx.int_const("m0_half"),
        2 * ctx.int_const("m1_half"),
        2 * ctx.int_const("m2_half"),
        2 * ctx.int_const("m3_half"),
        2 * ctx.int_const("p0_half"),
        2 * ctx.int_const("p1_half"),
        2 * ctx.int_const("p2_half")
    };

    solver.add(k >= 4);
    for (const z3::expr& label : labels) {
        solver.add(label >= 2 && label <= k);
    }
    if (selected_cut == 0) {
        // Canonical one-minus/two-plus cut {0,4,5}.  Sort only inside
        // the selected and complementary sign-class blocks.
        solver.add(labels[1] <= labels[2]);
        solver.add(labels[2] <= labels[3]);
        solver.add(labels[4] <= labels[5]);
        constexpr std::array<int, 4> representatives{0, 3, 5, 6};
        solver.add(
            labels[static_cast<std::size_t>(
                representatives[static_cast<std::size_t>(shallow_orbit)]
            )] == shallow
        );
    } else {
        // Canonical three-minus cut {1,2,3}, whose omitted minus is 0.
        solver.add(labels[1] <= labels[2]);
        solver.add(labels[2] <= labels[3]);
        solver.add(labels[4] <= labels[5]);
        solver.add(labels[5] <= labels[6]);
        constexpr std::array<int, 3> representatives{3, 0, 6};
        solver.add(
            labels[static_cast<std::size_t>(
                representatives[static_cast<std::size_t>(shallow_orbit)]
            )] == shallow
        );
    }
    for (std::size_t minus = 0; minus < 4U; ++minus) {
        for (std::size_t plus = 4U; plus < 7U; ++plus) {
            solver.add(labels[minus] != labels[plus]);
        }
    }

    const std::vector<Cut> cuts = negative_cuts();
    std::vector<z3::expr> ranks;
    ranks.reserve(cuts.size());
    z3::expr demand = ctx.int_val(0);
    for (const Cut& cut : cuts) {
        const std::array<int, 4> rest = complement(cut);
        const z3::expr active = fusion(
            k,
            labels[static_cast<std::size_t>(cut[0])],
            labels[static_cast<std::size_t>(cut[1])],
            labels[static_cast<std::size_t>(cut[2])]
        );
        const z3::expr rank = z3::ite(
            active,
            interval_rank(
                ctx,
                k,
                labels[static_cast<std::size_t>(rest[0])],
                labels[static_cast<std::size_t>(rest[1])],
                labels[static_cast<std::size_t>(rest[2])],
                labels[static_cast<std::size_t>(rest[3])]
            ),
            ctx.int_val(0)
        );
        ranks.push_back(rank);
        demand = demand + rank;
    }

    const std::size_t selected =
        static_cast<std::size_t>(selected_cut);
    solver.add(ranks[selected] > 0);
    for (const z3::expr& rank : ranks) {
        solver.add(rank <= ranks[selected]);
    }

    const z3::expr local = selected_local(
        ctx, k, shallow, labels, cuts[selected], level_parity
    );
    const z3::expr reservoir = pair_reservoir(
        ctx, k, shallow, labels, level_parity
    );
    solver.add(local + reservoir < demand);

    const z3::check_result result = solver.check();
    if (result == z3::unsat) {
        return {true, {}};
    }
    std::ostringstream diagnostic;
    diagnostic << "cut=" << selected_cut
               << " parity=" << level_parity
               << " shallow_orbit=" << shallow_orbit
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
    constexpr std::uint64_t kib_per_worker = 256U * 1024U;
    const unsigned memory_limit = available_kib == 0
        ? hardware
        : static_cast<unsigned>(std::max<std::uint64_t>(
            1U, available_kib / kib_per_worker
        ));
    return std::max(1U, std::min(cpu_limit, memory_limit));
}

}  // namespace

int main(int argc, char** argv) {
    constexpr int parity_count = 2;
    constexpr int shallow_orbits = 7;
    constexpr int task_count = shallow_orbits * parity_count;

    if (argc == 4) {
        const int cut = std::atoi(argv[1]);
        const int parity = std::atoi(argv[2]);
        const int shallow_orbit = std::atoi(argv[3]);
        if ((cut != 0 && cut != 12)
            || parity < 0 || parity >= parity_count
            || shallow_orbit < 0
            || (cut == 0 && shallow_orbit >= 4)
            || (cut == 12 && shallow_orbit >= 3)) {
            std::cerr
                << "usage: verify_su2_four_minus_shallow_z3 "
                << "[CUT(0|12) LEVEL_PARITY(0..1) "
                << "SHALLOW_ORBIT(0..3|0..2)]\n";
            return EXIT_FAILURE;
        }
        const QueryResult result = verify_query(
            cut, parity, shallow_orbit
        );
        if (!result.passed) {
            std::cerr << "SU2_FOUR_MINUS_SHALLOW_Z3 FAIL "
                      << result.diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_FOUR_MINUS_SHALLOW_Z3"
                  << " cut=" << cut
                  << " parity=" << parity
                  << " shallow_orbit=" << shallow_orbit
                  << " counterexample=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    }
    if (argc != 1) {
        std::cerr
            << "usage: verify_su2_four_minus_shallow_z3 "
            << "[CUT(0|12) LEVEL_PARITY(0..1) "
            << "SHALLOW_ORBIT(0..3|0..2)]\n";
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
                const int parity = task % parity_count;
                const int orbit = task / parity_count;
                const int cut = orbit < 4 ? 0 : 12;
                const int shallow_orbit =
                    cut == 0 ? orbit : orbit - 4;
                const QueryResult result = verify_query(
                    cut, parity, shallow_orbit
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
        std::cerr << "SU2_FOUR_MINUS_SHALLOW_Z3 FAIL"
                  << " completed=" << completed.load()
                  << " diagnostic=" << diagnostic << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "SU2_FOUR_MINUS_SHALLOW_Z3"
              << " parities=2 selected_orbits=2 shallow_orbits=7"
              << " queries=" << task_count
              << " workers=" << workers
              << " counterexamples=UNSAT result=PASS\n";
    return EXIT_SUCCESS;
}
