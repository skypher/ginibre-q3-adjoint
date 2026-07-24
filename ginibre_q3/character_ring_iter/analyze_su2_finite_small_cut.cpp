#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Labels2 = std::array<int, 2>;
using Labels5 = std::array<int, 5>;

bool use_cache = false;
thread_local int cached_level = -1;
thread_local std::unordered_map<std::uint64_t, std::int64_t> cache;

template <class Function>
void fusion_outputs(int k, int a, int b, Function function) {
    const int upper = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= upper; c += 2) {
        function(c);
    }
}

std::int64_t multiplicity(int k, const std::vector<int>& labels) {
    std::vector<int> ordered = labels;
    std::sort(ordered.begin(), ordered.end());
    std::uint64_t key = ordered.size();
    const auto base = static_cast<std::uint64_t>(k + 2);
    for (const int label : ordered) {
        key = key * base + static_cast<std::uint64_t>(label + 1);
    }
    if (use_cache) {
        if (cached_level != k) {
            cache.clear();
            cached_level = k;
        }
        const auto found = cache.find(key);
        if (found != cache.end()) {
            return found->second;
        }
    }
    std::array<std::int64_t, 256> current{};
    std::array<std::int64_t, 256> next{};
    if (k >= static_cast<int>(current.size())) {
        throw std::runtime_error("level exceeds internal bound");
    }
    current[0] = 1;
    for (const int label : labels) {
        next.fill(0);
        for (int source = 0; source <= k; ++source) {
            const auto value = current[static_cast<std::size_t>(source)];
            if (value == 0) {
                continue;
            }
            fusion_outputs(k, source, label, [&](int output) {
                next[static_cast<std::size_t>(output)] += value;
            });
        }
        current = next;
    }
    if (use_cache) {
        cache.emplace(key, current[0]);
    }
    return current[0];
}

bool disjoint(const Labels2& minus, const Labels5& plus) {
    for (const int q : minus) {
        for (const int p : plus) {
            if (q == p) {
                return false;
            }
        }
    }
    return true;
}

struct Statistics {
    std::int64_t n = 0;
    std::int64_t p_equal = 0;
    std::int64_t u_one = 0;
    std::int64_t u_two = 0;
    std::int64_t t = 0;
    int active = 0;
    int d = 0;
};

Statistics statistics(int k, const Labels2& minus, const Labels5& plus) {
    Statistics result;
    result.n = multiplicity(
        k, {minus[0], minus[1], plus[0], plus[1], plus[2], plus[3], plus[4]}
    );
    if (minus[0] == minus[1]) {
        result.p_equal += multiplicity(
            k, {plus[0], plus[1], plus[2], plus[3], plus[4]}
        );
    }
    for (int i = 0; i < 5; ++i) {
        std::vector<int> rest;
        for (int j = 0; j < 5; ++j) {
            if (j != i) {
                rest.push_back(plus[static_cast<std::size_t>(j)]);
            }
        }
        result.u_one += multiplicity(
            k, {minus[0], minus[1], plus[static_cast<std::size_t>(i)]}
        ) * multiplicity(k, rest);
    }
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            std::vector<int> rest{minus[0], minus[1]};
            std::vector<int> complement;
            for (int h = 0; h < 5; ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[static_cast<std::size_t>(h)]);
                    complement.push_back(plus[static_cast<std::size_t>(h)]);
                }
            }
            if (plus[static_cast<std::size_t>(i)]
                == plus[static_cast<std::size_t>(j)]) {
                result.p_equal += multiplicity(k, rest);
            }
            result.u_two += multiplicity(
                k,
                {
                    minus[0], minus[1],
                    plus[static_cast<std::size_t>(i)],
                    plus[static_cast<std::size_t>(j)]
                }
            ) * multiplicity(k, complement);
            for (int orientation = 0; orientation < 2; ++orientation) {
                const auto triple = multiplicity(
                    k,
                    {
                        minus[static_cast<std::size_t>(orientation)],
                        plus[static_cast<std::size_t>(i)],
                        plus[static_cast<std::size_t>(j)]
                    }
                );
                if (triple == 0) {
                    continue;
                }
                std::vector<int> block{
                    minus[static_cast<std::size_t>(1 - orientation)]
                };
                block.insert(block.end(), complement.begin(), complement.end());
                const auto cut = triple * multiplicity(k, block);
                if (cut > 0) {
                    ++result.active;
                    result.d = std::max(result.d, static_cast<int>(cut));
                    result.t += cut;
                }
            }
        }
    }
    return result;
}

struct Witness {
    bool set = false;
    std::int64_t value = 0;
    int k = 0;
    Labels2 minus{};
    Labels5 plus{};
    Statistics stats{};
};

void consider(
    Witness& witness,
    std::int64_t value,
    int k,
    const Labels2& minus,
    const Labels5& plus,
    const Statistics& stats
) {
    if (!witness.set || value < witness.value) {
        witness = {true, value, k, minus, plus, stats};
    }
}

void print(const char* name, const Witness& witness) {
    std::cout << name << " value=" << witness.value << " k=" << witness.k
              << " minus=[" << witness.minus[0] << ',' << witness.minus[1]
              << "] plus=[";
    for (std::size_t i = 0; i < witness.plus.size(); ++i) {
        std::cout << (i == 0 ? "" : ",") << witness.plus[i];
    }
    const auto& s = witness.stats;
    std::cout << "] d=" << s.d << " c=" << s.active
              << " N=" << s.n << " P=" << s.p_equal
              << " U1=" << s.u_one << " U2=" << s.u_two
              << " T=" << s.t << '\n';
}

void merge_witnesses(
    std::array<Witness, 5>& raw,
    std::array<Witness, 5>& paid,
    std::array<Witness, 5>& u_needed,
    const std::array<Witness, 5>& local_raw,
    const std::array<Witness, 5>& local_paid,
    const std::array<Witness, 5>& local_u_needed
) {
    for (std::size_t d = 1; d <= 4; ++d) {
        if (local_raw[d].set) {
            consider(
                raw[d], local_raw[d].value, local_raw[d].k,
                local_raw[d].minus, local_raw[d].plus,
                local_raw[d].stats
            );
        }
        if (local_paid[d].set) {
            consider(
                paid[d], local_paid[d].value, local_paid[d].k,
                local_paid[d].minus, local_paid[d].plus,
                local_paid[d].stats
            );
        }
        if (local_u_needed[d].set) {
            consider(
                u_needed[d], local_u_needed[d].value,
                local_u_needed[d].k, local_u_needed[d].minus,
                local_u_needed[d].plus, local_u_needed[d].stats
            );
        }
    }
}

void record(
    std::array<Witness, 5>& raw,
    std::array<Witness, 5>& paid,
    std::array<Witness, 5>& u_needed,
    int k,
    const Labels2& minus,
    const Labels5& plus
) {
    if (!disjoint(minus, plus)) {
        return;
    }
    const auto s = statistics(k, minus, plus);
    if (s.d < 1 || s.d > 4) {
        return;
    }
    const auto raw_margin = s.n + s.p_equal - s.t;
    const auto paid_margin = raw_margin + s.u_one + s.u_two;
    consider(
        raw[static_cast<std::size_t>(s.d)],
        raw_margin, k, minus, plus, s
    );
    consider(
        paid[static_cast<std::size_t>(s.d)],
        paid_margin, k, minus, plus, s
    );
    if (raw_margin < 0) {
        consider(
            u_needed[static_cast<std::size_t>(s.d)],
            paid_margin, k, minus, plus, s
        );
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 9) {
        const int k = std::atoi(argv[1]);
        Labels2 minus{std::atoi(argv[2]), std::atoi(argv[3])};
        Labels5 plus{
            std::atoi(argv[4]), std::atoi(argv[5]), std::atoi(argv[6]),
            std::atoi(argv[7]), std::atoi(argv[8])
        };
        std::sort(minus.begin(), minus.end());
        std::sort(plus.begin(), plus.end());
        Witness witness{
            true, 0, k, minus, plus, statistics(k, minus, plus)
        };
        print("case", witness);
        return EXIT_SUCCESS;
    }
    if (argc == 4 && std::string(argv[1]) == "--exhaustive") {
        const int minimum_level = std::atoi(argv[2]);
        const int maximum_level = std::atoi(argv[3]);
        if (minimum_level < 2 || maximum_level < minimum_level
            || maximum_level >= 256) {
            std::cerr << "invalid exhaustive level range\n";
            return EXIT_FAILURE;
        }
        use_cache = true;
        std::array<Witness, 5> raw;
        std::array<Witness, 5> paid;
        std::array<Witness, 5> u_needed;
        for (int k = minimum_level; k <= maximum_level; ++k) {
#pragma omp parallel
            {
                std::array<Witness, 5> local_raw;
                std::array<Witness, 5> local_paid;
                std::array<Witness, 5> local_u_needed;
#pragma omp for schedule(dynamic)
                for (int q = 1; q <= k; ++q) {
                    for (int a = q; a <= k; ++a) {
                        const Labels2 minus{q, a};
                        for (int p1 = 1; p1 <= k; ++p1)
                        for (int p2 = p1; p2 <= k; ++p2)
                        for (int p3 = p2; p3 <= k; ++p3)
                        for (int p4 = p3; p4 <= k; ++p4)
                        for (int p5 = p4; p5 <= k; ++p5) {
                            const Labels5 plus{p1, p2, p3, p4, p5};
                            record(
                                local_raw, local_paid, local_u_needed,
                                k, minus, plus
                            );
                        }
                    }
                }
#pragma omp critical
                merge_witnesses(
                    raw, paid, u_needed,
                    local_raw, local_paid, local_u_needed
                );
            }
            std::cerr << "finished k=" << k << '\n';
        }
        for (std::size_t d = 1; d <= 4; ++d) {
            std::cout << "d=" << d << '\n';
            print("  raw", raw[d]);
            print("  paid", paid[d]);
            if (u_needed[d].set) {
                print("  u-needed", u_needed[d]);
            }
        }
        return EXIT_SUCCESS;
    }
    if (argc == 4 && std::string(argv[1]) == "--q1-four-odd") {
        const int minimum_level = std::atoi(argv[2]);
        const int maximum_level = std::atoi(argv[3]);
        if (minimum_level < 3 || maximum_level < minimum_level
            || maximum_level >= 256) {
            std::cerr << "invalid q1-four-odd level range\n";
            return EXIT_FAILURE;
        }
        use_cache = true;
        std::array<std::array<Witness, 3>, 5> raw;
        std::array<std::array<Witness, 3>, 5> paid;
        for (int k = minimum_level; k <= maximum_level; ++k) {
            const Labels2 minus{1, 1};
#pragma omp parallel
            {
            std::array<std::array<Witness, 3>, 5> local_raw;
            std::array<std::array<Witness, 3>, 5> local_paid;
#pragma omp for schedule(dynamic)
            for (int p1 = 2; p1 <= k; ++p1)
            for (int p2 = p1; p2 <= k; ++p2)
            for (int p3 = p2; p3 <= k; ++p3)
            for (int p4 = p3; p4 <= k; ++p4)
            for (int p5 = p4; p5 <= k; ++p5) {
                const Labels5 plus{p1, p2, p3, p4, p5};
                int odd = 0;
                for (const int p : plus) {
                    odd += p % 2;
                }
                if (odd != 4) {
                    continue;
                }
                const auto s = statistics(k, minus, plus);
                if (s.t == 0) {
                    continue;
                }
                const auto active = static_cast<std::size_t>(s.active / 2);
                const auto d = static_cast<std::size_t>(s.d);
                consider(
                    local_raw[active][d], s.n + s.p_equal - s.t,
                    k, minus, plus, s
                );
                consider(
                    local_paid[active][d],
                    s.n + s.p_equal + s.u_one + s.u_two - s.t,
                    k, minus, plus, s
                );
            }
#pragma omp critical
            {
                for (std::size_t active = 1; active <= 4; ++active) {
                    for (std::size_t d = 1; d <= 2; ++d) {
                        if (local_raw[active][d].set) {
                            consider(
                                raw[active][d], local_raw[active][d].value,
                                local_raw[active][d].k,
                                local_raw[active][d].minus,
                                local_raw[active][d].plus,
                                local_raw[active][d].stats
                            );
                        }
                        if (local_paid[active][d].set) {
                            consider(
                                paid[active][d], local_paid[active][d].value,
                                local_paid[active][d].k,
                                local_paid[active][d].minus,
                                local_paid[active][d].plus,
                                local_paid[active][d].stats
                            );
                        }
                    }
                }
            }
            }
        }
        for (std::size_t active = 1; active <= 4; ++active) {
            for (std::size_t d = 1; d <= 2; ++d) {
                if (!raw[active][d].set) {
                    continue;
                }
                std::cout << "active=" << active << " d=" << d << '\n';
                print("  raw", raw[active][d]);
                print("  paid", paid[active][d]);
            }
        }
        return EXIT_SUCCESS;
    }
    if (argc == 4 && (
        std::string(argv[1]) == "--q1-deficits"
        || std::string(argv[1]) == "--q1-tight"
    )) {
        const bool include_tight = std::string(argv[1]) == "--q1-tight";
        const char* q1_name =
            include_tight ? "q1-tight" : "q1-deficit";
        const int minimum_level = std::atoi(argv[2]);
        const int maximum_level = std::atoi(argv[3]);
        if (minimum_level < 2 || maximum_level < minimum_level
            || maximum_level >= 256) {
            std::cerr << "invalid q1 level range\n";
            return EXIT_FAILURE;
        }
        use_cache = true;
        std::uint64_t deficits = 0;
        Witness orientation_minimum;
        for (int k = minimum_level; k <= maximum_level; ++k) {
            const Labels2 minus{1, 1};
            for (int p1 = 2; p1 <= k; ++p1)
            for (int p2 = p1; p2 <= k; ++p2)
            for (int p3 = p2; p3 <= k; ++p3)
            for (int p4 = p3; p4 <= k; ++p4)
            for (int p5 = p4; p5 <= k; ++p5) {
                const Labels5 plus{p1, p2, p3, p4, p5};
                const auto s = statistics(k, minus, plus);
                consider(
                    orientation_minimum, s.n - s.t / 2,
                    k, minus, plus, s
                );
                if (s.t == 0) {
                    continue;
                }
                if (
                    include_tight
                        ? s.n + s.p_equal > s.t
                        : s.n + s.p_equal >= s.t
                ) {
                    continue;
                }
                ++deficits;
                Witness witness{
                    true, s.n + s.p_equal - s.t,
                    k, minus, plus, s
                };
                print(q1_name, witness);
            }
        }
        print("q1-orientation-minimum", orientation_minimum);
        std::cout << q1_name << "-count=" << deficits << '\n';
        return EXIT_SUCCESS;
    }
    if (argc != 3 && argc != 4) {
        std::cerr
            << "usage: analyze_su2_finite_small_cut MAX_LEVEL SAMPLES"
               " [MIN_LEVEL]\n"
            << "   or: analyze_su2_finite_small_cut --exhaustive"
               " MIN_LEVEL MAX_LEVEL\n"
            << "   or: analyze_su2_finite_small_cut"
               " --q1-deficits MIN_LEVEL MAX_LEVEL\n"
            << "   or: analyze_su2_finite_small_cut"
               " --q1-tight MIN_LEVEL MAX_LEVEL\n"
            << "   or: analyze_su2_finite_small_cut"
               " --q1-four-odd MIN_LEVEL MAX_LEVEL\n"
            << "   or: analyze_su2_finite_small_cut k q a p1 p2 p3 p4 p5\n";
        return EXIT_FAILURE;
    }
    const int maximum_level = std::atoi(argv[1]);
    const std::uint64_t samples = std::strtoull(argv[2], nullptr, 10);
    const int minimum_level = argc == 4 ? std::atoi(argv[3]) : 2;
    if (minimum_level < 2 || maximum_level < minimum_level
        || maximum_level >= 256 || samples == 0) {
        std::cerr << "invalid argument\n";
        return EXIT_FAILURE;
    }

    std::array<Witness, 5> raw;
    std::array<Witness, 5> paid;
    std::array<Witness, 5> u_needed;

#pragma omp parallel
    {
        std::array<Witness, 5> local_raw;
        std::array<Witness, 5> local_paid;
        std::array<Witness, 5> local_u_needed;
#pragma omp for schedule(static)
        for (std::uint64_t sample = 0; sample < samples; ++sample) {
            std::uint64_t state = sample + 0x9e3779b97f4a7c15ULL;
            auto random = [&state]() {
                state += 0x9e3779b97f4a7c15ULL;
                std::uint64_t z = state;
                z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
                z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
                return z ^ (z >> 31U);
            };
            const int k = minimum_level + static_cast<int>(
                random() % static_cast<std::uint64_t>(
                    maximum_level - minimum_level + 1
                )
            );
            Labels2 minus{};
            Labels5 plus{};
            for (int& value : minus) {
                value = 1 + static_cast<int>(
                    random() % static_cast<std::uint64_t>(k)
                );
            }
            for (int& value : plus) {
                value = 1 + static_cast<int>(
                    random() % static_cast<std::uint64_t>(k)
                );
            }
            std::sort(minus.begin(), minus.end());
            std::sort(plus.begin(), plus.end());
            if (!disjoint(minus, plus)) {
                continue;
            }
            const auto s = statistics(k, minus, plus);
            if (s.d < 1 || s.d > 4) {
                continue;
            }
            const auto raw_margin = s.n + s.p_equal - s.t;
            const auto paid_margin =
                raw_margin + s.u_one + s.u_two;
            consider(
                local_raw[static_cast<std::size_t>(s.d)],
                raw_margin, k, minus, plus, s
            );
            consider(
                local_paid[static_cast<std::size_t>(s.d)],
                paid_margin, k, minus, plus, s
            );
            if (raw_margin < 0) {
                consider(
                    local_u_needed[static_cast<std::size_t>(s.d)],
                    s.u_one + s.u_two + raw_margin,
                    k, minus, plus, s
                );
            }
        }
#pragma omp critical
        {
            for (std::size_t d = 1; d <= 4; ++d) {
                if (local_raw[d].set) {
                    consider(
                        raw[d], local_raw[d].value, local_raw[d].k,
                        local_raw[d].minus, local_raw[d].plus,
                        local_raw[d].stats
                    );
                }
                if (local_paid[d].set) {
                    consider(
                        paid[d], local_paid[d].value, local_paid[d].k,
                        local_paid[d].minus, local_paid[d].plus,
                        local_paid[d].stats
                    );
                }
                if (local_u_needed[d].set) {
                    consider(
                        u_needed[d], local_u_needed[d].value,
                        local_u_needed[d].k, local_u_needed[d].minus,
                        local_u_needed[d].plus, local_u_needed[d].stats
                    );
                }
            }
        }
    }

    for (std::size_t d = 1; d <= 4; ++d) {
        std::cout << "d=" << d << '\n';
        print("  raw", raw[d]);
        print("  paid", paid[d]);
        if (u_needed[d].set) {
            print("  u-needed", u_needed[d]);
        }
    }
}
