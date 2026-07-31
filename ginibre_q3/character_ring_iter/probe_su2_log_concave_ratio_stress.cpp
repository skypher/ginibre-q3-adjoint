#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

Integer at(const std::vector<Integer>& p, int i) {
    if (i < 0 || i >= static_cast<int>(p.size())) {
        return 0;
    }
    return p[static_cast<std::size_t>(i)];
}

std::vector<Integer> transform(const std::vector<Integer>& p, int a) {
    std::vector<Integer> result(p.size() + static_cast<std::size_t>(a));
    for (int i = 0; i < static_cast<int>(result.size()); ++i) {
        for (int h = std::abs(i - a); h <= i + a; ++h) {
            result[static_cast<std::size_t>(i)] += at(p, h);
        }
    }
    return result;
}

std::vector<Integer> suffix_products(
    const std::vector<Integer>& left,
    const std::vector<Integer>& right,
    int size
) {
    std::vector<Integer> suffix(static_cast<std::size_t>(size + 1));
    for (int i = size - 1; i >= 0; --i) {
        suffix[static_cast<std::size_t>(i)] =
            suffix[static_cast<std::size_t>(i + 1)] + at(left, i) * at(right, i);
    }
    return suffix;
}

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

unsigned int adaptive_threads() {
    const unsigned int hardware = std::max(1U, std::thread::hardware_concurrency());
    double load[1] = {0.0};
    const int read = getloadavg(load, 1);
    const unsigned int occupied =
        read == 1
        ? static_cast<unsigned int>(std::ceil(std::max(0.0, load[0])))
        : 0U;
    return std::max(
        1U,
        hardware > occupied + 2U ? hardware - occupied - 2U : 1U
    );
}

struct Failure {
    int sample = -1;
    int q = 0;
    int a = 0;
    int cutoff = 0;
    int denominator = 0;
    std::vector<int> numerators;
    Integer determinant = 0;
    std::vector<Integer> profile;
};

struct GapFailure {
    int sample = -1;
    int q = 0;
    int a = 0;
    int gap = 0;
    int denominator = 0;
    std::vector<int> numerators;
    Integer value = 0;
    Integer current = 0;
    std::vector<Integer> profile;
    std::vector<Integer> energies;
};

int replay_complete_wall_gap() {
    const std::vector<int> numerators{
        84, 82, 81, 72, 65, 64, 56, 45, 20, 15
    };
    const int denominator = 42;
    const std::vector<Integer> profile{
        Integer{"17080198121677824"},
        Integer{"34160396243355648"},
        Integer{"66694106951313408"},
        Integer{"128624349120390144"},
        Integer{"220498884206383104"},
        Integer{"341248273176545280"},
        Integer{"519997368649973760"},
        Integer{"693329824866631680"},
        Integer{"742853383785676800"},
        Integer{"353739706564608000"},
        Integer{"126335609487360000"}
    };
    for (std::size_t i = 1; i < profile.size(); ++i) {
        if (
            profile[i] * denominator
            != profile[i - 1] * numerators[i - 1]
        ) {
            std::cerr << "fixed witness ratio mismatch\n";
            return EXIT_FAILURE;
        }
    }

    const int q = 1;
    const int a = 7;
    const auto q_transform = transform(profile, q);
    const auto a_transform = transform(profile, a);
    const int size = std::max({
        static_cast<int>(profile.size()),
        static_cast<int>(q_transform.size()),
        static_cast<int>(a_transform.size())
    });
    const auto pp = suffix_products(profile, profile, size);
    const auto pq = suffix_products(profile, q_transform, size);
    const auto pa = suffix_products(profile, a_transform, size);
    const auto qa = suffix_products(q_transform, a_transform, size);
    const Integer current =
        pp.front() * qa.front() - pq.front() * pa.front();
    std::vector<Integer> energies(static_cast<std::size_t>(size));
    for (int gap = 1; gap < size; ++gap) {
        for (int left = 0; left + gap < size; ++left) {
            const int right = left + gap;
            const Integer q_wedge =
                at(profile, left) * at(q_transform, right)
                - at(profile, right) * at(q_transform, left);
            const Integer a_wedge =
                at(profile, left) * at(a_transform, right)
                - at(profile, right) * at(a_transform, left);
            energies[static_cast<std::size_t>(gap)] += q_wedge * a_wedge;
        }
    }
    const Integer expected_gap{
        "-204839305586958588388030857811270968817700301343897173126550860595200"
    };
    const Integer expected_current{
        "2943060361665841378980816615817878511514771508655754812817654026598350848"
    };
    Integer energy_sum = 0;
    for (const Integer& energy : energies) {
        energy_sum += energy;
    }
    const bool pass =
        energies[9] == expected_gap
        && current == expected_current
        && energy_sum == current;

    std::cout
        << "SU2_LOG_CONCAVE_RATIO_STRESS replay_complete_wall_gap"
        << " q=" << q
        << " a=" << a
        << " gap=9"
        << " denominator=" << denominator
        << " ratio_numerators={";
    for (std::size_t i = 0; i < numerators.size(); ++i) {
        if (i != 0) {
            std::cout << ',';
        }
        std::cout << numerators[i];
    }
    std::cout
        << "} value=" << energies[9]
        << " current=" << current
        << " energy_sum=" << energy_sum
        << " gap_energies={";
    for (std::size_t gap = 1; gap < energies.size(); ++gap) {
        if (gap != 1U) {
            std::cout << ',';
        }
        std::cout << energies[gap];
    }
    std::cout
        << "} result=" << (pass ? "PASS" : "FAIL")
        << '\n';
    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}

int replay_complete_wall_gap_tail() {
    const std::vector<int> numerators{
        21, 20, 20, 18, 17, 16, 14, 13, 4, 3, 2
    };
    const int denominator = 11;
    const std::vector<Integer> profile{
        Integer{"285311670611"},
        Integer{"544685916621"},
        Integer{"990338030220"},
        Integer{"1800614600400"},
        Integer{"2946460255200"},
        Integer{"4553620394400"},
        Integer{"6623447846400"},
        Integer{"8429842713600"},
        Integer{"9962541388800"},
        Integer{"3622742323200"},
        Integer{"988020633600"},
        Integer{"179640115200"}
    };
    for (std::size_t i = 1; i < profile.size(); ++i) {
        if (
            profile[i] * denominator
            != profile[i - 1] * numerators[i - 1]
        ) {
            throw std::runtime_error(
                "complete-wall gap-tail witness ratio mismatch"
            );
        }
    }

    constexpr int q = 1;
    constexpr int a = 7;
    constexpr int minimum_gap = 8;
    const auto q_transform = transform(profile, q);
    const auto a_transform = transform(profile, a);
    const int size = std::max({
        static_cast<int>(profile.size()),
        static_cast<int>(q_transform.size()),
        static_cast<int>(a_transform.size())
    });
    const auto pp = suffix_products(profile, profile, size);
    const auto pq = suffix_products(profile, q_transform, size);
    const auto pa = suffix_products(profile, a_transform, size);
    const auto qa = suffix_products(q_transform, a_transform, size);
    const Integer current =
        pp.front() * qa.front() - pq.front() * pa.front();
    std::vector<Integer> energies(static_cast<std::size_t>(size));
    for (int gap = 1; gap < size; ++gap) {
        for (int left = 0; left + gap < size; ++left) {
            const int right = left + gap;
            const Integer q_wedge =
                at(profile, left) * at(q_transform, right)
                - at(profile, right) * at(q_transform, left);
            const Integer a_wedge =
                at(profile, left) * at(a_transform, right)
                - at(profile, right) * at(a_transform, left);
            energies[static_cast<std::size_t>(gap)] +=
                q_wedge * a_wedge;
        }
    }
    std::vector<Integer> tails(energies.size());
    Integer running = 0;
    for (int gap = size - 1; gap >= 1; --gap) {
        running += energies[static_cast<std::size_t>(gap)];
        tails[static_cast<std::size_t>(gap)] = running;
    }
    const Integer expected_tail{
        "-16682814543009861315632046030279412129824776371200"
    };
    const Integer expected_current{
        "94153010945848127878367241684587477315345747381424000"
    };
    if (
        tails[minimum_gap] != expected_tail
        || current != expected_current
        || tails[1] != current
    ) {
        throw std::runtime_error(
            "complete-wall gap-tail obstruction replay mismatch"
        );
    }
    std::cout
        << "SU2_COMPLETE_WALL_GAP_TAIL_OBSTRUCTION"
        << " q=" << q
        << " a=" << a
        << " minimum_gap=" << minimum_gap
        << " value=" << tails[minimum_gap]
        << " current=" << current
        << " result=PASS_EXACT\n";
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    if (
        argc == 2
        && std::string{argv[1]} == "--replay-complete-wall-gap"
    ) {
        return replay_complete_wall_gap();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-complete-wall-gap-tail"
    ) {
        return replay_complete_wall_gap_tail();
    }
    int samples = 20000;
    int max_label = 7;
    int maximum_length = 14;
    int maximum_ratio = 2;
    try {
        if (argc >= 2) {
            samples = parse_positive(argv[1], "samples");
        }
        if (argc >= 3) {
            max_label = parse_positive(argv[2], "max_label");
        }
        if (argc >= 4) {
            maximum_length = parse_positive(argv[3], "maximum_length");
            if (maximum_length < 2) {
                throw std::invalid_argument(
                    "maximum_length must be at least two"
                );
            }
        }
        if (argc >= 5) {
            maximum_ratio = parse_positive(argv[4], "maximum_ratio");
        }
        if (argc > 5) {
            throw std::invalid_argument(
                "usage: probe_su2_log_concave_ratio_stress "
                "[samples] [max_label] [maximum_length] [maximum_ratio] "
                "| --replay-complete-wall-gap "
                "| --replay-complete-wall-gap-tail"
            );
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }

    const unsigned int threads = adaptive_threads();
    std::atomic<int> next_sample{0};
    std::atomic<unsigned long long> determinants{0};
    std::atomic<unsigned long long> counterexamples{0};
    std::atomic<unsigned long long> complete_wall_gap_energies{0};
    std::atomic<unsigned long long> complete_wall_gap_failures{0};
    std::atomic<unsigned long long> complete_wall_gap_tails{0};
    std::atomic<unsigned long long> complete_wall_gap_tail_failures{0};
    std::atomic<unsigned long long> complete_wall_gap_prefixes{0};
    std::atomic<unsigned long long> complete_wall_gap_prefix_failures{0};
    std::mutex failure_mutex;
    Failure failure;
    GapFailure gap_failure;

    auto worker = [&]() {
        while (true) {
            const int sample =
                next_sample.fetch_add(1, std::memory_order_relaxed);
            if (sample >= samples) {
                return;
            }
            std::mt19937_64 generator(
                0x9e3779b97f4a7c15ULL
                ^ (
                    static_cast<unsigned long long>(sample)
                    * 0xbf58476d1ce4e5b9ULL
                )
            );
            const int length =
                2
                + static_cast<int>(
                    generator()
                    % static_cast<unsigned long long>(maximum_length - 1)
                );
            const int denominator = 2 + static_cast<int>(generator() % 63ULL);
            const int shift = static_cast<int>(generator() % 5ULL);
            std::vector<int> numerators(
                static_cast<std::size_t>(length - 1)
            );
            for (int& numerator : numerators) {
                numerator =
                    1
                    + static_cast<int>(
                        generator()
                        % static_cast<unsigned long long>(
                            maximum_ratio * denominator
                        )
                    );
            }
            std::sort(
                numerators.begin(),
                numerators.end(),
                std::greater<int>()
            );

            std::vector<Integer> core(static_cast<std::size_t>(length));
            core[0] = 1;
            for (int power = 1; power < length; ++power) {
                core[0] *= denominator;
            }
            for (int i = 1; i < length; ++i) {
                core[static_cast<std::size_t>(i)] =
                    core[static_cast<std::size_t>(i - 1)]
                    * numerators[static_cast<std::size_t>(i - 1)]
                    / denominator;
            }
            std::vector<Integer> profile(
                static_cast<std::size_t>(shift),
                Integer{0}
            );
            profile.insert(profile.end(), core.begin(), core.end());

            std::vector<std::vector<Integer>> transforms(
                static_cast<std::size_t>(max_label + 1)
            );
            for (int label = 1; label <= max_label; ++label) {
                transforms[static_cast<std::size_t>(label)] =
                    transform(profile, label);
            }
            for (int q = 1; q <= max_label; ++q) {
                for (int a = 1; a <= max_label; ++a) {
                    const int size = std::max({
                        static_cast<int>(profile.size()),
                        static_cast<int>(
                            transforms[static_cast<std::size_t>(q)].size()
                        ),
                        static_cast<int>(
                            transforms[static_cast<std::size_t>(a)].size()
                        )
                    });
                    const auto pp = suffix_products(profile, profile, size);
                    const auto pq = suffix_products(
                        profile,
                        transforms[static_cast<std::size_t>(q)],
                        size
                    );
                    const auto pa = suffix_products(
                        profile,
                        transforms[static_cast<std::size_t>(a)],
                        size
                    );
                    const auto qa = suffix_products(
                        transforms[static_cast<std::size_t>(q)],
                        transforms[static_cast<std::size_t>(a)],
                        size
                    );
                    std::vector<Integer> gap_energy(
                        static_cast<std::size_t>(size)
                    );
                    int first_negative_gap = -1;
                    Integer first_negative_gap_value = 0;
                    for (int gap = 1; gap < size; ++gap) {
                        Integer energy = 0;
                        for (int left = 0; left + gap < size; ++left) {
                            const int right = left + gap;
                            const Integer q_wedge =
                                at(profile, left)
                                    * at(
                                        transforms[
                                            static_cast<std::size_t>(q)
                                        ],
                                        right
                                    )
                                - at(profile, right)
                                    * at(
                                        transforms[
                                            static_cast<std::size_t>(q)
                                        ],
                                        left
                                    );
                            const Integer a_wedge =
                                at(profile, left)
                                    * at(
                                        transforms[
                                            static_cast<std::size_t>(a)
                                        ],
                                        right
                                    )
                                - at(profile, right)
                                    * at(
                                        transforms[
                                            static_cast<std::size_t>(a)
                                        ],
                                        left
                                    );
                            energy += q_wedge * a_wedge;
                        }
                        gap_energy[static_cast<std::size_t>(gap)] = energy;
                        complete_wall_gap_energies.fetch_add(
                            1,
                            std::memory_order_relaxed
                        );
                        if (energy < 0) {
                            complete_wall_gap_failures.fetch_add(
                                1,
                                std::memory_order_relaxed
                            );
                            if (first_negative_gap < 0) {
                                first_negative_gap = gap;
                                first_negative_gap_value = energy;
                            }
                        }
                    }
                    if (first_negative_gap >= 0) {
                        const Integer complete_current =
                            pp.front() * qa.front() -
                            pq.front() * pa.front();
                        std::lock_guard<std::mutex> lock(failure_mutex);
                        if (
                            gap_failure.sample < 0
                            || sample < gap_failure.sample
                        ) {
                            gap_failure = {
                                sample,
                                q,
                                a,
                                first_negative_gap,
                                denominator,
                                numerators,
                                first_negative_gap_value,
                                complete_current,
                                profile,
                                gap_energy
                            };
                        }
                    }
                    Integer complete_gap_tail = 0;
                    for (int gap = size - 1; gap >= 1; --gap) {
                        complete_gap_tail +=
                            gap_energy[static_cast<std::size_t>(gap)];
                        complete_wall_gap_tails.fetch_add(
                            1,
                            std::memory_order_relaxed
                        );
                        if (complete_gap_tail < 0) {
                            complete_wall_gap_tail_failures.fetch_add(
                                1,
                                std::memory_order_relaxed
                            );
                        }
                    }
                    Integer complete_gap_prefix = 0;
                    for (int gap = 1; gap < size; ++gap) {
                        complete_gap_prefix +=
                            gap_energy[static_cast<std::size_t>(gap)];
                        complete_wall_gap_prefixes.fetch_add(
                            1,
                            std::memory_order_relaxed
                        );
                        if (complete_gap_prefix < 0) {
                            complete_wall_gap_prefix_failures.fetch_add(
                                1,
                                std::memory_order_relaxed
                            );
                        }
                    }
                    for (int cutoff = 0; cutoff <= size; ++cutoff) {
                        const Integer determinant =
                            pp[static_cast<std::size_t>(cutoff)]
                                * qa[static_cast<std::size_t>(cutoff)]
                            - pq[static_cast<std::size_t>(cutoff)]
                                * pa[static_cast<std::size_t>(cutoff)];
                        determinants.fetch_add(1, std::memory_order_relaxed);
                        if (determinant < 0) {
                            counterexamples.fetch_add(
                                1,
                                std::memory_order_relaxed
                            );
                            std::lock_guard<std::mutex> lock(failure_mutex);
                            if (
                                failure.sample < 0
                                || sample < failure.sample
                            ) {
                                failure = {
                                    sample,
                                    q,
                                    a,
                                    cutoff,
                                    denominator,
                                    numerators,
                                    determinant,
                                    profile
                                };
                            }
                        }
                    }
                }
            }
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned int thread = 0; thread < threads; ++thread) {
        workers.emplace_back(worker);
    }
    for (std::thread& thread : workers) {
        thread.join();
    }

    if (counterexamples.load() != 0) {
        std::cout
            << "SU2_LOG_CONCAVE_RATIO_STRESS counterexample"
            << " sample=" << failure.sample
            << " q=" << failure.q
            << " a=" << failure.a
            << " cutoff=" << failure.cutoff
            << " denominator=" << failure.denominator
            << " ratio_numerators={";
        for (std::size_t i = 0; i < failure.numerators.size(); ++i) {
            if (i != 0) {
                std::cout << ',';
            }
            std::cout << failure.numerators[i];
        }
        std::cout
            << "}"
            << " determinant=" << failure.determinant
            << " profile={";
        for (std::size_t i = 0; i < failure.profile.size(); ++i) {
            if (i != 0) {
                std::cout << ',';
            }
            std::cout << failure.profile[i];
        }
        std::cout
            << "} determinants=" << determinants.load()
            << " counterexamples=" << counterexamples.load()
            << " complete_wall_gap_energies="
            << complete_wall_gap_energies.load()
            << " complete_wall_gap_failures="
            << complete_wall_gap_failures.load()
            << " complete_wall_gap_tails="
            << complete_wall_gap_tails.load()
            << " complete_wall_gap_tail_failures="
            << complete_wall_gap_tail_failures.load()
            << " complete_wall_gap_prefixes="
            << complete_wall_gap_prefixes.load()
            << " complete_wall_gap_prefix_failures="
            << complete_wall_gap_prefix_failures.load()
            << " threads=" << threads
            << " maximum_length=" << maximum_length
            << " maximum_ratio=" << maximum_ratio
            << " result=COUNTEREXAMPLE\n";
        return EXIT_SUCCESS;
    }
    std::cout
        << "SU2_LOG_CONCAVE_RATIO_STRESS"
        << " samples=" << samples
        << " determinants=" << determinants.load()
        << " max_label=" << max_label
        << " complete_wall_gap_energies="
        << complete_wall_gap_energies.load()
        << " complete_wall_gap_failures="
        << complete_wall_gap_failures.load()
        << " complete_wall_gap_tails="
        << complete_wall_gap_tails.load()
        << " complete_wall_gap_tail_failures="
        << complete_wall_gap_tail_failures.load()
        << " complete_wall_gap_prefixes="
        << complete_wall_gap_prefixes.load()
        << " complete_wall_gap_prefix_failures="
        << complete_wall_gap_prefix_failures.load()
        << " maximum_length=" << maximum_length
        << " maximum_ratio=" << maximum_ratio;
    if (gap_failure.sample >= 0) {
        std::cout
            << " first_complete_wall_gap={sample=" << gap_failure.sample
            << ",q=" << gap_failure.q
            << ",a=" << gap_failure.a
            << ",gap=" << gap_failure.gap
            << ",denominator=" << gap_failure.denominator
            << ",ratio_numerators={";
        for (std::size_t i = 0; i < gap_failure.numerators.size(); ++i) {
            if (i != 0) {
                std::cout << ',';
            }
            std::cout << gap_failure.numerators[i];
        }
        std::cout
            << "},value=" << gap_failure.value
            << ",current=" << gap_failure.current
            << ",profile={";
        for (std::size_t i = 0; i < gap_failure.profile.size(); ++i) {
            if (i != 0) {
                std::cout << ',';
            }
            std::cout << gap_failure.profile[i];
        }
        std::cout << "},gap_energies={";
        for (std::size_t gap = 1; gap < gap_failure.energies.size(); ++gap) {
            if (gap != 1U) {
                std::cout << ',';
            }
            std::cout << gap_failure.energies[gap];
        }
        std::cout << "}}";
    }
    std::cout
        << " threads=" << threads
        << " result="
        << (
            gap_failure.sample < 0
            ? "NO_COUNTEREXAMPLE"
            : "COMPLETE_WALL_GAP_COUNTEREXAMPLE"
        )
        << '\n';
    return EXIT_SUCCESS;
}
