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

}  // namespace

int main(int argc, char** argv) {
    int samples = 20000;
    int max_label = 7;
    try {
        if (argc >= 2) {
            samples = parse_positive(argv[1], "samples");
        }
        if (argc >= 3) {
            max_label = parse_positive(argv[2], "max_label");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_log_concave_ratio_stress [samples] [max_label]"
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
    std::mutex failure_mutex;
    Failure failure;

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
            const int length = 2 + static_cast<int>(generator() % 13ULL);
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
                        % static_cast<unsigned long long>(2 * denominator)
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
            << " threads=" << threads
            << " result=COUNTEREXAMPLE\n";
        return EXIT_SUCCESS;
    }
    std::cout
        << "SU2_LOG_CONCAVE_RATIO_STRESS"
        << " samples=" << samples
        << " determinants=" << determinants.load()
        << " max_label=" << max_label
        << " threads=" << threads
        << " result=NO_COUNTEREXAMPLE\n";
    return EXIT_SUCCESS;
}
