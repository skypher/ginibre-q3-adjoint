#include <boost/multiprecision/cpp_int.hpp>

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

namespace {

using Integer = boost::multiprecision::cpp_int;

Integer at(const std::vector<Integer>& profile, int index) {
    return
        index >= 0 && index < static_cast<int>(profile.size())
        ? profile[static_cast<std::size_t>(index)]
        : Integer{0};
}

std::vector<Integer> transform(
    const std::vector<Integer>& profile,
    int label
) {
    std::vector<Integer> result(
        profile.size() + static_cast<std::size_t>(label)
    );
    for (int target = 0; target < static_cast<int>(result.size()); ++target) {
        for (int source = std::abs(target - label);
             source <= target + label;
             ++source) {
            result[static_cast<std::size_t>(target)] += at(profile, source);
        }
    }
    return result;
}

Integer inner(
    const std::vector<Integer>& left,
    const std::vector<Integer>& right
) {
    Integer result = 0;
    const int size = std::max(left.size(), right.size());
    for (int index = 0; index < size; ++index) {
        result += at(left, index) * at(right, index);
    }
    return result;
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
    const unsigned int hardware =
        std::max(1U, std::thread::hardware_concurrency());
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
    Integer determinant = 0;
    std::vector<Integer> profile;
};

std::string show(const std::vector<Integer>& profile) {
    std::string result = "{";
    for (std::size_t index = 0; index < profile.size(); ++index) {
        if (index != 0U) {
            result.push_back(',');
        }
        result += profile[index].convert_to<std::string>();
    }
    result.push_back('}');
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int samples = 100000;
        int maximum_label = 10;
        if (argc >= 2) {
            samples = parse_positive(argv[1], "samples");
        }
        if (argc >= 3) {
            maximum_label =
                parse_positive(argv[2], "maximum_label");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_log_concave_full_covariance"
                " [samples] [maximum_label]"
            );
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
                    0xd1b54a32d192ed03ULL
                    ^ (
                        static_cast<unsigned long long>(sample)
                        * 0x9e3779b97f4a7c15ULL
                    )
                );
                const int length =
                    2 + static_cast<int>(generator() % 19ULL);
                const int denominator =
                    2 + static_cast<int>(generator() % 127ULL);
                std::vector<int> numerators(
                    static_cast<std::size_t>(length - 1)
                );
                for (int& numerator : numerators) {
                    numerator =
                        1
                        + static_cast<int>(
                            generator()
                            % static_cast<unsigned long long>(
                                3 * denominator
                            )
                        );
                }
                std::sort(
                    numerators.begin(),
                    numerators.end(),
                    std::greater<int>()
                );
                std::vector<Integer> profile(
                    static_cast<std::size_t>(length)
                );
                profile[0] = 1;
                for (int power = 1; power < length; ++power) {
                    profile[0] *= denominator;
                }
                for (int index = 1; index < length; ++index) {
                    profile[static_cast<std::size_t>(index)] =
                        profile[static_cast<std::size_t>(index - 1)]
                        * numerators[static_cast<std::size_t>(index - 1)]
                        / denominator;
                }

                std::vector<std::vector<Integer>> transforms(
                    static_cast<std::size_t>(maximum_label + 1)
                );
                for (int label = 1; label <= maximum_label; ++label) {
                    transforms[static_cast<std::size_t>(label)] =
                        transform(profile, label);
                }
                const Integer norm = inner(profile, profile);
                std::vector<Integer> anchored(
                    static_cast<std::size_t>(maximum_label + 1)
                );
                for (int label = 1; label <= maximum_label; ++label) {
                    anchored[static_cast<std::size_t>(label)] =
                        inner(
                            profile,
                            transforms[static_cast<std::size_t>(label)]
                        );
                }
                for (int q = 1; q <= maximum_label; ++q) {
                    for (int a = 1; a <= maximum_label; ++a) {
                        const Integer determinant =
                            norm
                            * inner(
                                transforms[static_cast<std::size_t>(q)],
                                transforms[static_cast<std::size_t>(a)]
                            )
                            - anchored[static_cast<std::size_t>(q)]
                                * anchored[static_cast<std::size_t>(a)];
                        determinants.fetch_add(
                            1,
                            std::memory_order_relaxed
                        );
                        if (determinant < 0) {
                            counterexamples.fetch_add(
                                1,
                                std::memory_order_relaxed
                            );
                            std::lock_guard<std::mutex> lock(
                                failure_mutex
                            );
                            if (
                                failure.sample < 0
                                || sample < failure.sample
                            ) {
                                failure = {
                                    sample,
                                    q,
                                    a,
                                    determinant,
                                    profile
                                };
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

        std::cout
            << "SU2_LOG_CONCAVE_FULL_COVARIANCE"
            << " samples=" << samples
            << " maximum_label=" << maximum_label
            << " determinants=" << determinants.load()
            << " counterexamples=" << counterexamples.load()
            << " threads=" << threads;
        if (failure.sample >= 0) {
            std::cout
                << " first={sample=" << failure.sample
                << ",q=" << failure.q
                << ",a=" << failure.a
                << ",determinant=" << failure.determinant
                << ",profile=" << show(failure.profile)
                << "}";
        }
        std::cout
            << " result="
            << (
                failure.sample < 0
                ? "NO_COUNTEREXAMPLE"
                : "COUNTEREXAMPLE"
            )
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
