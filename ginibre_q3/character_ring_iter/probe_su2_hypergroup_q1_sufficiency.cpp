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

Integer at(const std::vector<Integer>& sequence, int index) {
    return
        index >= 0 && index < static_cast<int>(sequence.size())
        ? sequence[static_cast<std::size_t>(index)]
        : Integer{0};
}

Integer determinant(
    const std::vector<Integer>& sequence,
    int q,
    int a
) {
    Integer translated = 0;
    for (int label = std::abs(q - a); label <= q + a; ++label) {
        translated += at(sequence, label);
    }
    return
        sequence[0] * translated
        - at(sequence, q) * at(sequence, a);
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

std::string show(const std::vector<Integer>& sequence) {
    std::string result = "{";
    for (std::size_t index = 0; index < sequence.size(); ++index) {
        if (index != 0U) {
            result.push_back(',');
        }
        result += sequence[index].convert_to<std::string>();
    }
    result.push_back('}');
    return result;
}

struct Failure {
    int sample = -1;
    int q = 0;
    int a = 0;
    Integer value = 0;
    std::vector<Integer> sequence;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        int samples = 100000;
        int maximum_label = 12;
        if (argc >= 2) {
            samples = parse_positive(argv[1], "samples");
        }
        if (argc >= 3) {
            maximum_label = parse_positive(argv[2], "maximum_label");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_hypergroup_q1_sufficiency"
                " [samples] [maximum_label]"
            );
        }

        const unsigned int threads = adaptive_threads();
        std::atomic<int> next_sample{0};
        std::atomic<unsigned long long> diagonal_profiles{0};
        std::atomic<unsigned long long> tested_determinants{0};
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
                    0xa0761d6478bd642fULL
                    ^ (
                        static_cast<unsigned long long>(sample)
                        * 0xe7037ed1a0b428dbULL
                    )
                );
                const int length =
                    2 + static_cast<int>(generator() % 19ULL);
                constexpr int denominator = 1024;
                std::vector<int> numerators(
                    static_cast<std::size_t>(length - 1)
                );
                for (int& numerator : numerators) {
                    const int exponent =
                        static_cast<int>(generator() % 21ULL);
                    const int mantissa =
                        1 + static_cast<int>(generator() % 3ULL);
                    numerator = mantissa * (1 << exponent);
                }
                std::sort(
                    numerators.begin(),
                    numerators.end(),
                    std::greater<int>()
                );
                std::vector<Integer> sequence(
                    static_cast<std::size_t>(length)
                );
                sequence[0] = 1;
                for (int index = 1; index < length; ++index) {
                    sequence[0] *= denominator;
                }
                for (int index = 1; index < length; ++index) {
                    sequence[static_cast<std::size_t>(index)] =
                        sequence[static_cast<std::size_t>(index - 1)]
                        * numerators[static_cast<std::size_t>(index - 1)]
                        / denominator;
                }

                bool diagonals_pass = true;
                for (int q = 1; q <= maximum_label; ++q) {
                    if (determinant(sequence, q, q) < 0) {
                        diagonals_pass = false;
                        break;
                    }
                }
                if (!diagonals_pass) {
                    continue;
                }
                diagonal_profiles.fetch_add(1, std::memory_order_relaxed);
                for (int q = 1; q <= maximum_label; ++q) {
                    for (int a = 1; a <= maximum_label; ++a) {
                        if (q == a) {
                            continue;
                        }
                        const Integer value =
                            determinant(sequence, q, a);
                        tested_determinants.fetch_add(
                            1,
                            std::memory_order_relaxed
                        );
                        if (value < 0) {
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
                                    value,
                                    sequence
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
            << "SU2_HYPERGROUP_Q1_SUFFICIENCY"
            << " samples=" << samples
            << " diagonal_profiles=" << diagonal_profiles.load()
            << " tested_determinants=" << tested_determinants.load()
            << " threads=" << threads;
        if (failure.sample >= 0) {
            std::cout
                << " first={sample=" << failure.sample
                << ",q=" << failure.q
                << ",a=" << failure.a
                << ",value=" << failure.value
                << ",sequence=" << show(failure.sequence)
                << "} result=COUNTEREXAMPLE_TO_REDUCTION\n";
        } else {
            std::cout << " result=NO_COUNTEREXAMPLE\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
