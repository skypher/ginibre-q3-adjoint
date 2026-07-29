#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

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

Integer at(const std::vector<Integer>& sequence, int index) {
    return
        0 <= index && index < static_cast<int>(sequence.size())
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

long double character_density(
    const std::vector<long double>& sequence,
    long double x
) {
    long double result = sequence[0];
    if (sequence.size() == 1U) {
        return result;
    }
    long double previous = 1.0L;
    long double current = 2.0L * x;
    int degree = 1;
    for (std::size_t label = 1U; label < sequence.size(); ++label) {
        long double even_character = 0.0L;
        while (degree < 2 * static_cast<int>(label)) {
            const long double next =
                2.0L * x * current - previous;
            previous = current;
            current = next;
            ++degree;
        }
        even_character = current;
        result += sequence[label] * even_character;
    }
    return result;
}

long double sampled_density_minimum(
    const std::vector<Integer>& sequence,
    int grid
) {
    Integer scale = 0;
    for (const Integer& value : sequence) {
        scale = std::max(scale, value);
    }
    std::vector<long double> normalized(sequence.size());
    const long double denominator = scale.convert_to<long double>();
    for (std::size_t index = 0; index < sequence.size(); ++index) {
        normalized[index] =
            sequence[index].convert_to<long double>() / denominator;
    }

    long double minimum = std::numeric_limits<long double>::infinity();
    for (int point = 0; point <= grid; ++point) {
        const long double x =
            static_cast<long double>(point)
            / static_cast<long double>(grid);
        minimum = std::min(
            minimum,
            character_density(normalized, x)
        );
    }
    return minimum;
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
    Integer determinant = 0;
    long double density_minimum = 0.0L;
    std::vector<Integer> sequence;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        int samples = 200000;
        int maximum_label = 14;
        int grid = 4096;
        if (argc >= 2) {
            samples = parse_positive(argv[1], "samples");
        }
        if (argc >= 3) {
            maximum_label = parse_positive(argv[2], "maximum_label");
        }
        if (argc >= 4) {
            grid = parse_positive(argv[3], "grid");
        }
        if (argc > 4) {
            throw std::invalid_argument(
                "usage: probe_su2_positive_density_covariance"
                " [samples] [maximum_label] [grid]"
            );
        }

        const unsigned int threads = adaptive_threads();
        std::atomic<int> next_sample{0};
        std::atomic<unsigned long long> negative_profiles{0};
        std::atomic<unsigned long long>
            endpoint_nonnegative_negative_profiles{0};
        std::atomic<unsigned long long> sampled_nonnegative_profiles{0};
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
                    0x243f6a8885a308d3ULL
                    ^ (
                        static_cast<unsigned long long>(sample)
                        * 0x9e3779b97f4a7c15ULL
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

                Integer first_negative = 0;
                int first_q = 0;
                int first_a = 0;
                for (int q = 1; q <= maximum_label; ++q) {
                    for (int a = 1; a <= maximum_label; ++a) {
                        const Integer value =
                            determinant(sequence, q, a);
                        if (value < 0) {
                            first_negative = value;
                            first_q = q;
                            first_a = a;
                            break;
                        }
                    }
                    if (first_negative < 0) {
                        break;
                    }
                }
                if (first_negative >= 0) {
                    continue;
                }
                negative_profiles.fetch_add(1, std::memory_order_relaxed);

                std::vector<Integer> endpoint_sequence = sequence;
                Integer alternating_value = 0;
                for (std::size_t index = 0;
                     index < endpoint_sequence.size();
                     ++index) {
                    alternating_value +=
                        index % 2U == 0U
                        ? endpoint_sequence[index]
                        : -endpoint_sequence[index];
                }
                if (alternating_value >= 0) {
                    endpoint_nonnegative_negative_profiles.fetch_add(
                        1,
                        std::memory_order_relaxed
                    );
                }

                const long double minimum =
                    sampled_density_minimum(sequence, grid);
                if (minimum < -1.0e-12L) {
                    continue;
                }
                sampled_nonnegative_profiles.fetch_add(
                    1,
                    std::memory_order_relaxed
                );
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (failure.sample < 0 || sample < failure.sample) {
                    failure = {
                        sample,
                        first_q,
                        first_a,
                        first_negative,
                        minimum,
                        sequence
                    };
                }
            }
        };

        std::vector<std::thread> workers;
        workers.reserve(threads);
        for (unsigned int thread = 0; thread < threads; ++thread) {
            workers.emplace_back(worker);
        }
        for (std::thread& worker_thread : workers) {
            worker_thread.join();
        }

        std::cout
            << "SU2_POSITIVE_DENSITY_COVARIANCE"
            << " samples=" << samples
            << " maximum_label=" << maximum_label
            << " grid=" << grid
            << " negative_profiles=" << negative_profiles.load()
            << " endpoint_nonnegative_negative_profiles="
            << endpoint_nonnegative_negative_profiles.load()
            << " sampled_nonnegative_profiles="
            << sampled_nonnegative_profiles.load()
            << " threads=" << threads;
        if (failure.sample >= 0) {
            std::cout
                << " first={sample=" << failure.sample
                << ",q=" << failure.q
                << ",a=" << failure.a
                << ",determinant=" << failure.determinant
                << ",sampled_density_minimum="
                << static_cast<double>(failure.density_minimum)
                << ",sequence=" << show(failure.sequence)
                << "} result=CANDIDATE_COUNTEREXAMPLE\n";
        } else {
            std::cout << " result=NO_CANDIDATE\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
