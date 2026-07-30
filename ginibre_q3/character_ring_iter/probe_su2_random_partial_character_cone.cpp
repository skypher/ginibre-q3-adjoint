#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

namespace {

using boost::multiprecision::cpp_int;

class SplitMix64 {
public:
    explicit SplitMix64(std::uint64_t state) : state_(state) {}

    std::uint64_t operator()() {
        std::uint64_t value =
            (state_ += UINT64_C(0x9e3779b97f4a7c15));
        value = (value ^ (value >> 30U))
            * UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27U))
            * UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31U);
    }

    int uniform_int(int lower, int upper) {
        if (lower > upper) {
            throw std::runtime_error("invalid random interval");
        }
        const std::uint64_t width =
            static_cast<std::uint64_t>(upper - lower + 1);
        return lower + static_cast<int>((*this)() % width);
    }

private:
    std::uint64_t state_;
};

template <class Function>
void for_each_output(int level, int left, int right, Function function) {
    const int lower = std::abs(left - right);
    const int upper = std::min(
        left + right,
        2 * level - left - right
    );
    for (int output = lower; output <= upper; output += 2) {
        function(output);
    }
}

struct Failure {
    std::uint64_t sample = 0U;
    int level = 0;
    int target = 0;
    cpp_int value = 0;
    std::vector<int> signed_labels;
};

bool inspect(
    int level,
    const std::vector<int>& signed_labels,
    Failure& failure
) {
    const int width = level + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(width * width),
        0
    );
    current[0] = 1;
    for (const int signed_label : signed_labels) {
        const int label = std::abs(signed_label);
        const bool minus = signed_label < 0;
        std::vector<cpp_int> next(current.size(), 0);
        for (int left = 0; left <= level; ++left) {
            for (int right = 0; right <= level; ++right) {
                const cpp_int& coefficient = current[
                    static_cast<std::size_t>(left * width + right)
                ];
                if (coefficient == 0) {
                    continue;
                }
                for_each_output(
                    level,
                    left,
                    label,
                    [&](int output) {
                        next[static_cast<std::size_t>(
                            output * width + right
                        )] += coefficient;
                    }
                );
                for_each_output(
                    level,
                    right,
                    label,
                    [&](int output) {
                        cpp_int& destination =
                            next[static_cast<std::size_t>(
                                left * width + output
                            )];
                        destination += minus
                            ? -coefficient
                            : coefficient;
                    }
                );
            }
        }
        current = std::move(next);
    }
    for (int target = 0; target <= level; ++target) {
        const cpp_int& value = current[
            static_cast<std::size_t>(target * width)
        ];
        if (value < 0) {
            failure.level = level;
            failure.target = target;
            failure.value = value;
            failure.signed_labels = signed_labels;
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error(
                "usage: probe_su2_random_partial_character_cone "
                "SAMPLES MAXIMUM_LEVEL MAXIMUM_FACTORS"
            );
        }
        const std::uint64_t samples = std::stoull(argv[1]);
        const int maximum_level = std::stoi(argv[2]);
        const int maximum_factors = std::stoi(argv[3]);
        if (
            samples == 0U
            || maximum_level < 1
            || maximum_factors < 1
        ) {
            throw std::runtime_error("invalid probe bounds");
        }

        std::atomic<bool> failed{false};
        std::atomic<std::uint64_t> completed{0U};
        Failure first_failure;
        std::mutex failure_mutex;

#pragma omp parallel for schedule(dynamic, 16)
        for (
            std::int64_t raw_sample = 0;
            raw_sample < static_cast<std::int64_t>(samples);
            ++raw_sample
        ) {
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }
            const std::uint64_t sample =
                static_cast<std::uint64_t>(raw_sample);
            SplitMix64 random(
                UINT64_C(0x243f6a8885a308d3)
                ^ (sample * UINT64_C(0x9e3779b97f4a7c15))
            );
            const int level = random.uniform_int(1, maximum_level);
            const int minimum_factors =
                (sample & 1U) == 0U
                    ? 1
                    : std::max(1, maximum_factors / 2);
            const int factors =
                random.uniform_int(minimum_factors, maximum_factors);
            const int regime = static_cast<int>(sample % 4U);
            std::vector<int> labels;
            labels.reserve(static_cast<std::size_t>(factors));
            for (int index = 0; index < factors; ++index) {
                int label = 1;
                if (regime == 0) {
                    label = random.uniform_int(1, level);
                } else if (regime == 1) {
                    const int depth = std::min(level - 1, 8);
                    label = level - random.uniform_int(0, depth);
                } else if (regime == 2) {
                    const int center = (level + 1) / 2;
                    const int radius = std::min(level / 2, 6);
                    label = std::clamp(
                        center + random.uniform_int(-radius, radius),
                        1,
                        level
                    );
                } else {
                    label = random.uniform_int(1, std::min(level, 8));
                }
                labels.push_back(label);
            }

            std::vector<int> signs(
                static_cast<std::size_t>(level + 1),
                0
            );
            for (const int label : labels) {
                int& sign = signs[static_cast<std::size_t>(label)];
                if (sign == 0) {
                    sign = (random() & 1U) == 0U ? 1 : -1;
                }
            }
            int first_positive = 0;
            int first_negative = 0;
            int distinct = 0;
            for (int label = 1; label <= level; ++label) {
                const int sign = signs[static_cast<std::size_t>(label)];
                if (sign == 0) {
                    continue;
                }
                ++distinct;
                if (sign > 0 && first_positive == 0) {
                    first_positive = label;
                }
                if (sign < 0 && first_negative == 0) {
                    first_negative = label;
                }
            }
            if (distinct >= 2 && first_positive == 0) {
                for (int label = 1; label <= level; ++label) {
                    if (signs[static_cast<std::size_t>(label)] != 0) {
                        signs[static_cast<std::size_t>(label)] = 1;
                        break;
                    }
                }
            } else if (distinct >= 2 && first_negative == 0) {
                for (int label = 1; label <= level; ++label) {
                    if (signs[static_cast<std::size_t>(label)] != 0) {
                        signs[static_cast<std::size_t>(label)] = -1;
                        break;
                    }
                }
            }
            std::vector<int> signed_labels;
            signed_labels.reserve(labels.size());
            for (const int label : labels) {
                signed_labels.push_back(
                    signs[static_cast<std::size_t>(label)] * label
                );
            }

            Failure local_failure;
            local_failure.sample = sample;
            if (!inspect(level, signed_labels, local_failure)) {
                bool expected = false;
                if (
                    failed.compare_exchange_strong(
                        expected,
                        true,
                        std::memory_order_relaxed
                    )
                ) {
                    std::lock_guard<std::mutex> lock(failure_mutex);
                    first_failure = std::move(local_failure);
                }
            }
            completed.fetch_add(1U, std::memory_order_relaxed);
        }

        if (failed.load(std::memory_order_relaxed)) {
            std::cout
                << "SU2_RANDOM_PARTIAL_CHARACTER_CONE"
                << " result=FAIL"
                << " sample=" << first_failure.sample
                << " level=" << first_failure.level
                << " target=" << first_failure.target
                << " coefficient=" << first_failure.value
                << " word=[";
            for (std::size_t index = 0U;
                 index < first_failure.signed_labels.size();
                 ++index) {
                if (index != 0U) {
                    std::cout << ',';
                }
                const int value = first_failure.signed_labels[index];
                std::cout << (value < 0 ? '-' : '+') << std::abs(value);
            }
            std::cout << "]\n";
            return EXIT_FAILURE;
        }

        std::cout
            << "SU2_RANDOM_PARTIAL_CHARACTER_CONE"
            << " samples=" << completed.load()
            << " maximum_level=" << maximum_level
            << " maximum_factors=" << maximum_factors
            << " threads=" << omp_get_max_threads()
            << " sign_supports=disjoint"
            << " regimes=uniform,top,mid,low"
            << " result=PASS_NO_COUNTEREXAMPLE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_RANDOM_PARTIAL_CHARACTER_CONE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
