#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using Integer = std::int64_t;
using Vector = std::vector<Integer>;
using Matrix = std::vector<Vector>;

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0
        || parsed > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

unsigned int adaptive_threads() {
    const unsigned int hardware
        = std::max(1U, std::thread::hardware_concurrency());
    double load[1] = {0.0};
    const int read = getloadavg(load, 1);
    const unsigned int occupied = read == 1
        ? static_cast<unsigned int>(
            std::ceil(std::max(0.0, load[0])))
        : 0U;
    return std::max(
        1U,
        hardware > occupied + 2U
            ? hardware - occupied - 2U
            : 1U);
}

Integer at(const Vector& values, int index) {
    if (index < 0 || index >= static_cast<int>(values.size())) {
        return 0;
    }
    return values[static_cast<std::size_t>(index)];
}

Integer slope(const Vector& values, int index) {
    return at(values, index) - at(values, index + 1);
}

Vector reserve_coordinates(const Vector& values, int level) {
    const int period = 2 * level + 2;
    const int maximum = static_cast<int>(values.size()) - 1;
    Vector reserve(static_cast<std::size_t>(level + 1), 0);
    for (int label = 0; label < level; ++label) {
        for (int image = 1; image * period + label <= maximum; ++image) {
            reserve[static_cast<std::size_t>(label)]
                += slope(values, image * period + label)
                  + slope(
                        values,
                        image * period + period - label - 2);
        }
    }
    for (int image = 1; image * period + level <= maximum; ++image) {
        reserve[static_cast<std::size_t>(level)]
            += slope(values, image * period + level);
    }
    return reserve;
}

Matrix fusion_matrix(int level, int factor) {
    Matrix matrix(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0));
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper
            = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            matrix[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return matrix;
}

Matrix matrix_product(const Matrix& left, const Matrix& right) {
    Matrix output(
        left.size(),
        Vector(left.size(), 0));
    for (std::size_t row = 0U; row < left.size(); ++row) {
        for (std::size_t middle = 0U; middle < left.size(); ++middle) {
            for (std::size_t column = 0U;
                 column < left.size();
                 ++column) {
                output[row][column]
                    += left[row][middle] * right[middle][column];
            }
        }
    }
    return output;
}

Matrix reserve_transform(int level, int factor) {
    const Matrix fusion = fusion_matrix(level, factor);
    const Matrix square = matrix_product(fusion, fusion);
    Matrix transform(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0));
    for (int source = 0; source <= level; ++source) {
        for (int row = 0; row <= level; ++row) {
            Integer value = 0;
            for (int column = 0; column <= source; ++column) {
                value += square[static_cast<std::size_t>(row)]
                               [static_cast<std::size_t>(column)];
            }
            if (row < level) {
                for (int column = 0; column <= source; ++column) {
                    value -= square[static_cast<std::size_t>(row + 1)]
                                   [static_cast<std::size_t>(column)];
                }
            }
            transform[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(source)] = value;
        }
    }
    return transform;
}

std::string render(const Vector& values) {
    std::string result = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(values[index]);
    }
    return result + ']';
}

struct Failure {
    int level = -1;
    int factor = -1;
    int label = -1;
    Integer value = 0;
    Integer boundary = 0;
    Vector profile;
    Vector reserve;
    Vector caps;
};

void record_failure(
    Failure& failure,
    int level,
    int factor,
    int label,
    Integer value,
    Integer boundary,
    const Vector& profile,
    const Vector& reserve,
    const Vector& caps) {
    if (failure.level >= 0) {
        return;
    }
    failure = {
        level,
        factor,
        label,
        value,
        boundary,
        profile,
        reserve,
        caps};
}

std::string render_failure(const Failure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " factor=" + std::to_string(failure.factor)
        + " label=" + std::to_string(failure.label)
        + " value=" + std::to_string(failure.value)
        + " boundary=" + std::to_string(failure.boundary)
        + " profile=" + render(failure.profile)
        + " reserve=" + render(failure.reserve)
        + " caps=" + render(failure.caps);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 3;
        const int maximum_coefficient = argc >= 3
            ? positive_argument(argv[2], "maximum_coefficient")
            : 12;
        const int maximum_length = argc >= 4
            ? positive_argument(argv[3], "maximum_length")
            : 20;
        if (argc > 4 || maximum_level < 2 || maximum_length < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim3_exhaustive "
                "[maximum_half_level] [maximum_coefficient] "
                "[maximum_length]");
        }

        std::vector<std::pair<int, int>> prefixes;
        for (int first = 1; first <= maximum_coefficient; ++first) {
            const int second_upper = std::min(
                maximum_coefficient,
                3 * first);
            for (int second = 1; second <= second_upper; ++second) {
                prefixes.emplace_back(first, second);
            }
        }

        const unsigned int threads = std::min(
            adaptive_threads(),
            static_cast<unsigned int>(prefixes.size()));
        std::atomic<std::size_t> next_prefix{0U};
        std::atomic<std::uint64_t> profiles{0U};
        std::atomic<std::uint64_t> admissible_profiles{0U};
        std::atomic<std::uint64_t> cap_checks{0U};
        std::atomic<std::uint64_t> cap_failures{0U};
        std::atomic<std::uint64_t> payment_checks{0U};
        std::atomic<std::uint64_t> payment_failures{0U};
        std::atomic<std::uint64_t> stable_checks{0U};
        std::atomic<std::uint64_t> stable_failures{0U};
        std::mutex failure_mutex;
        Failure first_cap;
        Failure first_payment;
        Failure first_stable;

        auto inspect = [&](const Vector& profile) {
            profiles.fetch_add(1U, std::memory_order_relaxed);
            for (int level = 2; level <= maximum_level; ++level) {
                const int period = 2 * level + 2;
                const Vector reserve
                    = reserve_coordinates(profile, level);
                if (std::any_of(
                        reserve.begin(),
                        reserve.end(),
                        [](Integer value) {
                            return value < 0;
                        })) {
                    continue;
                }
                admissible_profiles.fetch_add(
                    1U,
                    std::memory_order_relaxed);
                Vector suffix(
                    static_cast<std::size_t>(level + 1),
                    0);
                Integer cumulative = 0;
                for (int source = level; source >= 0; --source) {
                    cumulative += reserve[
                        static_cast<std::size_t>(source)];
                    suffix[static_cast<std::size_t>(source)]
                        = cumulative;
                }
                Vector caps(
                    static_cast<std::size_t>(level + 1),
                    0);
                for (int source = 0; source <= level; ++source) {
                    caps[static_cast<std::size_t>(source)]
                        = at(profile, period + source)
                          - at(
                                profile,
                                2 * period - source - 1)
                          + at(
                                profile,
                                2 * period + source);
                    const Integer margin
                        = caps[static_cast<std::size_t>(source)]
                          - suffix[static_cast<std::size_t>(source)];
                    cap_checks.fetch_add(
                        1U,
                        std::memory_order_relaxed);
                    if (margin < 0) {
                        cap_failures.fetch_add(
                            1U,
                            std::memory_order_relaxed);
                        std::lock_guard<std::mutex> lock(
                            failure_mutex);
                        record_failure(
                            first_cap,
                            level,
                            -1,
                            source,
                            margin,
                            caps[static_cast<std::size_t>(source)],
                            profile,
                            reserve,
                            caps);
                    }
                }
                for (int factor = 1;
                     factor <= level / 2;
                     ++factor) {
                    const Matrix transform
                        = reserve_transform(level, factor);
                    for (int label = 0;
                         label < 2 * factor;
                         ++label) {
                        Integer boundary = 0;
                        for (int index
                                 = period - 2 * factor + label;
                             index
                                 <= period + 2 * factor - label - 1;
                             ++index) {
                            boundary += at(profile, index);
                        }
                        Integer payment = boundary;
                        for (int source = 0;
                             source <= level;
                             ++source) {
                            const Integer coefficient
                                = transform[
                                    static_cast<std::size_t>(label)]
                                    [static_cast<std::size_t>(source)];
                            if (coefficient < 0) {
                                payment += coefficient
                                    * reserve[
                                        static_cast<std::size_t>(
                                            source)];
                            }
                        }
                        payment_checks.fetch_add(
                            1U,
                            std::memory_order_relaxed);
                        if (payment < 0) {
                            payment_failures.fetch_add(
                                1U,
                                std::memory_order_relaxed);
                            std::lock_guard<std::mutex> lock(
                                failure_mutex);
                            record_failure(
                                first_payment,
                                level,
                                factor,
                                label,
                                payment,
                                boundary,
                                profile,
                                reserve,
                                caps);
                        }
                        Integer margin = boundary;
                        for (int source = 0;
                             source <= level;
                             ++source) {
                            int weight = 0;
                            if (label < factor) {
                                if (source >= label + 1
                                    && source
                                        <= 2 * factor - label - 1) {
                                    weight = 2;
                                } else if (
                                    source >= 2 * factor - label
                                    && source
                                        <= 2 * factor + label + 1) {
                                    weight = 1;
                                }
                            } else if (
                                source
                                    >= 3 * label - 2 * factor + 2
                                && source
                                    <= 2 * factor + label + 1) {
                                weight = 1;
                            }
                            margin -= weight
                                * caps[
                                    static_cast<std::size_t>(source)];
                        }
                        stable_checks.fetch_add(
                            1U,
                            std::memory_order_relaxed);
                        if (margin < 0) {
                            stable_failures.fetch_add(
                                1U,
                                std::memory_order_relaxed);
                            std::lock_guard<std::mutex> lock(
                                failure_mutex);
                            record_failure(
                                first_stable,
                                level,
                                factor,
                                label,
                                margin,
                                boundary,
                                profile,
                                reserve,
                                caps);
                        }
                    }
                }
            }
        };

        auto worker = [&]() {
            while (true) {
                const std::size_t prefix_index = next_prefix.fetch_add(
                    1U,
                    std::memory_order_relaxed);
                if (prefix_index >= prefixes.size()) {
                    return;
                }
                const auto [first, second] = prefixes[prefix_index];
                Vector profile{
                    static_cast<Integer>(first),
                    static_cast<Integer>(second)};
                auto extend = [&](auto&& self) -> void {
                    inspect(profile);
                    if (static_cast<int>(profile.size())
                        >= maximum_length) {
                        return;
                    }
                    const std::size_t size = profile.size();
                    const Integer previous = profile[size - 1U];
                    const Integer before_previous
                        = profile[size - 2U];
                    const Integer log_concave_upper
                        = previous * previous / before_previous;
                    const Integer index
                        = static_cast<Integer>(size - 1U);
                    const Integer normalized_upper
                        = (2 * index + 3) * previous
                          / (2 * index + 1);
                    const Integer upper = std::min<Integer>(
                        maximum_coefficient,
                        std::min(
                            log_concave_upper,
                            normalized_upper));
                    for (Integer next = 1; next <= upper; ++next) {
                        profile.push_back(next);
                        self(self);
                        profile.pop_back();
                    }
                };
                extend(extend);
            }
        };

        std::vector<std::thread> workers;
        workers.reserve(threads);
        for (unsigned int thread = 0U; thread < threads; ++thread) {
            workers.emplace_back(worker);
        }
        for (std::thread& thread : workers) {
            thread.join();
        }

        std::cout
            << "SU2_AIM3_EXHAUSTIVE"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_coefficient=" << maximum_coefficient
            << " maximum_length=" << maximum_length
            << " threads=" << threads
            << " profiles=" << profiles.load()
            << " admissible_profiles=" << admissible_profiles.load()
            << " cap_checks=" << cap_checks.load()
            << " cap_failures=" << cap_failures.load()
            << " payment_checks=" << payment_checks.load()
            << " payment_failures=" << payment_failures.load()
            << " stable_checks=" << stable_checks.load()
            << " stable_failures=" << stable_failures.load()
            << '\n'
            << "FIRST_CAP_FAILURE "
            << render_failure(first_cap)
            << '\n'
            << "FIRST_PAYMENT_FAILURE "
            << render_failure(first_payment)
            << '\n'
            << "FIRST_STABLE_FAILURE "
            << render_failure(first_stable)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
