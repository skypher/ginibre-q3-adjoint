#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;
using SmallVector = std::vector<long long>;
using SmallMatrix = std::vector<SmallVector>;

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

cpp_int at(const Vector& values, int index) {
    if (index < 0 || index >= static_cast<int>(values.size())) {
        return 0;
    }
    return values[static_cast<std::size_t>(index)];
}

cpp_int greatest_common_divisor(cpp_int left, cpp_int right) {
    while (right != 0) {
        const cpp_int remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

cpp_int slope(const Vector& values, int index) {
    return at(values, index)
        - at(values, index + 1);
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

SmallMatrix fusion_matrix(int level, int factor) {
    SmallMatrix matrix(
        static_cast<std::size_t>(level + 1),
        SmallVector(static_cast<std::size_t>(level + 1), 0));
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

SmallMatrix matrix_product(
    const SmallMatrix& left,
    const SmallMatrix& right) {
    SmallMatrix output(
        left.size(),
        SmallVector(left.size(), 0));
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

SmallMatrix reserve_transform(int level, int factor) {
    const SmallMatrix fusion = fusion_matrix(level, factor);
    const SmallMatrix square = matrix_product(fusion, fusion);
    SmallMatrix transform(
        static_cast<std::size_t>(level + 1),
        SmallVector(static_cast<std::size_t>(level + 1), 0));
    for (int source = 0; source <= level; ++source) {
        for (int row = 0; row <= level; ++row) {
            long long value = 0;
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
        result += values[index].convert_to<std::string>();
    }
    return result + ']';
}

std::string render(const SmallVector& values) {
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
    int sample = -1;
    int level = -1;
    int factor = -1;
    int label = -1;
    cpp_int value = 0;
    cpp_int boundary = 0;
    Vector profile;
    Vector reserve;
    SmallVector transform_row;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        const int samples = argc >= 2
            ? positive_argument(argv[1], "samples")
            : 100000;
        const int maximum_level = argc >= 3
            ? positive_argument(argv[2], "maximum_half_level")
            : 10;
        const int maximum_length = argc >= 4
            ? positive_argument(argv[3], "maximum_length")
            : 40;
        const int maximum_denominator = argc >= 5
            ? positive_argument(argv[4], "maximum_denominator")
            : 64;
        const bool normalized_mode = argc >= 6
            && std::string(argv[5]) == "normalized";
        if (argc > 6 || (argc >= 6 && !normalized_mode)
            || maximum_level < 2
            || maximum_length < 2 || maximum_denominator < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_aim1_abstract_payment "
                "[samples] [maximum_half_level] "
                "[maximum_length] [maximum_denominator] [normalized]");
        }

        const unsigned int threads = adaptive_threads();
        std::atomic<int> next_sample{0};
        std::atomic<std::uint64_t> admissible_profiles{0U};
        std::atomic<std::uint64_t> payment_checks{0U};
        std::atomic<std::uint64_t> payment_failures{0U};
        std::atomic<std::uint64_t> normalized_profiles{0U};
        std::atomic<std::uint64_t> normalized_payment_checks{0U};
        std::atomic<std::uint64_t> normalized_payment_failures{0U};
        std::atomic<std::uint64_t> all_three_term_cap_checks{0U};
        std::atomic<std::uint64_t> all_three_term_cap_failures{0U};
        std::atomic<std::uint64_t> three_term_cap_checks{0U};
        std::atomic<std::uint64_t> three_term_cap_failures{0U};
        std::atomic<std::uint64_t> three_term_bound_checks{0U};
        std::atomic<std::uint64_t> three_term_bound_failures{0U};
        std::atomic<std::uint64_t> stable_bound_checks{0U};
        std::atomic<std::uint64_t> stable_bound_failures{0U};
        std::mutex failure_mutex;
        Failure first;
        Failure first_normalized;
        Failure first_all_three_term_cap;
        Failure first_three_term_cap;
        Failure first_three_term_bound;
        Failure first_stable_bound;

        auto worker = [&]() {
            while (true) {
                const int sample
                    = next_sample.fetch_add(1, std::memory_order_relaxed);
                if (sample >= samples) {
                    return;
                }
                std::mt19937_64 generator(
                    UINT64_C(0x6a09e667f3bcc909)
                    ^ (static_cast<std::uint64_t>(sample)
                       * UINT64_C(0xbb67ae8584caa73b)));
                const int length
                    = 2 + static_cast<int>(
                        generator()
                        % static_cast<std::uint64_t>(
                            maximum_length - 1));
                const int denominator
                    = 2 + static_cast<int>(
                        generator()
                        % static_cast<std::uint64_t>(
                            maximum_denominator - 1));
                std::vector<int> numerators(
                    static_cast<std::size_t>(length - 1));
                for (int& numerator : numerators) {
                    numerator = 1 + static_cast<int>(
                        generator()
                        % static_cast<std::uint64_t>(2 * denominator));
                }
                std::sort(
                    numerators.begin(),
                    numerators.end(),
                    std::greater<int>());

                Vector profile(static_cast<std::size_t>(length), 0);
                if (normalized_mode) {
                    std::vector<int> normalized(
                        static_cast<std::size_t>(length),
                        0);
                    for (int& value : normalized) {
                        value = 1 + static_cast<int>(
                            generator()
                            % static_cast<std::uint64_t>(
                                maximum_denominator));
                    }
                    std::sort(
                        normalized.begin(),
                        normalized.end(),
                        std::greater<int>());
                    for (int index = 0; index < length; ++index) {
                        profile[static_cast<std::size_t>(index)]
                            = (2 * index + 1)
                              * normalized[
                                  static_cast<std::size_t>(index)];
                    }
                } else {
                    profile[0] = 1;
                    for (int power = 1; power < length; ++power) {
                        profile[0] *= denominator;
                    }
                    for (int index = 1; index < length; ++index) {
                        profile[static_cast<std::size_t>(index)]
                            = profile[
                                  static_cast<std::size_t>(index - 1)]
                              * numerators[
                                  static_cast<std::size_t>(index - 1)]
                              / denominator;
                    }
                }
                cpp_int divisor = 0;
                for (const cpp_int& value : profile) {
                    divisor = greatest_common_divisor(divisor, value);
                }
                if (divisor > 1) {
                    for (cpp_int& value : profile) {
                        value /= divisor;
                    }
                }
                bool dimension_normalized = true;
                for (int index = 0; index + 1 < length; ++index) {
                    if ((2 * index + 3)
                            * profile[static_cast<std::size_t>(index)]
                        < (2 * index + 1)
                            * profile[
                                static_cast<std::size_t>(index + 1)]) {
                        dimension_normalized = false;
                        break;
                    }
                }

                for (int level = 2; level <= maximum_level; ++level) {
                    const Vector reserve
                        = reserve_coordinates(profile, level);
                    const int period = 2 * level + 2;
                    if (std::any_of(
                            reserve.begin(),
                            reserve.end(),
                            [](const cpp_int& value) {
                                return value < 0;
                            })) {
                        continue;
                    }
                    admissible_profiles.fetch_add(
                        1U,
                        std::memory_order_relaxed);
                    for (int source = 0; source <= level; ++source) {
                        cpp_int suffix = 0;
                        for (int coordinate = source;
                             coordinate <= level;
                             ++coordinate) {
                            suffix += reserve[
                                static_cast<std::size_t>(coordinate)];
                        }
                        const cpp_int cap
                            = at(profile, period + source)
                              - at(
                                  profile,
                                  2 * period - source - 1)
                              + at(
                                  profile,
                                  2 * period + source);
                        const cpp_int margin = cap - suffix;
                        all_three_term_cap_checks.fetch_add(
                            1U,
                            std::memory_order_relaxed);
                        if (margin < 0) {
                            all_three_term_cap_failures.fetch_add(
                                1U,
                                std::memory_order_relaxed);
                            std::lock_guard<std::mutex> lock(
                                failure_mutex);
                            if (first_all_three_term_cap.sample < 0
                                || sample
                                    < first_all_three_term_cap.sample) {
                                first_all_three_term_cap = {
                                    sample,
                                    level,
                                    -1,
                                    source,
                                    margin,
                                    cap,
                                    profile,
                                    reserve,
                                    {}};
                            }
                        }
                    }
                    if (dimension_normalized) {
                        normalized_profiles.fetch_add(
                            1U,
                            std::memory_order_relaxed);
                        for (int source = 0;
                             source <= level;
                             ++source) {
                            cpp_int suffix = 0;
                            for (int coordinate = source;
                                 coordinate <= level;
                                 ++coordinate) {
                                suffix += reserve[
                                    static_cast<std::size_t>(
                                        coordinate)];
                            }
                            const cpp_int cap
                                = at(profile, period + source)
                                  - at(
                                      profile,
                                      2 * period - source - 1)
                                  + at(
                                      profile,
                                      2 * period + source);
                            const cpp_int margin = cap - suffix;
                            three_term_cap_checks.fetch_add(
                                1U,
                                std::memory_order_relaxed);
                            if (margin < 0) {
                                three_term_cap_failures.fetch_add(
                                    1U,
                                    std::memory_order_relaxed);
                                std::lock_guard<std::mutex> lock(
                                    failure_mutex);
                                if (first_three_term_cap.sample < 0
                                    || sample
                                        < first_three_term_cap.sample) {
                                    first_three_term_cap = {
                                        sample,
                                        level,
                                        -1,
                                        source,
                                        margin,
                                        cap,
                                        profile,
                                        reserve,
                                        {}};
                                }
                            }
                        }
                    }
                    for (int factor = 1;
                         factor <= level / 2;
                         ++factor) {
                        const SmallMatrix transform
                            = reserve_transform(level, factor);
                        for (int label = 0; label <= level; ++label) {
                            cpp_int boundary = 0;
                            if (label < 2 * factor) {
                                for (int index
                                         = period - 2 * factor + label;
                                     index
                                         <= period + 2 * factor - label - 1;
                                     ++index) {
                                    boundary += at(profile, index);
                                }
                            }
                            if (dimension_normalized
                                && label < 2 * factor) {
                                cpp_int bound = boundary;
                                cpp_int stable_bound = boundary;
                                long long previous_load = 0;
                                for (int source = 0;
                                     source <= level;
                                     ++source) {
                                    const long long load = std::max(
                                        0LL,
                                        -transform[
                                            static_cast<std::size_t>(
                                                label)]
                                            [static_cast<std::size_t>(
                                                source)]);
                                    const cpp_int cap
                                        = at(
                                            profile,
                                            period + source)
                                          - at(
                                              profile,
                                              2 * period - source - 1)
                                          + at(
                                              profile,
                                              2 * period + source);
                                    bound -= (load - previous_load) * cap;
                                    long long stable_weight = 0;
                                    if (label < factor) {
                                        if (source >= label + 1
                                            && source
                                                <= 2 * factor
                                                   - label - 1) {
                                            stable_weight = 2;
                                        } else if (
                                            source
                                                >= 2 * factor - label
                                            && source
                                                <= 2 * factor
                                                   + label + 1) {
                                            stable_weight = 1;
                                        }
                                    } else if (
                                        source
                                            >= 3 * label
                                               - 2 * factor + 2
                                        && source
                                            <= 2 * factor
                                               + label + 1) {
                                        stable_weight = 1;
                                    }
                                    stable_bound -= stable_weight * cap;
                                    previous_load = load;
                                }
                                three_term_bound_checks.fetch_add(
                                    1U,
                                    std::memory_order_relaxed);
                                if (bound < 0) {
                                    three_term_bound_failures.fetch_add(
                                        1U,
                                        std::memory_order_relaxed);
                                    std::lock_guard<std::mutex> lock(
                                        failure_mutex);
                                    if (first_three_term_bound.sample < 0
                                        || sample
                                            < first_three_term_bound.sample) {
                                        first_three_term_bound = {
                                            sample,
                                            level,
                                            factor,
                                            label,
                                            bound,
                                            boundary,
                                            profile,
                                            reserve,
                                            transform[
                                                static_cast<std::size_t>(
                                                    label)]};
                                    }
                                }
                                stable_bound_checks.fetch_add(
                                    1U,
                                    std::memory_order_relaxed);
                                if (stable_bound < 0) {
                                    stable_bound_failures.fetch_add(
                                        1U,
                                        std::memory_order_relaxed);
                                    std::lock_guard<std::mutex> lock(
                                        failure_mutex);
                                    if (first_stable_bound.sample < 0
                                        || sample
                                            < first_stable_bound.sample) {
                                        first_stable_bound = {
                                            sample,
                                            level,
                                            factor,
                                            label,
                                            stable_bound,
                                            boundary,
                                            profile,
                                            reserve,
                                            transform[
                                                static_cast<std::size_t>(
                                                    label)]};
                                    }
                                }
                            }
                            cpp_int payment = boundary;
                            for (int source = 0;
                                 source <= level;
                                 ++source) {
                                const long long coefficient
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
                            if (dimension_normalized) {
                                normalized_payment_checks.fetch_add(
                                    1U,
                                    std::memory_order_relaxed);
                            }
                            if (payment < 0) {
                                payment_failures.fetch_add(
                                    1U,
                                    std::memory_order_relaxed);
                                std::lock_guard<std::mutex> lock(
                                    failure_mutex);
                                if (first.sample < 0
                                    || sample < first.sample) {
                                    first = {
                                        sample,
                                        level,
                                        factor,
                                        label,
                                        payment,
                                        boundary,
                                        profile,
                                        reserve,
                                        transform[
                                            static_cast<std::size_t>(
                                                label)]};
                                }
                                if (dimension_normalized) {
                                    normalized_payment_failures.fetch_add(
                                        1U,
                                        std::memory_order_relaxed);
                                    if (first_normalized.sample < 0
                                        || sample
                                            < first_normalized.sample) {
                                        first_normalized = {
                                            sample,
                                            level,
                                            factor,
                                            label,
                                            payment,
                                            boundary,
                                            profile,
                                            reserve,
                                            transform[
                                                static_cast<std::size_t>(
                                                    label)]};
                                    }
                                }
                            }
                        }
                    }
                }
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

        Vector stable_without_reserve_profile(20U, 0);
        for (std::size_t index = 0U;
             index < stable_without_reserve_profile.size();
             ++index) {
            stable_without_reserve_profile[index] = 2 * index + 1;
        }
        const Vector stable_without_reserve_reserve
            = reserve_coordinates(stable_without_reserve_profile, 3);
        Vector stable_without_reserve_caps(4U, 0);
        for (int source = 0; source <= 3; ++source) {
            stable_without_reserve_caps[
                static_cast<std::size_t>(source)]
                = at(stable_without_reserve_profile, 8 + source)
                  - at(
                        stable_without_reserve_profile,
                        16 - source - 1)
                  + at(
                        stable_without_reserve_profile,
                        16 + source);
        }
        cpp_int stable_without_reserve_boundary = 0;
        for (int index = 6; index <= 9; ++index) {
            stable_without_reserve_boundary
                += at(stable_without_reserve_profile, index);
        }
        const cpp_int stable_without_reserve_margin
            = stable_without_reserve_boundary
              - 2 * stable_without_reserve_caps[1]
              - stable_without_reserve_caps[2]
              - stable_without_reserve_caps[3];
        const Vector expected_reserve{-6, -6, -6, 37};
        const Vector expected_caps{19, 25, 31, 37};
        if (stable_without_reserve_boundary != 64
            || stable_without_reserve_margin != -54
            || stable_without_reserve_reserve != expected_reserve
            || stable_without_reserve_caps != expected_caps) {
            throw std::runtime_error(
                "stable-without-reserve counterexample mismatch");
        }

        std::cout
            << "SU2_AIM1_ABSTRACT_PAYMENT"
            << " samples=" << samples
            << " maximum_level=" << 2 * maximum_level
            << " maximum_length=" << maximum_length
            << " maximum_denominator=" << maximum_denominator
            << " mode="
            << (normalized_mode ? "normalized" : "log_concave")
            << " threads=" << threads
            << " admissible_profiles=" << admissible_profiles.load()
            << " payment_checks=" << payment_checks.load()
            << " payment_failures=" << payment_failures.load()
            << " normalized_profiles=" << normalized_profiles.load()
            << " normalized_payment_checks="
            << normalized_payment_checks.load()
            << " normalized_payment_failures="
            << normalized_payment_failures.load()
            << " all_three_term_cap_checks="
            << all_three_term_cap_checks.load()
            << " all_three_term_cap_failures="
            << all_three_term_cap_failures.load()
            << " three_term_cap_checks="
            << three_term_cap_checks.load()
            << " three_term_cap_failures="
            << three_term_cap_failures.load()
            << " three_term_bound_checks="
            << three_term_bound_checks.load()
            << " three_term_bound_failures="
            << three_term_bound_failures.load()
            << " stable_bound_checks="
            << stable_bound_checks.load()
            << " stable_bound_failures="
            << stable_bound_failures.load()
            << '\n'
            << "FIRST_PAYMENT_FAILURE"
            << " sample=" << first.sample
            << " level="
            << (first.level < 0 ? -1 : 2 * first.level)
            << " factor=" << first.factor
            << " label=" << first.label
            << " value=" << first.value
            << " boundary=" << first.boundary
            << " profile=" << render(first.profile)
            << " reserve=" << render(first.reserve)
            << " transform_row=" << render(first.transform_row)
            << '\n'
            << "FIRST_NORMALIZED_PAYMENT_FAILURE"
            << " sample=" << first_normalized.sample
            << " level="
            << (first_normalized.level < 0
                    ? -1
                    : 2 * first_normalized.level)
            << " factor=" << first_normalized.factor
            << " label=" << first_normalized.label
            << " value=" << first_normalized.value
            << " boundary=" << first_normalized.boundary
            << " profile=" << render(first_normalized.profile)
            << " reserve=" << render(first_normalized.reserve)
            << " transform_row="
            << render(first_normalized.transform_row)
            << '\n'
            << "FIRST_ALL_THREE_TERM_CAP_FAILURE"
            << " sample=" << first_all_three_term_cap.sample
            << " level="
            << (first_all_three_term_cap.level < 0
                    ? -1
                    : 2 * first_all_three_term_cap.level)
            << " factor=" << first_all_three_term_cap.factor
            << " label=" << first_all_three_term_cap.label
            << " value=" << first_all_three_term_cap.value
            << " boundary=" << first_all_three_term_cap.boundary
            << " profile=" << render(first_all_three_term_cap.profile)
            << " reserve=" << render(first_all_three_term_cap.reserve)
            << " transform_row="
            << render(first_all_three_term_cap.transform_row)
            << '\n'
            << "FIRST_THREE_TERM_CAP_FAILURE"
            << " sample=" << first_three_term_cap.sample
            << " level="
            << (first_three_term_cap.level < 0
                    ? -1
                    : 2 * first_three_term_cap.level)
            << " factor=" << first_three_term_cap.factor
            << " label=" << first_three_term_cap.label
            << " value=" << first_three_term_cap.value
            << " boundary=" << first_three_term_cap.boundary
            << " profile=" << render(first_three_term_cap.profile)
            << " reserve=" << render(first_three_term_cap.reserve)
            << " transform_row="
            << render(first_three_term_cap.transform_row)
            << '\n'
            << "FIRST_THREE_TERM_BOUND_FAILURE"
            << " sample=" << first_three_term_bound.sample
            << " level="
            << (first_three_term_bound.level < 0
                    ? -1
                    : 2 * first_three_term_bound.level)
            << " factor=" << first_three_term_bound.factor
            << " label=" << first_three_term_bound.label
            << " value=" << first_three_term_bound.value
            << " boundary=" << first_three_term_bound.boundary
            << " profile=" << render(first_three_term_bound.profile)
            << " reserve=" << render(first_three_term_bound.reserve)
            << " transform_row="
            << render(first_three_term_bound.transform_row)
            << '\n'
            << "FIRST_STABLE_BOUND_FAILURE"
            << " sample=" << first_stable_bound.sample
            << " level="
            << (first_stable_bound.level < 0
                    ? -1
                    : 2 * first_stable_bound.level)
            << " factor=" << first_stable_bound.factor
            << " label=" << first_stable_bound.label
            << " value=" << first_stable_bound.value
            << " boundary=" << first_stable_bound.boundary
            << " profile=" << render(first_stable_bound.profile)
            << " reserve=" << render(first_stable_bound.reserve)
            << " transform_row="
            << render(first_stable_bound.transform_row)
            << '\n'
            << "STABLE_WITHOUT_RESERVE_COUNTEREXAMPLE"
            << " level=6"
            << " factor=1"
            << " label=0"
            << " value=" << stable_without_reserve_margin
            << " boundary=" << stable_without_reserve_boundary
            << " profile="
            << render(stable_without_reserve_profile)
            << " reserve="
            << render(stable_without_reserve_reserve)
            << " caps=" << render(stable_without_reserve_caps)
            << '\n'
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
