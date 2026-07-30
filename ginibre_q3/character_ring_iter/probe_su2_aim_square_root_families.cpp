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
#include <unistd.h>
#include <utility>
#include <vector>

#define main probe_su2_aim_log_concave_roots_embedded_main
#include "probe_su2_aim_log_concave_roots.cpp"
#undef main

namespace {

struct RootPaymentFailure {
    int level = -1;
    int shell = -1;
    int factor = -1;
    int radius = -1;
    cpp_int value = 0;
    cpp_int boundary = 0;
    Vector root;
    Vector square;
    Vector reserve;
};

unsigned int root_family_adaptive_threads() {
    const unsigned int hardware
        = std::max(1U, std::thread::hardware_concurrency());
    double load[1] = {0.0};
    const int read = getloadavg(load, 1);
    const unsigned int occupied = read == 1
        ? static_cast<unsigned int>(
            std::ceil(std::max(0.0, load[0])))
        : 0U;
    const unsigned int load_limit = hardware > occupied + 2U
        ? hardware - occupied - 2U
        : 1U;

    const long pages = sysconf(_SC_AVPHYS_PAGES);
    const long page_size = sysconf(_SC_PAGESIZE);
    constexpr std::uint64_t bytes_per_worker
        = UINT64_C(64) * UINT64_C(1024) * UINT64_C(1024);
    unsigned int memory_limit = hardware;
    if (pages > 0 && page_size > 0) {
        const std::uint64_t available
            = static_cast<std::uint64_t>(pages)
              * static_cast<std::uint64_t>(page_size);
        memory_limit = static_cast<unsigned int>(
            std::max<std::uint64_t>(
                1U,
                std::min<std::uint64_t>(
                    hardware,
                    available / bytes_per_worker)));
    }
    return std::max(1U, std::min(load_limit, memory_limit));
}

Vector shifted(const Vector& values, int shift) {
    Vector result(static_cast<std::size_t>(shift), 0);
    result.insert(result.end(), values.begin(), values.end());
    return result;
}

void add_root(
    std::vector<Vector>& roots,
    const Vector& values,
    int maximum_length) {
    const int shift_limit = std::min(
        3,
        maximum_length - static_cast<int>(values.size()));
    for (int shift = 0; shift <= shift_limit; ++shift) {
        const Vector candidate = shifted(values, shift);
        if (!log_concave_interval(candidate)) {
            throw std::runtime_error(
                "constructed root is not log concave");
        }
        roots.push_back(candidate);
    }
}

std::string render_root_failure(const RootPaymentFailure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " shell=" + std::to_string(failure.shell)
        + " factor=" + std::to_string(failure.factor)
        + " radius=" + std::to_string(failure.radius)
        + " value=" + failure.value.convert_to<std::string>()
        + " boundary=" + failure.boundary.convert_to<std::string>()
        + " root=" + render(failure.root)
        + " square=" + render(failure.square)
        + " reserve=" + render(failure.reserve);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_length = argc >= 2
            ? positive_argument(argv[1], "maximum_length")
            : 60;
        const int maximum_level = argc >= 3
            ? positive_argument(argv[2], "maximum_half_level")
            : 30;
        if (argc > 3 || maximum_length < 2 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_aim_square_root_families "
                "[maximum_length] [maximum_half_level]");
        }

        std::vector<Vector> roots;
        for (int length = 1; length <= maximum_length; ++length) {
            Vector constant(static_cast<std::size_t>(length), 1);
            Vector rising(static_cast<std::size_t>(length), 0);
            Vector falling(static_cast<std::size_t>(length), 0);
            Vector tent(static_cast<std::size_t>(length), 0);
            for (int index = 0; index < length; ++index) {
                rising[static_cast<std::size_t>(index)] = index + 1;
                falling[static_cast<std::size_t>(index)]
                    = length - index;
                tent[static_cast<std::size_t>(index)]
                    = std::min(index + 1, length - index);
            }
            add_root(roots, constant, maximum_length);
            add_root(roots, rising, maximum_length);
            add_root(roots, falling, maximum_length);
            add_root(roots, tent, maximum_length);

            std::vector<int> heights{
                1,
                2,
                3,
                std::max(1, length / 4),
                std::max(1, length / 3),
                (length + 1) / 2};
            std::sort(heights.begin(), heights.end());
            heights.erase(
                std::unique(heights.begin(), heights.end()),
                heights.end());
            for (const int height : heights) {
                Vector trapezoid(
                    static_cast<std::size_t>(length),
                    0);
                for (int index = 0; index < length; ++index) {
                    trapezoid[static_cast<std::size_t>(index)]
                        = std::min(
                            height,
                            std::min(index + 1, length - index));
                }
                add_root(roots, trapezoid, maximum_length);
            }

            Vector binomial(static_cast<std::size_t>(length), 1);
            for (int index = 1; index < length; ++index) {
                binomial[static_cast<std::size_t>(index)]
                    = binomial[static_cast<std::size_t>(index - 1)]
                      * (length - index) / index;
            }
            add_root(roots, binomial, maximum_length);

            for (int base = 2; base <= 3; ++base) {
                Vector increasing_power(
                    static_cast<std::size_t>(length),
                    1);
                Vector decreasing_power(
                    static_cast<std::size_t>(length),
                    1);
                for (int index = 1; index < length; ++index) {
                    increasing_power[
                        static_cast<std::size_t>(index)]
                        = increasing_power[
                              static_cast<std::size_t>(index - 1)]
                          * base;
                    decreasing_power[
                        static_cast<std::size_t>(
                            length - index - 1)]
                        = decreasing_power[
                              static_cast<std::size_t>(
                                  length - index)]
                          * base;
                }
                add_root(
                    roots,
                    increasing_power,
                    maximum_length);
                add_root(
                    roots,
                    decreasing_power,
                    maximum_length);
            }
        }
        std::sort(roots.begin(), roots.end());
        roots.erase(
            std::unique(roots.begin(), roots.end()),
            roots.end());

        std::vector<std::vector<Matrix>> transforms(
            static_cast<std::size_t>(maximum_level + 1));
        for (int level = 2; level <= maximum_level; ++level) {
            transforms[static_cast<std::size_t>(level)].resize(
                static_cast<std::size_t>(level / 2 + 1));
            for (int factor = 1;
                 factor <= level / 2;
                 ++factor) {
                transforms[static_cast<std::size_t>(level)]
                          [static_cast<std::size_t>(factor)]
                    = reserve_transform(level, factor);
            }
        }

        const unsigned int threads = std::min(
            root_family_adaptive_threads(),
            static_cast<unsigned int>(roots.size()));
        std::atomic<std::size_t> next_root{0U};
        std::atomic<std::uint64_t> suffixes{0U};
        std::atomic<std::uint64_t> reserve_admissible_suffixes{0U};
        std::atomic<std::uint64_t> payment_checks{0U};
        std::atomic<std::uint64_t> payment_failures{0U};
        std::mutex failure_mutex;
        RootPaymentFailure first_failure;

        auto worker = [&]() {
            while (true) {
                const std::size_t root_index = next_root.fetch_add(
                    1U,
                    std::memory_order_relaxed);
                if (root_index >= roots.size()) {
                    return;
                }
                const Vector& root = roots[root_index];
                const Vector square = ordinary_square(root);
                const int maximum_square
                    = static_cast<int>(square.size()) - 1;
                for (int level = 2;
                     level <= maximum_level;
                     ++level) {
                    const int period = 2 * level + 2;
                    for (int shell = 1;
                         shell * period
                             <= maximum_square + level + 1;
                         ++shell) {
                        suffixes.fetch_add(
                            1U,
                            std::memory_order_relaxed);
                        const Vector suffix
                            = image_suffix(square, level, shell);
                        Vector reserve(
                            static_cast<std::size_t>(level + 1),
                            0);
                        for (int radius = 0;
                             radius < level;
                             ++radius) {
                            reserve[
                                static_cast<std::size_t>(radius)]
                                = suffix[
                                    static_cast<std::size_t>(radius)]
                                  - suffix[
                                    static_cast<std::size_t>(
                                        radius + 1)];
                        }
                        reserve[static_cast<std::size_t>(level)]
                            = suffix[
                                static_cast<std::size_t>(level)];
                        if (std::any_of(
                                reserve.begin(),
                                reserve.end(),
                                [](const cpp_int& value) {
                                    return value < 0;
                                })) {
                            continue;
                        }
                        reserve_admissible_suffixes.fetch_add(
                            1U,
                            std::memory_order_relaxed);
                        for (int factor = 1;
                             factor <= level / 2;
                             ++factor) {
                            const Matrix& transform
                                = transforms[
                                    static_cast<std::size_t>(level)]
                                    [static_cast<std::size_t>(
                                        factor)];
                            for (int radius = 0;
                                 radius <= level;
                                 ++radius) {
                                cpp_int boundary = 0;
                                if (radius < 2 * factor) {
                                    for (int index
                                             = shell * period
                                               - 2 * factor + radius;
                                         index
                                             <= shell * period
                                                + 2 * factor
                                                - radius - 1;
                                         ++index) {
                                        boundary += at(square, index);
                                    }
                                }
                                cpp_int payment = boundary;
                                for (int source = 0;
                                     source <= level;
                                     ++source) {
                                    const long long coefficient
                                        = transform[
                                            static_cast<std::size_t>(
                                                radius)]
                                            [static_cast<std::size_t>(
                                                source)];
                                    if (coefficient < 0) {
                                        payment += coefficient
                                            * reserve[
                                                static_cast<
                                                    std::size_t>(
                                                    source)];
                                    }
                                }
                                payment_checks.fetch_add(
                                    1U,
                                    std::memory_order_relaxed);
                                if (payment >= 0) {
                                    continue;
                                }
                                payment_failures.fetch_add(
                                    1U,
                                    std::memory_order_relaxed);
                                std::lock_guard<std::mutex> lock(
                                    failure_mutex);
                                if (first_failure.level < 0) {
                                    first_failure = {
                                        level,
                                        shell,
                                        factor,
                                        radius,
                                        payment,
                                        boundary,
                                        root,
                                        square,
                                        reserve};
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

        std::cout
            << "SU2_AIM_SQUARE_ROOT_FAMILIES"
            << " maximum_length=" << maximum_length
            << " maximum_level=" << 2 * maximum_level
            << " threads=" << threads
            << " roots=" << roots.size()
            << " suffixes=" << suffixes.load()
            << " reserve_admissible_suffixes="
            << reserve_admissible_suffixes.load()
            << " payment_checks=" << payment_checks.load()
            << " payment_failures=" << payment_failures.load()
            << '\n'
            << "FIRST_PAYMENT_FAILURE "
            << render_root_failure(first_failure)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
