#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

namespace {

long double eigenvalue(const int level, const int label, const int mode) {
    const long double pi = std::acos(-1.0L);
    const long double theta = static_cast<long double>(mode + 1) * pi
        / static_cast<long double>(level + 2);
    return std::sin(static_cast<long double>(label + 1) * theta)
        / std::sin(theta);
}

long double eigenvector(const int level, const int label, const int mode) {
    const long double pi = std::acos(-1.0L);
    const long double theta = static_cast<long double>(mode + 1) * pi
        / static_cast<long double>(level + 2);
    return std::sqrt(2.0L / static_cast<long double>(level + 2))
        * std::sin(static_cast<long double>(label + 1) * theta);
}

long double wedge_coordinate(
    const int level,
    const int first,
    const int second,
    const int left_mode,
    const int right_mode
) {
    return eigenvector(level, first, left_mode)
            * eigenvector(level, second, right_mode)
        - eigenvector(level, first, right_mode)
            * eigenvector(level, second, left_mode);
}

using Orbit = std::pair<int, int>;

Orbit canonical_orbit(const int level, const int first, const int second) {
    Orbit answer{level + 1, level + 1};
    for (const int left : {first, level - first}) {
        for (const int right : {second, level - second}) {
            const Orbit candidate{
                std::min(left, right), std::max(left, right)
            };
            answer = std::min(answer, candidate);
        }
    }
    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_two_odd_spectral_orbits MAXIMUM_LEVEL"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        if (maximum_level < 3) {
            throw std::runtime_error("maximum level must be at least three");
        }
        constexpr long double tolerance = 1.0e-14L;
        long long cases = 0;
        long long orbits = 0;
        bool found_negative_weight = false;
        bool found_negative_coordinate = false;
        for (int level = 3; level <= maximum_level; ++level) {
            for (int q = 3; q <= level; q += 2) {
                for (int r = q; r <= level; r += 2) {
                    for (int a = 1; a <= level; a += 2) {
                        if (a == q || a == r) {
                            continue;
                        }
                        if (level < q + r + a + 1) {
                            continue;
                        }
                        if (a + 1 < q + r) {
                            continue;
                        }
                        ++cases;
                        std::map<Orbit, long double> weights;
                        for (int first = 0; first <= level; ++first) {
                            for (int second = first + 1;
                                 second <= level;
                                 ++second) {
                                const long double weight = wedge_coordinate(
                                    level, a, 0, first, second
                                ) * wedge_coordinate(
                                    level, 1, 0, first, second
                                ) * (eigenvalue(level, q, first)
                                    + eigenvalue(level, q, second))
                                  * (eigenvalue(level, r, first)
                                    + eigenvalue(level, r, second));
                                weights[canonical_orbit(
                                    level, first, second
                                )] += weight;
                            }
                        }
                        for (const auto& [orbit, weight] : weights) {
                            ++orbits;
                            if (weight < -tolerance
                                && !found_negative_weight) {
                                found_negative_weight = true;
                                std::cout << std::setprecision(18)
                                          << "NEGATIVE_ORBIT_WEIGHT level="
                                          << level << " q=" << q << " r=" << r
                                          << " a=" << a << " orbit=("
                                          << orbit.first << ',' << orbit.second
                                          << ") weight=" << weight << '\n';
                            }
                            if (std::abs(weight) <= tolerance) {
                                continue;
                            }
                            for (int even = 2; even <= level; even += 2) {
                                const long double coordinate =
                                    eigenvalue(level, even, orbit.first)
                                    + eigenvalue(level, even, orbit.second);
                                if (coordinate < -tolerance
                                    && !found_negative_coordinate) {
                                    found_negative_coordinate = true;
                                    std::cout << std::setprecision(18)
                                              << "ACTIVE_NEGATIVE_COORDINATE level="
                                              << level << " q=" << q
                                              << " r=" << r << " a=" << a
                                              << " orbit=(" << orbit.first << ','
                                              << orbit.second << ") even=" << even
                                              << " orbit_weight=" << weight
                                              << " coordinate=" << coordinate
                                              << '\n';
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_TWO_ODD_SPECTRAL_ORBITS cases=" << cases
                  << " orbits=" << orbits
                  << " nonnegative_orbit_weights="
                  << (found_negative_weight ? "false" : "true")
                  << " nonnegative_active_coordinates="
                  << (found_negative_coordinate ? "false" : "true")
                  << " maximum_level=" << maximum_level << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
