#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right); output <= left + right; output += 2) {
        function(output);
    }
}

int ray_entry(int left, int right, int first, int second) {
    const auto on_oriented_ray = [](int row, int column, int start_row, int start_column) {
        if (row < start_row || column < start_column) {
            return false;
        }
        return row - start_row == column - start_column
            && (row - start_row) % 2 == 0;
    };
    int value = on_oriented_ray(left, right, first, second) ? 1 : 0;
    if (first != second && on_oriented_ray(left, right, second, first)) {
        ++value;
    }
    return value;
}

int defect_update_entry(
    int plus_label,
    int left,
    int right,
    int ray_first,
    int ray_second
) {
    int value = 0;
    for_each_output(plus_label, left, [&](int output) {
        value += ray_entry(output, right, ray_first, ray_second);
    });
    for_each_output(plus_label, right, [&](int output) {
        value += ray_entry(left, output, ray_first, ray_second);
    });
    for_each_output(left, right, [&](int output) {
        value -= ray_entry(output, plus_label, ray_first, ray_second);
    });
    return value;
}

int tls_update_entry(
    int plus_label,
    int left,
    int right,
    int ray_first,
    int ray_second
) {
    const int current = defect_update_entry(
        plus_label, left, right, ray_first, ray_second
    );
    if (left == 1 || right == 1) {
        return current;
    }
    return current - defect_update_entry(
        plus_label, left - 2, right - 2, ray_first, ray_second
    );
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: verify_su2_tls_operator MAXIMUM_LABEL");
        }
        const int maximum = std::stoi(argv[1]);
        if (maximum < 1) {
            throw std::runtime_error("MAXIMUM_LABEL must be positive");
        }
        const int maximum_ray_start = 2 * maximum;
        std::uint64_t tested = 0;
        for (int plus_label = 1; plus_label <= maximum; ++plus_label) {
            for (int left = 1; left <= maximum; ++left) {
                for (int right = left; right <= maximum; ++right) {
                    if (plus_label == left || plus_label == right) {
                        continue;
                    }
                    for (int ray_first = 1;
                         ray_first <= maximum_ray_start;
                         ++ray_first) {
                        for (int ray_second = ray_first;
                             ray_second <= maximum_ray_start;
                             ++ray_second) {
                            const int coefficient = tls_update_entry(
                                plus_label,
                                left,
                                right,
                                ray_first,
                                ray_second
                            );
                            ++tested;
                            if (coefficient < 0) {
                                std::cout
                                    << "SU2_TLS_OPERATOR FAIL plus=" << plus_label
                                    << " target=(" << left << ',' << right << ')'
                                    << " ray=(" << ray_first << ',' << ray_second << ')'
                                    << " coefficient=" << coefficient << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_TLS_OPERATOR PASS tested=" << tested
                  << " maximum_label=" << maximum << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
