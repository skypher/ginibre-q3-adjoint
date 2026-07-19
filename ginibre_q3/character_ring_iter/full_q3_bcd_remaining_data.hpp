#ifndef GINIBRE_Q3_FULL_Q3_BCD_REMAINING_DATA_HPP
#define GINIBRE_Q3_FULL_Q3_BCD_REMAINING_DATA_HPP

#include <array>
#include <cstddef>

namespace full_q3_bcd_remaining {

enum class TailMethod : unsigned char {
    polynomial,
    rational_cap,
    directed_interval
};

struct RowCutoff {
    char family;
    int rank;
    int tail_onset;
    int moment_through;
    TailMethod tail_method;
};

// The stable prefix ends at 2r+1 in types B/C and at r-1 in type D.
// Exact hierarchy values are checked through tail_onset-2.  The moment
// endpoint can be larger only when an exact polynomial tail consumes it.
inline constexpr std::array<RowCutoff, 70> row_cutoffs{{
    {'B', 2, 47, 45, TailMethod::rational_cap},
    {'B', 3, 37, 35, TailMethod::polynomial},
    {'B', 4, 41, 42, TailMethod::polynomial},
    {'B', 5, 55, 53, TailMethod::polynomial},
    {'B', 6, 51, 49, TailMethod::rational_cap},
    {'B', 7, 53, 51, TailMethod::directed_interval},
    {'B', 8, 51, 49, TailMethod::directed_interval},
    {'B', 9, 53, 51, TailMethod::directed_interval},
    {'B', 10, 57, 55, TailMethod::directed_interval},
    {'B', 11, 57, 55, TailMethod::directed_interval},
    {'B', 12, 61, 59, TailMethod::directed_interval},
    {'B', 13, 61, 59, TailMethod::directed_interval},
    {'B', 14, 45, 43, TailMethod::directed_interval},
    {'B', 15, 45, 43, TailMethod::directed_interval},
    {'B', 16, 45, 43, TailMethod::directed_interval},
    {'B', 17, 45, 43, TailMethod::directed_interval},

    {'C', 2, 47, 45, TailMethod::rational_cap},
    {'C', 3, 37, 35, TailMethod::polynomial},
    {'C', 4, 61, 59, TailMethod::directed_interval},
    {'C', 5, 55, 53, TailMethod::directed_interval},
    {'C', 6, 39, 37, TailMethod::directed_interval},
    {'C', 7, 39, 37, TailMethod::directed_interval},
    {'C', 8, 49, 47, TailMethod::directed_interval},
    {'C', 9, 49, 47, TailMethod::directed_interval},
    {'C', 10, 53, 51, TailMethod::directed_interval},
    {'C', 11, 55, 53, TailMethod::directed_interval},
    {'C', 12, 53, 51, TailMethod::directed_interval},
    {'C', 13, 57, 55, TailMethod::directed_interval},
    {'C', 14, 59, 57, TailMethod::directed_interval},
    {'C', 15, 55, 53, TailMethod::directed_interval},
    {'C', 16, 55, 53, TailMethod::directed_interval},
    {'C', 17, 57, 55, TailMethod::directed_interval},
    {'C', 18, 45, 43, TailMethod::directed_interval},
    {'C', 19, 45, 43, TailMethod::directed_interval},
    {'C', 20, 45, 43, TailMethod::directed_interval},
    {'C', 21, 47, 45, TailMethod::directed_interval},
    {'C', 22, 47, 45, TailMethod::directed_interval},
    {'C', 23, 49, 47, TailMethod::directed_interval},
    {'C', 24, 51, 49, TailMethod::directed_interval},
    {'C', 25, 53, 51, TailMethod::directed_interval},
    {'C', 26, 55, 53, TailMethod::directed_interval},
    {'C', 27, 57, 55, TailMethod::directed_interval},
    {'C', 28, 59, 57, TailMethod::directed_interval},

    {'D', 4, 47, 45, TailMethod::polynomial},
    {'D', 5, 33, 31, TailMethod::polynomial},
    {'D', 6, 49, 47, TailMethod::rational_cap},
    {'D', 7, 53, 51, TailMethod::rational_cap},
    {'D', 8, 57, 55, TailMethod::directed_interval},
    {'D', 9, 61, 59, TailMethod::rational_cap},
    {'D', 10, 51, 49, TailMethod::directed_interval},
    {'D', 11, 57, 55, TailMethod::directed_interval},
    {'D', 12, 45, 43, TailMethod::directed_interval},
    {'D', 13, 45, 43, TailMethod::directed_interval},
    {'D', 14, 45, 43, TailMethod::directed_interval},
    {'D', 15, 45, 43, TailMethod::directed_interval},
    {'D', 16, 43, 41, TailMethod::directed_interval},
    {'D', 17, 43, 41, TailMethod::directed_interval},
    {'D', 18, 43, 41, TailMethod::directed_interval},
    {'D', 19, 45, 43, TailMethod::directed_interval},
    {'D', 20, 45, 43, TailMethod::directed_interval},
    {'D', 21, 45, 43, TailMethod::directed_interval},
    {'D', 22, 45, 43, TailMethod::directed_interval},
    {'D', 23, 45, 43, TailMethod::directed_interval},
    {'D', 24, 47, 45, TailMethod::directed_interval},
    {'D', 25, 45, 43, TailMethod::directed_interval},
    {'D', 26, 47, 45, TailMethod::directed_interval},
    {'D', 27, 47, 45, TailMethod::directed_interval},
    {'D', 28, 49, 47, TailMethod::directed_interval},
    {'D', 29, 47, 45, TailMethod::directed_interval},
    {'D', 30, 49, 47, TailMethod::directed_interval},
}};

inline const RowCutoff* find_row(char family, int rank) {
    for (const RowCutoff& row : row_cutoffs) {
        if (row.family == family && row.rank == rank) return &row;
    }
    return nullptr;
}

inline constexpr int required_maximum_moment = 59;
inline constexpr std::size_t polynomial_rows = 6;
inline constexpr std::size_t rational_cap_rows = 6;
inline constexpr std::size_t directed_interval_rows = 58;
inline constexpr std::size_t full_residual_pairs = 12993;

}  // namespace full_q3_bcd_remaining

#endif
