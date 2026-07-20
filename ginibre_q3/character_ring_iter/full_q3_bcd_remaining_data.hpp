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
inline constexpr std::array<RowCutoff, 114> row_cutoffs{{
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
    {'B', 18, 47, 45, TailMethod::directed_interval},
    {'B', 19, 45, 43, TailMethod::directed_interval},
    {'B', 20, 45, 43, TailMethod::directed_interval},
    {'B', 21, 45, 43, TailMethod::directed_interval},

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
    {'D', 31, 49, 47, TailMethod::directed_interval},
    {'D', 32, 51, 49, TailMethod::directed_interval},
    {'D', 33, 51, 49, TailMethod::directed_interval},
    {'D', 34, 51, 49, TailMethod::directed_interval},
    {'D', 35, 51, 49, TailMethod::directed_interval},
    {'D', 36, 53, 51, TailMethod::directed_interval},
    {'D', 37, 53, 51, TailMethod::directed_interval},
    {'D', 38, 55, 53, TailMethod::directed_interval},
    {'D', 39, 53, 51, TailMethod::directed_interval},
    {'D', 40, 55, 53, TailMethod::directed_interval},
    {'D', 41, 55, 53, TailMethod::directed_interval},
    {'D', 42, 57, 55, TailMethod::directed_interval},
    {'D', 43, 55, 53, TailMethod::directed_interval},
    {'D', 44, 57, 55, TailMethod::directed_interval},
    {'D', 45, 57, 55, TailMethod::directed_interval},
    {'D', 46, 59, 57, TailMethod::directed_interval},
    {'D', 47, 59, 57, TailMethod::directed_interval},
    {'D', 48, 61, 59, TailMethod::directed_interval},
    {'D', 49, 61, 59, TailMethod::directed_interval},
    {'D', 50, 61, 59, TailMethod::directed_interval},
    {'D', 51, 61, 59, TailMethod::directed_interval},
    {'D', 52, 63, 61, TailMethod::directed_interval},
    {'D', 53, 61, 59, TailMethod::directed_interval},
    {'D', 54, 63, 61, TailMethod::directed_interval},
    {'D', 55, 63, 61, TailMethod::directed_interval},
    {'D', 56, 65, 63, TailMethod::directed_interval},
    {'D', 57, 63, 61, TailMethod::directed_interval},
    {'D', 58, 65, 63, TailMethod::directed_interval},
    {'D', 59, 65, 63, TailMethod::directed_interval},
    {'D', 60, 67, 65, TailMethod::directed_interval},
    {'D', 61, 67, 65, TailMethod::directed_interval},
    {'D', 62, 69, 67, TailMethod::directed_interval},
    {'D', 63, 67, 65, TailMethod::directed_interval},
    {'D', 64, 69, 67, TailMethod::directed_interval},
    {'D', 65, 69, 67, TailMethod::directed_interval},
    {'D', 66, 71, 69, TailMethod::directed_interval},
    {'D', 67, 71, 69, TailMethod::directed_interval},
    {'D', 68, 71, 69, TailMethod::directed_interval},
    {'D', 69, 71, 69, TailMethod::directed_interval},
    {'D', 70, 73, 71, TailMethod::directed_interval},
}};

inline const RowCutoff* find_row(char family, int rank) {
    for (const RowCutoff& row : row_cutoffs) {
        if (row.family == family && row.rank == rank) return &row;
    }
    return nullptr;
}

inline constexpr int required_maximum_moment = 71;
inline constexpr std::size_t polynomial_rows = 6;
inline constexpr std::size_t rational_cap_rows = 6;
inline constexpr std::size_t directed_interval_rows = 102;
inline constexpr std::size_t full_residual_pairs = 17862;

// The directed-MPFR construction was proved for the original low-tail
// subledger.  Later B18--B21 and D31--D70 rows are determinant-only rows and
// must not silently enlarge that verifier's analytic scope.
inline constexpr int mpfr_b_rank_end = 17;
inline constexpr int mpfr_c_rank_end = 28;
inline constexpr int mpfr_d_rank_end = 30;
inline constexpr std::size_t mpfr_directed_interval_rows = 58;

}  // namespace full_q3_bcd_remaining

#endif
