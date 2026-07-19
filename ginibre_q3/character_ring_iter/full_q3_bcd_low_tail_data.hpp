#ifndef GINIBRE_Q3_FULL_Q3_BCD_LOW_TAIL_DATA_HPP
#define GINIBRE_Q3_FULL_Q3_BCD_LOW_TAIL_DATA_HPP

#include <array>

namespace full_q3_bcd_low_tail {

inline constexpr int b_first_rank = 18;
inline constexpr std::array<int, 4> b_onsets{{47, 45, 45, 45}};

inline constexpr int d_first_rank = 31;
inline constexpr std::array<int, 40> d_onsets{{
    49, 51, 51, 51, 51, 53, 53, 55, 53, 55,
    55, 57, 55, 57, 57, 59, 59, 61, 61, 61,
    61, 63, 61, 63, 63, 65, 63, 65, 65, 67,
    67, 69, 67, 69, 69, 71, 71, 71, 71, 73,
}};

}  // namespace full_q3_bcd_low_tail

#endif
