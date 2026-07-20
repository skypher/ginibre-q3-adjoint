#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <gmpxx.h>
#include <omp.h>
#include "frontier_runtime_scope.hpp"

#define main entry
namespace f09 {
#include "post_m29_bc_ninth_frontier_gmp.cpp"
} namespace f10 {
#include "post_m29_bc_tenth_frontier_gmp.cpp"
} namespace f11 {
#include "post_m29_bc_eleventh_frontier_gmp.cpp"
} namespace f12 {
#include "post_m29_bc_twelfth_frontier_gmp.cpp"
} namespace f13 {
#include "post_m29_bc_thirteenth_frontier_gmp.cpp"
} namespace f14 {
#include "post_m29_bc_fourteenth_frontier_gmp.cpp"
} namespace f15 {
#include "post_m29_bc_fifteenth_frontier_gmp.cpp"
} namespace f16 {
#include "post_m29_bc_sixteenth_frontier_gmp.cpp"
} namespace f17 {
#include "post_m29_bc_seventeenth_frontier_gmp.cpp"
} namespace f18 {
#include "post_m29_bc_eighteenth_frontier_gmp.cpp"
} namespace f19 {
#include "post_m29_bc_nineteenth_frontier_gmp.cpp"
} namespace f20 {
#include "post_m29_bc_twentieth_frontier_gmp.cpp"
} namespace f21 {
#include "post_m29_bc_twentyfirst_frontier_gmp.cpp"
} namespace f22 {
#include "post_m29_bc_twentysecond_frontier_gmp.cpp"
} namespace f23 {
#include "post_m29_bc_twentythird_frontier_gmp.cpp"
} namespace f24 {
#include "post_m29_bc_twentyfourth_frontier_gmp.cpp"
} namespace f25 {
#include "post_m29_bc_twentyfifth_frontier_gmp.cpp"
} namespace f26 {
#include "post_m29_bc_twentysixth_frontier_gmp.cpp"
} namespace f27 {
#include "post_m29_bc_twentyseventh_frontier_gmp.cpp"
}
#undef main

int main(int argc, char** argv) {
    if (argc != 1) return 2;
    const std::string name = argv[0];
    if (name.ends_with("ninth_frontier_gmp")) return f09::entry();
    if (name.ends_with("tenth_frontier_gmp")) return f10::entry();
    if (name.ends_with("eleventh_frontier_gmp")) return f11::entry();
    if (name.ends_with("twelfth_frontier_gmp")) return f12::entry();
    if (name.ends_with("thirteenth_frontier_gmp")) return f13::entry();
    if (name.ends_with("fourteenth_frontier_gmp")) return f14::entry();
    if (name.ends_with("fifteenth_frontier_gmp")) return f15::entry();
    if (name.ends_with("sixteenth_frontier_gmp")) return f16::entry();
    if (name.ends_with("seventeenth_frontier_gmp")) return f17::entry();
    if (name.ends_with("eighteenth_frontier_gmp")) return f18::entry();
    if (name.ends_with("nineteenth_frontier_gmp")) return f19::entry();
    if (name.ends_with("twentieth_frontier_gmp")) return f20::entry();
    if (name.ends_with("twentyfirst_frontier_gmp")) return f21::entry();
    if (name.ends_with("twentysecond_frontier_gmp")) return f22::entry();
    if (name.ends_with("twentythird_frontier_gmp")) return f23::entry();
    if (name.ends_with("twentyfourth_frontier_gmp")) return f24::entry();
    if (name.ends_with("twentyfifth_frontier_gmp")) return f25::entry();
    if (name.ends_with("twentysixth_frontier_gmp")) return f26::entry();
    if (name.ends_with("twentyseventh_frontier_gmp")) return f27::entry();
    return 2;
}
