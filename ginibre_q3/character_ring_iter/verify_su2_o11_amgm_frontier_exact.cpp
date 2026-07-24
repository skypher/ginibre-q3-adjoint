#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <atomic>
#include <initializer_list>
#include <mutex>
#include <thread>
#include <array>
#include <cctype>
#include <future>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;
using Rat = cpp_rational;

struct Interval {
    Rat lo;
    Rat hi;

    Interval() : lo(0), hi(0) {}
    explicit Interval(const Rat& value) : lo(value), hi(value) {}
    Interval(const Rat& lower, const Rat& upper) : lo(lower), hi(upper) {
        if (hi < lo) {
            throw std::runtime_error("invalid interval");
        }
    }
};

Interval operator+(const Interval& left, const Interval& right) {
    return Interval(left.lo + right.lo, left.hi + right.hi);
}

Interval operator-(const Interval& value) {
    return Interval(-value.hi, -value.lo);
}

Interval operator-(const Interval& left, const Interval& right) {
    return left + (-right);
}

Interval operator*(const Interval& left, const Interval& right) {
    const std::array<Rat, 4> values{
        left.lo * right.lo,
        left.lo * right.hi,
        left.hi * right.lo,
        left.hi * right.hi,
    };
    return Interval(
        *std::min_element(values.begin(), values.end()),
        *std::max_element(values.begin(), values.end())
    );
}

Interval operator*(const Interval& value, const Rat& scalar) {
    if (scalar >= 0) {
        return Interval(value.lo * scalar, value.hi * scalar);
    }
    return Interval(value.hi * scalar, value.lo * scalar);
}

Interval operator*(const Rat& scalar, const Interval& value) {
    return value * scalar;
}

Interval operator/(const Interval& value, const Rat& scalar) {
    if (scalar == 0) {
        throw std::runtime_error("division by zero");
    }
    return value * (Rat(1) / scalar);
}

Interval square(const Interval& value) {
    if (value.lo <= 0 && value.hi >= 0) {
        return Interval(0, std::max(value.lo * value.lo, value.hi * value.hi));
    }
    const Rat first = value.lo * value.lo;
    const Rat second = value.hi * value.hi;
    return Interval(std::min(first, second), std::max(first, second));
}

Interval power(Interval base, unsigned exponent) {
    Interval result(Rat(1));
    while (exponent != 0U) {
        if ((exponent & 1U) != 0U) {
            result = result * base;
        }
        exponent >>= 1U;
        if (exponent != 0U) {
            base = base * base;
        }
    }
    return result;
}

Rat decimal_rational(const std::string& text) {
    bool negative = false;
    std::size_t position = 0;
    if (!text.empty() && text.front() == '-') {
        negative = true;
        position = 1;
    }
    cpp_int numerator = 0;
    cpp_int denominator = 1;
    bool after_decimal = false;
    for (; position < text.size(); ++position) {
        const char character = text[position];
        if (character == '.') {
            if (after_decimal) {
                throw std::runtime_error("invalid decimal");
            }
            after_decimal = true;
            continue;
        }
        if (!std::isdigit(static_cast<unsigned char>(character))) {
            throw std::runtime_error("invalid decimal character");
        }
        numerator *= 10;
        numerator += static_cast<unsigned>(character - '0');
        if (after_decimal) {
            denominator *= 10;
        }
    }
    if (negative) {
        numerator = -numerator;
    }
    return Rat(numerator, denominator);
}

using Polynomial = std::vector<Rat>;

void trim(Polynomial& polynomial) {
    while (polynomial.size() > 1U && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

Polynomial derivative(const Polynomial& polynomial) {
    if (polynomial.size() <= 1U) {
        return Polynomial{Rat(0)};
    }
    Polynomial result(polynomial.size() - 1U);
    for (std::size_t index = 1; index < polynomial.size(); ++index) {
        result[index - 1U] = Rat(index) * polynomial[index];
    }
    trim(result);
    return result;
}

std::pair<Polynomial, Polynomial> divide_with_remainder(
    Polynomial dividend,
    Polynomial divisor
) {
    trim(dividend);
    trim(divisor);
    if (divisor.size() == 1U && divisor.front() == 0) {
        throw std::runtime_error("polynomial division by zero");
    }
    Polynomial quotient(
        dividend.size() >= divisor.size()
            ? dividend.size() - divisor.size() + 1U
            : 1U,
        Rat(0)
    );
    while (!(dividend.size() == 1U && dividend.front() == 0)
           && dividend.size() >= divisor.size()) {
        const std::size_t degree = dividend.size() - divisor.size();
        const Rat coefficient = dividend.back() / divisor.back();
        quotient[degree] = coefficient;
        for (std::size_t index = 0; index < divisor.size(); ++index) {
            dividend[degree + index] -= coefficient * divisor[index];
        }
        trim(dividend);
    }
    trim(quotient);
    trim(dividend);
    return {quotient, dividend};
}

std::vector<Polynomial> sturm_sequence(const Polynomial& polynomial) {
    std::vector<Polynomial> result{polynomial, derivative(polynomial)};
    while (!(result.back().size() == 1U && result.back().front() == 0)) {
        auto division = divide_with_remainder(
            result[result.size() - 2U],
            result.back()
        );
        Polynomial remainder = std::move(division.second);
        for (Rat& coefficient : remainder) {
            coefficient = -coefficient;
        }
        trim(remainder);
        if (remainder.size() == 1U && remainder.front() == 0) {
            break;
        }
        result.push_back(std::move(remainder));
    }
    return result;
}

Rat evaluate(const Polynomial& polynomial, const Rat& point) {
    Rat result = 0;
    for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend(); ++iterator) {
        result *= point;
        result += *iterator;
    }
    return result;
}

int sign_variations(const std::vector<Polynomial>& sequence, const Rat& point) {
    int previous = 0;
    int variations = 0;
    for (const Polynomial& polynomial : sequence) {
        const Rat value = evaluate(polynomial, point);
        const int sign = value > 0 ? 1 : (value < 0 ? -1 : 0);
        if (sign == 0) {
            continue;
        }
        if (previous != 0 && previous != sign) {
            ++variations;
        }
        previous = sign;
    }
    return variations;
}

Interval polynomial_interval(
    const std::vector<long long>& coefficients,
    const Interval& value
) {
    Interval result(Rat(0));
    for (auto iterator = coefficients.rbegin(); iterator != coefficients.rend(); ++iterator) {
        result = result * value + Interval(Rat(*iterator));
    }
    return result;
}

std::array<Interval, 6> isolated_roots() {
    const Polynomial polynomial{
        Rat(1), Rat(-2), Rat(-7), Rat(6), Rat(5), Rat(-5), Rat(1)
    };
    const std::vector<Polynomial> sturm = sturm_sequence(polynomial);
    const std::array<std::pair<std::string, std::string>, 6> bounds{{
        {
            "2.7709120513064197918007510440301977572109968326950698036048076174535151140673",
            "2.7709120513064197918007510440301977572109968326950698036048076174535151140676"
        },
        {
            "2.1361294934623116050236151182550332490669851049070363389593250031351633580822",
            "2.13612949346231160502361511825503324906698510490703633895932500313516335808226"
        },
        {
            "1.2410733605106461066981353749050871645473623184551428095938976480335891600339",
            "1.2410733605106461066981353749050871645473623184551428095938976480335891600343"
        },
        {
            "0.2907902259149287480607242147999630513672891357724104931565212626479200447",
            "0.29079022591492874806072421479996305136728913577241049315652126264792000452"
        },
        {
            "-0.4970214963422021972691194027027676929031803516522681392624782946634147091",
            "-0.4970214963422021972691194027027676929031803516522681392624782946634147086"
        },
        {
            "-0.9418836348521040543139645525875784544997302114780071771752890529540834175203",
            "-0.94188363485210405431396455258757845444997302114780071771752890529540834175198"
        }
    }};
    std::array<Interval, 6> roots{};
    for (std::size_t index = 0; index < roots.size(); ++index) {
        const Rat lower = decimal_rational(bounds[index].first);
        const Rat upper = decimal_rational(bounds[index].second);
        if (evaluate(polynomial, lower) * evaluate(polynomial, upper) >= 0) {
            throw std::runtime_error("root interval does not bracket a root");
        }
        const int count = sign_variations(sturm, lower) - sign_variations(sturm, upper);
        if (count != 1) {
            throw std::runtime_error("root interval is not isolating");
        }
        roots[index] = Interval(lower, upper);
    }
    const int total = sign_variations(sturm, Rat(-2)) - sign_variations(sturm, Rat(3));
    if (total != 6) {
        throw std::runtime_error("root isolation is incomplete");
    }
    return roots;
}

struct Term {
    std::string pair;
    Interval coefficient;
    std::array<Interval, 5> bases;
};

struct Chamber {
    int support;
    int parity;
    std::array<int, 5> signs;
    std::array<unsigned, 5> minimum_power;
    std::vector<Term> positive;
    std::vector<Term> negative;
};

std::array<Interval, 5> orbit_characters(const Interval& x) {
    return {
       x,
        polynomial_interval({-1, -1, 1}, x),
        polynomial_interval({1, -1, -2, 1}, x),
        polynomial_interval({0, 3, 0, -3, 1}, x),
        polynomial_interval({-1, -2, 5, 2, -4, 1}, x),
    };
}

Chamber make_chamber(int support, int parity) {
    Chamber result{};
    result.support = support;
    result.parity = parity;
    int code = support;
    std::vector<int> active;
    for (int label = 0; label < 5; ++label) {
        const int digit = code % 3;
        code /= 3;
        result.signs[static_cast<std::size_t>(label)]
            = digit == 0 ? 0 : (digit == 1 ? 1 : -1);
        if (digit != 0) {
            active.push_back(label);
        }
    }
    for (std::size_t index = 0; index < active.size(); ++index) {
        result.minimum_power[static_cast<std::size_t>(active[index])]
            = ((parity >> index) & 1) != 0 ? 1U : 2U;
    }

    const std::array<Interval, 6> roots = isolated_roots();
    std::array<std::array<Interval, 5>, 6> values{};
    std::array<Interval, 6> trace_weights{};
    for (std::size_t node = 0; node < roots.size(); ++node) {
        values[node] = orbit_characters(roots[node]);
        trace_weights[node] = (Interval(Rat(3)) - roots[node]) / Rat(13);
    }
    for (int left = 0; left < 6; +*ìeft) {
        for (int right = left + 1; right < 6; ++right) {
            Interval coefficient = Rat(2)
                * trace_weights[static_cast<std::size_t>(left)]
                * trace_weights[static_cast<std::size_t>(right)];
            std::array<Interval, 5> bases{
                Interval(Rat(1)), Interval(Rat(1)), Interval(Rat(1)),
                Interval(Rat(1)), Interval(Rat(1))
            };
            bool negative = false;
            for (const int label : active) {
                const Interval factor
                    = values[static_cast<std::size_t>(left)][static_cast<std::size_t>(label)]
                    + Rat(result.signs[static_cast<std::size_t>(label)])
                        * values[static_cast<std::size_t>(right)][static_cast<std::size_t>(label)];
                if (!(factor.lo > 0 || factor.hi < 0)) {
                    throw std::runtime_error("unresolved factor sign");
                }
                const unsigned exponent
                    = result.minimum_power[static_cast<std::size_t>(label)];
                if (factor.hi < 0 && (exponent & 1U) != 0U) {
                    negative = !negative;
                }
                coefficient = coefficient * power(
                    factor.hi < 0 ? -factor : factor,
                    exponent
                );
                bases[static_cast<std::size_t>(label)] = square(factor);
            }
            const std::string pair = std::to_string(left) + std::to_string(right);
            Term term{pair, coefficient, bases};
            (negative ? result.negative : result.positive).push_back(std::move(term));
        }
    }
    return result;
}

struct Allocation {
    std::string negative;
    std::vector<std::pair<std::string, unsigned>> positive_weights;
};

struct Certificate {
    const char* name;
    int support;
    int parity;
    std::array<unsigned, 5> floors;
    std::vector<int> variables;
    std::vector<Allocation> allocations;
};

Allocation a(const char* negative, std::initializer_list<std::pair<const char*, unsigned>> entries) {
    Allocation result;
    result.negative = negative;
    for (const auto& entry : entries) {
        result.positive_weights.emplace_back(entry.first, entry.second);
    }
    return result;
}

const std::vector<Certificate> certificates{
    {"s86_p4_r7",86,4,{1,1,0,0,1},{0,1,4},{
        a("12",{{"04",100}}),
        a("13",{{"02",69},{"03",1},{"04",25},{"05",1},{"23",3},{"24",1}}),
        a("14",{{"02",11},{"03",17},{"04",59},{"24",13}}),
        a("15",{{"02",1},{"04",41},{"05",58}}),
        a("34",{{"03",99},{"23",1}}),
        a("35",{{"05",89},{"23",11}}),
    }},
    {"s86_p6_r7",86,6,{1,1,0,0,1},{0,1,4},{
        a("12",{{"01",53},{"04",47}}),
        a("13",{{"02",70},{"04",15},{"05",10},{"25",1},{"34",1},{"35",3}}),
        a("14",{{"02",20},{"03",5},{"04",48},{"05",27}}),
        a("15",{{"04",48},{"05",52}}),
        a("23",{{"03",3},{"05",62},{"25",3},{"34",6},{"35",26}}),
        a("24",{{"02",82},{"03",5},{"04",1},{"05",12}}),
    }},
    {"s89_p4_r5",89,4,{1,0,0,0,1},{0,4},{
        a("12",{{"04",100}}),a("13",{{"02",36},{"04",64}}),
        a("14",{{"04",47},{"25",53}}),a("15",{{"04",91},{"24",9}}),
        a("34",{{"05",100}}), a("35",{{"04",100}}),
    }},
    {"s89_p4_r7",89,4,{1,1,0,0,1},{0,1,4},{
        a("12",{{"04",100}}),a("13",{{"02",36},{"04",64}}),
        a("14",{{"03",28},{"04",38},{"25",34}}),
        a("15",{{"03",21},{"04",79}}),
        a("34",{{"04",100}}), a("35",{{"04",100}}),
    }},
    {"s89_p7_r5_strip",89,7,{1,0,0,0,1},{0},{
        a("12",{{"23",100}}),a("13",{{"04",100}}),
        a("14",{{"02",14},{"03",10},{"04",63},{"05",8},{"35",5}}),
        a("15",{{"02",21},{"05",79}}),a("24",{{"03",1},{"04",99}}),
        a("25",{{"02",19},{"03",1},{"04",45},{"05",9},{"35",26}}),
        a("45",{{"34",100}}),
    }},
    {"s89_p7_r5_tail",89,7,{1,0,0,0,2},{0,4},{
        a("12",{{"04",100}}), a("13",{{"02",36},{"04",64}}),
        a("14",{{"04",78},{"35",22}}),a("15",{{"04",94},{"35",6}}),
        a("24",{{"04",59},{"05",41}}), a("25",{{"04",100}}),a("45",{{"04",100}}),
    }},
    {"s89_p7_r7",89,7,{1,1,0,0,1},{0,1,4},{
        a("12",{{"04",55},{"34",45}}),a("13",{{"02",36},{"04",64}}),
        a("14",{{"03",76},{"04",22},{"05",2}}),a("15",{{"03",20},{"04",80}}),
        a("24",{{"04",72},{"35",28}}),a("25",{{"04",67},{"35",33}}),
        a("45",{{"03",19},{"04",81}}),
    }},
    {"s95_p8_r11",95,8,{1,1,0,0,1},{0,1,4},{
        a("12",{{"01",58},{"02",42}}),a("13",{{"02",75},{"04",25}}),
        a("14",{{"02",38},{"03",1},{"05",61}}),a("15",{{"04",60},{"05",40}}),
        a("34",{{"02",100}}),a("35",{{"02",100}}),
    }},
    {"s;h¢V(º¯{®wßŒŠmyÛh‘éì•çí