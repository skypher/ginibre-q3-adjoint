#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
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
Interval operator-(const Interval& value) { return Interval(-value.hi, -value.lo); }
Interval operator-(const Interval& left, const Interval& right) { return left + (-right); }
Interval operator*(const Interval& left, const Interval& right) {
    const std::array<Rat, 4> values{
        left.lo * right.lo, left.lo * right.hi,
        left.hi * right.lo, left.hi * right.hi,
    };
    return Interval(
        *std::min_element(values.begin(), values.end()),
        *std::max_element(values.begin(), values.end())
    );
}
Interval operator*(const Interval& value, const Rat& scalar) {
    if (scalar >= 0) return Interval(value.lo * scalar, value.hi * scalar);
    return Interval(value.hi * scalar, value.lo * scalar);
}
Interval operator*(const Rat& scalar, const Interval& value) { return value * scalar; }
Interval operator/(const Interval& value, const Rat& scalar) {
    if (scalar == 0) throw std::runtime_error("division by zero");
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
        if ((exponent & 1U) != 0U) result = result * base;
        exponent >>= 1U;
        if (exponent != 0U) base = base * base;
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
            if (after_decimal) throw std::runtime_error("invalid decimal");
            after_decimal = true;
            continue;
        }
        if (!std::isdigit(static_cast<unsigned char>(character))) {
            throw std::runtime_error("invalid decimal character");
        }
        numerator *= 10;
        numerator += static_cast<unsigned>(character - '0');
        if (after_decimal) denominator *= 10;
    }
    if (negative) numerator = -numerator;
    return Rat(numerator, denominator);
}

using Polynomial = std::vector<Rat>;
void trim(Polynomial& polynomial) {
    while (polynomial.size() > 1U && polynomial.back() == 0) polynomial.pop_back();
}
Polynomial derivative(const Polynomial& polynomial) {
    if (polynomial.size() <= 1U) return Polynomial{Rat(0)};
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
            ? dividend.size() - divisor.size() + 1U : 1U,
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
        auto division = divide_with_remainder(result[result.size() - 2U], result.back());
        Polynomial remainder = std::move(division.second);
        for (Rat& coefficient : remainder) coefficient = -coefficient;
        trim(remainder);
        if (remainder.size() == 1U && remainder.front() == 0) break;
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
        if (sign == 0) continue;
        if (previous != 0 && previous != sign) ++variations;
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
        {"2.7709120513064197918007510440301977572109968326950698036048076174535151140673",
         "2.7709120513064197918007510440301977572109968326950698036048076174535151140676"},
        {"2.1361294934623116050236151182550332490669851049070363389593250031351633580822",
         "2.1361294934623116050236151182550332490669851049070363389593250031351633580826"},
        {"1.2410733605106461066981353749050871645473623184551428095938976480335891600339",
         "1.2410733605106461066981353749050871645473623184551428095938976480335891600343"},
        {"0.2907902259149287480607242147999630513672891357724104931565212626264792000447",
         "0.2907902259149287480607242147999630513672891357724104931565212626264792000452"},
        {"-0.4970214963422021972692611994027027676929031803516522681392624782946634147091",
         "-0.4970214963422021972692611994027027676929031803516522681392624782946634147086"},
        {"-0.9418836348521040543139645525875784544997302114780071771752890529540834175203",
         "-0.9418836348521040543139645525875784544997302114780071771752890529540834175198"}
    }};
    std::array<Interval, 6> roots{};
    for (std::size_t index = 0; index < roots.size(); ++index) {
        const Rat lower = decimal_rational(bounds[index].first);
        const Rat upper = decimal_rational(bounds[index].second);
        if (evaluate(polynomial, lower) * evaluate(polynomial, upper) >= 0) {
            throw std::runtime_error("root interval does not bracket a root");
        }
        const int count = sign_variations(sturm, lower) - sign_variations(sturm, upper);
        if (count != 1) throw std::runtime_error("root interval is not isolating");
        roots[index] = Interval(lower, upper);
    }
    const int total = sign_variations(sturm, Rat(-2)) - sign_variations(sturm, Rat(3));
    if (total != 6) throw std::runtime_error("root isolation is incomplete");
    return roots;
}

struct SpectralTerm {
    std::string pair;
    Interval coefficient;
    Interval a_base;
    Interval c_base;
};
struct SpectralData {
    std::vector<SpectralTerm> positive;
    std::vector<SpectralTerm> negative;
};
SpectralData build_spectral_data() {
    const std::array<Interval, 6> roots = isolated_roots();
    const std::vector<long long> b5_coefficients{-1, -2, 5, 2, -4, 1};
    std::array<Interval, 6> b5{};
    std::array<Interval, 6> weights{};
    for (std::size_t index = 0; index < roots.size(); ++index) {
        b5[index] = polynomial_interval(b5_coefficients, roots[index]);
        weights[index] = (Interval(Rat(3)) - roots[index]) / Rat(13);
        if (weights[index].lo <= 0) throw std::runtime_error("nonpositive trace weight");
    }
    const std::vector<std::string> positive_order{
        "01", "45", "23", "03", "02", "05", "04", "24", "25"
    };
    const std::vector<std::string> negative_order{
        "15", "13", "35", "14", "34", "12"
    };
    std::map<std::string, SpectralTerm> terms;
    for (int left = 0; left < 6; ++left) {
        for (int right = left + 1; right < 6; ++right) {
            const Interval difference = roots[static_cast<std::size_t>(left)]
                - roots[static_cast<std::size_t>(right)];
            const Interval sum = b5[static_cast<std::size_t>(left)]
                + b5[static_cast<std::size_t>(right)];
            const Interval a_base = square(difference);
            const Interval c_base = square(sum);
            const Interval coefficient = Rat(2)
                * weights[static_cast<std::size_t>(left)]
                * weights[static_cast<std::size_t>(right)]
                * a_base * sum;
            const std::string pair = std::to_string(left) + std::to_string(right);
            terms.emplace(pair, SpectralTerm{pair, coefficient, a_base, c_base});
        }
    }
    SpectralData result;
    for (const std::string& pair : positive_order) {
        const SpectralTerm term = terms.at(pair);
        if (term.coefficient.lo <= 0) throw std::runtime_error("expected positive spectral coefficient");
        result.positive.push_back(term);
    }
    for (const std::string& pair : negative_order) {
        SpectralTerm term = terms.at(pair);
        if (term.coefficient.hi >= 0) throw std::runtime_error("expected negative spectral coefficient");
        term.coefficient = -term.coefficient;
        result.negative.push_back(term);
    }
    return result;
}

struct RegionCertificate {
    const char* name;
    unsigned p_floor;
    unsigned q_floor;
    bool check_c_base;
    std::array<std::array<unsigned, 9>, 6> weights;
};
constexpr std::array<std::array<unsigned, 9>, 6> interior_weights{{
    {{0,0,0,0,5,8,83,3,1}}, {{0,0,0,0,50,0,32,18,0}},
    {{0,0,0,0,42,0,42,0,16}}, {{0,0,0,1,0,26,12,0,61}},
    {{0,58,0,0,0,0,0,42,0}}, {{0,0,0,100,0,0,0,0,0}},
}};
constexpr std::array<std::array<unsigned, 9>, 6> q3_weights{{
    {{0,0,0,0,0,11,77,11,1}}, {{0,0,0,0,38,0,32,30,0}},
    {{0,0,0,0,1,21,0,11,67}}, {{0,0,0,0,1,7,42,0,50}},
    {{0,0,0,0,100,0,0,0,0}}, {{0,0,0,0,0,100,0,0,0}},
}};
constexpr std::array<std::array<unsigned, 9>, 6> q2_weights{{
    {{0,0,0,0,0,19,68,13,0}}, {{0,0,0,0,15,0,58,20,7}},
    {{0,0,0,2,45,3,22,10,18}}, {{0,0,0,2,0,78,0,20,0}},
    {{0,19,0,0,81,0,0,0,0}}, {{0,0,0,0,0,0,0,0,100}},
}};
constexpr std::array<std::array<unsigned, 9>, 6> q1_weights{{
    {{0,0,0,0,0,23,63,14,0}}, {{0,0,0,0,12,30,0,32,26}},
    {{0,0,0,18,3,6,2,11,60}}, {{0,0,0,0,2,16,30,3,49}},
    {{0,18,61,0,8,0,13,0,0}}, {{0,8,3,0,0,23,0,66,0}},
}};
constexpr std::array<std::array<unsigned, 9>, 6> q0_weights{{
    {{0,0,0,1,0,47,39,4,9}}, {{0,0,0,3,0,86,0,0,11}},
    {{0,0,2,11,10,3,63,9,2}}, {{0,0,0,3,1,59,1,9,27}},
    {{0,3,0,0,0,0,38,59,0}}, {{0,7,0,0,0,0,0,93,0}},
}};
const std::array<RegionCertificate, 5> certificates{{
    {"interior_q_ge_4",0U,4U,true,interior_weights},
    {"face_q_3",0U,3U,false,q3_weights},
    {"face_q_2_p_ge_1",1U,2U,false,q2_weights},
    {"face_q_1_p_ge_1",1U,1U,false,q1_weights},
    {"face_q_0_p_ge_2",2U,0U,false,q0_weights},
}};
Interval shifted_magnitude(const SpectralTerm& term,unsigned p_floor,unsigned q_floor) {
    return term.coefficient * power(term.a_base,p_floor) * power(term.c_base,q_floor);
}
void check_region(const SpectralData& data,const RegionCertificate& certificate) {
    constexpr unsigned denominator=100U;
    for (std::size_t negative=0;negative<data.negative.size();++negative) {
        unsigned total=0U;
        Interval a_product(Rat(1));
        Interval c_product(Rat(1));
        for (std::size_t positive=0;positive<data.positive.size();++positive) {
            const unsigned weight=certificate.weights[negative][positive];
            total+=weight;
            if (weight!=0U) {
                a_product=a_product*power(data.positive[positive].a_base,weight);
                if (certificate.check_c_base) {
                    c_product=c_product*power(data.positive[positive].c_base,weight);
                }
            }
        }
        if (total!=denominator) throw std::runtime_error(std::string(certificate.name)+": invalid AM-GM weights");
        if (a_product.lo<power(data.negative[negative].a_base,denominator).hi) {
            throw std::runtime_error(std::string(certificate.name)+": A-base inequality failed");
        }
        if (certificate.check_c_base
            && c_product.lo<power(data.negative[negative].c_base,denominator).hi) {
            throw std::runtime_error(std::string(certificate.name)+": C-base inequality failed");
        }
    }
    for (std::size_t positive=0;positive<data.positive.size();++positive) {
        Interval used(Rat(0));
        for (std::size_t negative=0;negative<data.negative.size();++negative) {
            const unsigned weight=certificate.weights[negative][positive];
            if (weight==0U) continue;
            used=used+shifted_magnitude(data.negative[negative],certificate.p_floor,certificate.q_floor)
                *Rat(weight)/Rat(denominator);
        }
        const Interval capacity=shifted_magnitude(data.positive[positive],certificate.p_floor,certificate.q_floor);
        if (used.hi>capacity.lo) throw std::runtime_error(std::string(certificate.name)+": capacity inequality failed");
    }
}

constexpr int orbit_rank=6;
constexpr int tensor_rank=orbit_rank*orbit_rank;
using State=std::array<cpp_int,tensor_rank>;
std::vector<int> fusion_outputs(int left,int right) {
    const int lower=std::abs(left-right);
    const int upper=std::min(left+right,11-left-right);
    std::vector<int> result;
    for (int output=lower;output<=upper;++output) result.push_back(output);
    return result;
}
State fusion_step(const State& state,int label,int sign) {
    State result{};
    for (int left=0;left<orbit_rank;++left) {
        for (int right=0;right<orbit_rank;++right) {
            const cpp_int coefficient=state[static_cast<std::size_t>(left*orbit_rank+right)];
            if (coefficient==0) continue;
            for (const int output:fusion_outputs(left,label)) {
                result[static_cast<std::size_t>(output*orbit_rank+right)]+=coefficient;
            }
            for (const int output:fusion_outputs(right,label)) {
                result[static_cast<std::size_t>(left*orbit_rank+output)]+=sign*coefficient;
            }
        }
    }
    return result;
}
State paired_fusion_step(const State& state,int label,int sign) {
    return fusion_step(fusion_step(state,label,sign),label,sign);
}
cpp_int exact_corner(unsigned p,unsigned q) {
    State state{};
    state.front()=1;
    state=paired_fusion_step(state,1,-1);
    state=fusion_step(state,5,1);
    for (unsigned index=0;index<p;++index) state=paired_fusion_step(state,1,-1);
    for (unsigned index=0;index<q;++index) state=paired_fusion_step(state,5,1);
    return state.front();
}

}  // namespace

int main() {
    try {
        const SpectralData data=build_spectral_data();
        std::vector<std::future<void>> workers;
        workers.reserve(certificates.size());
        for (const RegionCertificate& certificate:certificates) {
            workers.push_back(std::async(std::launch::async,[&data,&certificate]() {
                check_region(data,certificate);
            }));
        }
        for (std::future<void>& worker:workers) worker.get();
        const std::array<std::pair<unsigned,unsigned>,4> leaves{{
            {0U,0U},{1U,0U},{0U,1U},{0U,2U}
        }};
        for (const auto& leaf:leaves) {
            if (exact_corner(leaf.first,leaf.second)!=0) {
                throw std::runtime_error("unexpected boundary leaf value");
            }
        }
        std::cout<<"SU2_O11_FIRST_CHAMBER_EXACT PASS"
                 <<" roots=6 spectral_pairs=15 regions=5"
                 <<" amgm_denominator=100 exact_leaves=4 threads=5\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"error: "<<error.what()<<'\n';
        return 1;
    }
}
