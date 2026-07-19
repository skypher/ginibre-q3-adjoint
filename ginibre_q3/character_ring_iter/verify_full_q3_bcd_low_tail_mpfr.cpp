#define main post_m29_bc_layered_mgf_legacy_main
#include "post_m29_bc_layered_mgf_mpfr.cpp"
#undef main

#include "full_q3_bcd_remaining_data.hpp"

#include <set>

namespace {

struct FullTailResult {
    int onset;
    Real margin_lower;
    const char* moment_bound;
};

Real exact_full_q3_event_probability(
    const ParameterSet& parameters,
    int rank,
    const std::vector<mpz_class>& factorials
) {
    const Interval rank_interval = exact_ui(static_cast<unsigned long>(rank));
    const Interval c = mul(parameter(parameters.r), rank_interval);
    const Interval q = mul(
        parameter(parameters.square_alpha), sqrt_positive(rank_interval)
    );
    const Interval central_cutoff = parameter(parameters.character_cutoff);
    const Interval lower_gap = parameter(parameters.lower_gap);
    const Interval upper_start = parameter(parameters.upper_start);
    const Interval upper_step = parameter(parameters.upper_step);
    const Interval lower_fraction = parameter(parameters.lower_lambda_fraction);
    const Interval upper_fraction = parameter(parameters.upper_lambda_fraction);
    const Interval square_fraction = parameter(parameters.square_lambda_fraction);
    const Interval threshold = sqrt_positive(add(
        add(scale_ui(c, 2), scale_ui(central_cutoff, 2)), q
    ));
    const Interval saddle = add(threshold, lower_gap);

    Real minimum_pivot;
    mpfr_set_ui(minimum_pivot.value, 1, MPFR_RNDD);
    Interval defining_tail;
    if (parameters.family == 'B') {
        const ExactLayeredTail positive = exact_so_odd_layered_tail(
            rank,
            1,
            saddle,
            lower_gap,
            upper_start,
            upper_step,
            lower_fraction,
            upper_fraction,
            parameters.search_from_first_hit,
            factorials,
            minimum_pivot
        );
        defining_tail = positive.probability;
        const Interval negative_support = exact_ui(
            static_cast<unsigned long>(2 * rank - 1)
        );
        if (!certainly_positive(sub(threshold, negative_support))) {
            const ExactLayeredTail negative = exact_so_odd_layered_tail(
                rank,
                -1,
                saddle,
                lower_gap,
                upper_start,
                upper_step,
                lower_fraction,
                upper_fraction,
                parameters.search_from_first_hit,
                factorials,
                minimum_pivot
            );
            defining_tail = add(defining_tail, negative.probability);
        }
    } else if (parameters.family == 'C') {
        defining_tail = scale_ui(
            exact_o_even_minus_layered_tail(
                2 * rank + 2,
                saddle,
                lower_gap,
                upper_start,
                upper_step,
                lower_fraction,
                upper_fraction,
                parameters.search_from_first_hit,
                factorials,
                minimum_pivot
            ).probability,
            2
        );
    } else if (parameters.family == 'D') {
        defining_tail = scale_ui(
            exact_so_even_layered_tail(
                rank,
                saddle,
                lower_gap,
                upper_start,
                upper_step,
                factorials,
                minimum_pivot
            ).probability,
            2
        );
    } else {
        die("unknown exact full-Q3 event family");
    }

    const Interval q_minus_one = sub(q, exact_ui(1));
    const Interval square_lambda = mul(
        square_fraction,
        div_positive(q_minus_one, exact_ui(2))
    );
    Real dummy_fredholm;
    mpfr_set_ui(dummy_fredholm.value, 1, MPFR_RNDD);
    Real square_pivot;
    mpfr_set_ui(square_pivot.value, 1, MPFR_RNDD);
    const Interval square_tail = square_tail_upper(
        parameters.family,
        rank,
        q,
        square_lambda,
        factorials,
        dummy_fredholm,
        square_pivot,
        parameters.family != 'B'
    );
    const Interval probability = sub(defining_tail, square_tail);
    if (!certainly_positive(probability)) {
        die("exact full-Q3 event probability is not positive");
    }
    Real answer;
    mpfr_set(answer.value, probability.lo.value, MPFR_RNDD);
    return answer;
}

Interval exact_defining_log_mgf(
    char family,
    int rank,
    int sign,
    const Interval& argument,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    if (sign != 1 && sign != -1) die("invalid defining-trace MGF sign");
    if (family == 'B') {
        return log_so_odd_defining_mgf(
            rank, argument, sign, factorials, minimum_pivot
        );
    }
    if (family == 'C') {
        return log_o_even_minus_defining_mgf(
            2 * rank + 2, argument, factorials, minimum_pivot
        );
    }
    if (family == 'D') {
        return log_so_even_defining_mgf(
            rank, argument, factorials, minimum_pivot
        );
    }
    die("unknown defining-trace MGF family");
}

Interval exact_one_step_signed_tail(
    char family,
    int rank,
    int sign,
    const Interval& threshold,
    const Interval& h,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    const Interval saddle = add(threshold, h);
    const Interval upper_argument = add(saddle, h);
    const Interval log_mgf_s = exact_defining_log_mgf(
        family, rank, sign, saddle, factorials, minimum_pivot
    );
    const Interval log_mgf_lower = exact_defining_log_mgf(
        family, rank, sign, threshold, factorials, minimum_pivot
    );
    const Interval log_mgf_upper = exact_defining_log_mgf(
        family, rank, sign, upper_argument, factorials, minimum_pivot
    );
    const Interval lower_tilt_tail = exp_interval(add(
        mul(h, threshold), sub(log_mgf_lower, log_mgf_s)
    ));
    const Interval upper_endpoint = add(threshold, scale_ui(h, 2));
    const Interval support = exact_ui(static_cast<unsigned long>(
        family == 'B'
            ? (sign == 1 ? 2 * rank + 1 : 2 * rank - 1)
            : 2 * rank
    ));
    const bool reaches_support =
        mpfr_cmp(upper_endpoint.lo.value, support.hi.value) >= 0;
    const Interval upper_tilt_tail = reaches_support
        ? exact_ui(0)
        : exp_interval(add(
              neg(mul(h, upper_endpoint)), sub(log_mgf_upper, log_mgf_s)
          ));
    const Interval tilted_mass = one_minus(add(
        lower_tilt_tail, upper_tilt_tail
    ));
    if (!certainly_positive(tilted_mass)) {
        die("exact one-step tilted mass is not positive");
    }
    return mul(
        exp_interval(sub(
            log_mgf_s, mul(saddle, upper_endpoint)
        )),
        tilted_mass
    );
}

Real exact_one_step_full_q3_event_probability(
    char family,
    int rank,
    const RationalParameter& r_parameter,
    const RationalParameter& q_parameter,
    const RationalParameter& central_parameter,
    const RationalParameter& h_parameter,
    const std::vector<mpz_class>& factorials
) {
    const Interval rank_interval = exact_ui(static_cast<unsigned long>(rank));
    const Interval c = mul(parameter(r_parameter), rank_interval);
    const Interval q = parameter(q_parameter);
    const Interval central_cutoff = parameter(central_parameter);
    const Interval h = parameter(h_parameter);
    const Interval threshold = sqrt_positive(add(
        add(scale_ui(c, 2), scale_ui(central_cutoff, 2)), q
    ));
    Real minimum_pivot;
    mpfr_set_ui(minimum_pivot.value, 1, MPFR_RNDD);

    Interval defining_tail = exact_one_step_signed_tail(
        family, rank, 1, threshold, h, factorials, minimum_pivot
    );
    if (family == 'B') {
        const Interval negative_support = exact_ui(
            static_cast<unsigned long>(2 * rank - 1)
        );
        if (!certainly_positive(sub(threshold, negative_support))) {
            defining_tail = add(
                defining_tail,
                exact_one_step_signed_tail(
                    family, rank, -1, threshold, h,
                    factorials, minimum_pivot
                )
            );
        }
    } else {
        defining_tail = scale_ui(defining_tail, 2);
    }

    const Interval q_minus_one = sub(q, exact_ui(1));
    if (!certainly_positive(q_minus_one)) {
        die("one-step square-trace cutoff does not exceed one");
    }
    const Interval square_lambda = div_positive(q_minus_one, exact_ui(2));
    Real dummy_fredholm;
    mpfr_set_ui(dummy_fredholm.value, 1, MPFR_RNDD);
    Real square_pivot;
    mpfr_set_ui(square_pivot.value, 1, MPFR_RNDD);
    const Interval square_tail = square_tail_upper(
        family,
        rank,
        q,
        square_lambda,
        factorials,
        dummy_fredholm,
        square_pivot,
        family != 'B'
    );
    const Interval event_probability = sub(defining_tail, square_tail);
    if (!certainly_positive(event_probability)) {
        die("one-step full-Q3 event probability is not positive");
    }
    Real answer;
    mpfr_set(answer.value, event_probability.lo.value, MPFR_RNDD);
    return answer;
}

Interval point_from_lower_bound(const Real& lower_bound) {
    if (mpfr_sgn(lower_bound.value) <= 0) {
        die("nonpositive probability lower bound in full-Q3 tail");
    }
    Interval answer;
    mpfr_set(answer.lo.value, lower_bound.value, MPFR_RNDD);
    mpfr_set(answer.hi.value, lower_bound.value, MPFR_RNDU);
    return answer;
}

int first_odd_above(int value) {
    int answer = value + 1;
    if (answer % 2 == 0) ++answer;
    return answer;
}

Interval square_trace_log_mgf_upper(
    char family,
    int rank,
    const Interval& argument,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    if (family == 'B') return add(argument, square(argument));
    if (family != 'C' && family != 'D') {
        die("unknown square-trace MGF family");
    }
    const int m0 = family == 'D' ? (rank + 1) / 2 : (rank + 2) / 2;
    const int m1 = family == 'D' ? rank / 2 : (rank + 1) / 2;
    const Interval first = family == 'D'
        ? log_so_even_defining_mgf(m0, argument, factorials, minimum_pivot)
        : log_o_even_minus_defining_mgf(
              2 * m0, argument, factorials, minimum_pivot
          );
    const Interval second = log_so_odd_defining_mgf(
        m1, argument, -1, factorials, minimum_pivot
    );
    return add(argument, add(first, second));
}

Interval square_trace_log_moment_upper(
    char family,
    int rank,
    int degree,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    if (degree < 1) die("invalid square-trace moment degree");
    const Interval l = exact_ui(static_cast<unsigned long>(degree));
    const Interval v = div_positive(
        sub(
            sqrt_positive(add(exact_ui(1), scale_ui(l, 8))),
            exact_ui(1)
        ),
        exact_ui(4)
    );
    return add(
        square_trace_log_mgf_upper(
            family, rank, v, factorials, minimum_pivot
        ),
        mul(
            l,
            sub(
                sub(log_positive(l), exact_ui(1)),
                log_positive(v)
            )
        )
    );
}

FullTailResult full_tail_from_probability(
    char family,
    int rank,
    const RationalParameter& r_parameter,
    const RationalParameter& central_parameter,
    const Real& positive_probability_lower,
    const std::vector<mpz_class>& factorials
) {
    if (family != 'B' && family != 'C' && family != 'D') {
        die("unknown family in full-Q3 tail audit");
    }
    const int negative_endpoint =
        family == 'D' && rank % 2 != 0 ? rank - 2 : rank;
    const Interval rank_interval = exact_ui(static_cast<unsigned long>(rank));
    const Interval c = mul(parameter(r_parameter), rank_interval);
    const Interval twice_endpoint = exact_ui(
        static_cast<unsigned long>(2 * negative_endpoint)
    );
    if (!certainly_positive(sub(c, twice_endpoint))) {
        die("full-Q3 propagation base does not exceed twice the negative cap");
    }

    const Interval central_cutoff = parameter(central_parameter);
    const Interval central_probability = one_minus(
        reciprocal_positive(square(central_cutoff))
    );
    if (!certainly_positive(central_probability)) {
        die("nonpositive full-Q3 central probability");
    }
    const Interval base_probability = mul(
        point_from_lower_bound(positive_probability_lower),
        central_probability
    );

    const int stable_through = family == 'D' ? rank - 1 : 2 * rank + 1;
    const Interval support = exact_ui(static_cast<unsigned long>(
        family == 'B' ? 2 * rank + 1 : 2 * rank
    ));
    Real minimum_pivot;
    mpfr_set_ui(minimum_pivot.value, 1, MPFR_RNDD);

    const Interval log_probability = log_positive(base_probability);
    const Interval log_c = log_positive(c);
    const Interval log_support = log_positive(support);
    for (int degree = first_odd_above(stable_through);
         degree <= 10001;
         degree += 2) {
        const Interval degree_interval = exact_ui(
            static_cast<unsigned long>(degree)
        );
        const Interval support_log_moment = mul(degree_interval, log_support);
        const Interval mgf_log_moment = square_trace_log_moment_upper(
            family, rank, degree, factorials, minimum_pivot
        );
        const bool use_support =
            mpfr_cmp(support_log_moment.hi.value, mgf_log_moment.hi.value) <= 0;
        const Interval log_moment_upper =
            use_support ? support_log_moment : mgf_log_moment;
        const Interval margin = sub(
            add(log_probability, mul(degree_interval, log_c)),
            log_moment_upper
        );
        if (certainly_positive(margin)) {
            FullTailResult result;
            result.onset = degree;
            mpfr_set(result.margin_lower.value, margin.lo.value, MPFR_RNDD);
            result.moment_bound = use_support
                ? "support"
                : (family == 'B' ? "B_MGF" : "exact_square_MGF");
            return result;
        }
    }
    die("full-Q3 tail onset search exceeded 10001");
}

void print_tail(
    char family,
    int rank,
    const FullTailResult& tail,
    const Real& positive_probability
) {
    std::cout << "FULL_Q3_LOW_TAIL row=" << family << '_' << rank
              << " onset=" << tail.onset
              << " moment_bound=" << tail.moment_bound
              << " event_probability_lower="
              << format_real(positive_probability, 24)
              << " onset_log_margin_lower="
              << format_real(tail.margin_lower, 24) << '\n';
}

void record_directed_tail(
    char family,
    int rank,
    const ParameterSet& parameters,
    const FullTailResult& tail,
    const Real& positive_probability,
    std::set<std::pair<char, int>>& checked_rows
) {
    const full_q3_bcd_remaining::RowCutoff* cutoff =
        full_q3_bcd_remaining::find_row(family, rank);
    if (cutoff == nullptr
        || cutoff->tail_method
            != full_q3_bcd_remaining::TailMethod::directed_interval
        || cutoff->tail_onset != tail.onset) {
        die("full-Q3 directed-tail onset/row ledger mismatch");
    }
    if (!checked_rows.emplace(family, rank).second) {
        die("duplicate full-Q3 directed-tail row");
    }
    const auto print_rational = [](const RationalParameter& value) {
        return std::to_string(value.numerator) + '/'
            + std::to_string(value.denominator);
    };
    std::cout << "FULL_Q3_LOW_TAIL_SCHEDULE row=" << family << '_' << rank
              << " r=" << print_rational(parameters.r)
              << " square_alpha=" << print_rational(parameters.square_alpha)
              << " central=" << print_rational(parameters.character_cutoff)
              << " lower_gap=" << print_rational(parameters.lower_gap)
              << " upper_start=" << print_rational(parameters.upper_start)
              << " upper_step=" << print_rational(parameters.upper_step)
              << " lower_lambda_fraction="
              << print_rational(parameters.lower_lambda_fraction)
              << " upper_lambda_fraction="
              << print_rational(parameters.upper_lambda_fraction)
              << " square_lambda_fraction="
              << print_rational(parameters.square_lambda_fraction)
              << " threshold_anchored="
              << (parameters.search_from_first_hit ? 1 : 0) << '\n';
    print_tail(family, rank, tail, positive_probability);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout << "FULL_Q3_LOW_TAIL_CONFIG mpfr_precision=" << PRECISION
                  << " layers=" << LAYERS
                  << " bessel_series_updates=" << BESSEL_SERIES_TERMS
                  << " rounding=outward\n";
        constexpr int maximum = 592;
        std::vector<mpz_class> factorials(static_cast<std::size_t>(maximum + 1));
        factorials[0] = 1;
        for (int index = 1; index <= maximum; ++index) {
            factorials[static_cast<std::size_t>(index)] =
                factorials[static_cast<std::size_t>(index - 1)] * index;
        }
        if (argc != 1) {
            if ((argc == 6 || argc == 10 || argc == 12)
                && std::string(argv[1]) == "--probe-square-alpha") {
                const std::string family_text = argv[2];
                if (family_text != "B" && family_text != "C") {
                    die("full-Q3 square-alpha probe family must be B or C");
                }
                const char family = family_text[0];
                const int rank = std::stoi(argv[3]);
                ParameterSet parameters{};
                bool found = false;
                for (const LowBCExactSchedule& schedule : LOW_BC_EXACT_SCHEDULES) {
                    if (schedule.parameters.family == family
                        && schedule.parameters.rank_low == rank) {
                        parameters = schedule.parameters;
                        found = true;
                        break;
                    }
                }
                if (!found && family == 'B' && rank == 8) {
                    for (const LowBCExactSchedule& schedule : LOW_BC_EXACT_SCHEDULES) {
                        if (schedule.parameters.family == 'B'
                            && schedule.parameters.rank_low == 9) {
                            parameters = schedule.parameters;
                            parameters.label = "B8-full-Q3";
                            parameters.rank_low = 8;
                            parameters.rank_high = 8;
                            parameters.r = {550, 100};
                            found = true;
                            break;
                        }
                    }
                }
                if (!found && family == 'B' && (rank == 14 || rank == 15)) {
                    parameters = ParameterSet{
                        "B14-15-full-Q3", 'B', rank, rank,
                        {16, 5}, {529, 100}, {213, 50}, {143, 100},
                        {109, 100}, {13, 10}, {31, 1000}, {1, 1}, {1, 1},
                        {1, 1}, {509, 500}, 5, rank - 6, false
                    };
                    found = true;
                }
                if (!found && family == 'B' && (rank == 16 || rank == 17)) {
                    parameters = PARAMETER_SETS[static_cast<std::size_t>(rank - 16)];
                    found = true;
                }
                if (!found && family == 'C' && rank >= 20 && rank <= 28) {
                    parameters = PARAMETER_SETS[3];
                    parameters.rank_low = rank;
                    parameters.rank_high = rank;
                    found = true;
                }
                if (!found) die("full-Q3 square-alpha probe row has no schedule");
                if (family == 'C' && rank == 4) {
                    parameters.r = {33, 10};
                } else if (family == 'C' && rank == 13) {
                    parameters.r = {495, 100};
                } else if (family == 'C' && rank == 14) {
                    parameters.r = {485, 100};
                }
                parameters.square_alpha = {
                    std::stoul(argv[4]), std::stoul(argv[5])
                };
                if (argc >= 10) {
                    parameters.lower_gap = {
                        std::stoul(argv[6]), std::stoul(argv[7])
                    };
                    parameters.upper_start = {
                        std::stoul(argv[8]), std::stoul(argv[9])
                    };
                }
                if (argc == 12) {
                    parameters.r = {
                        std::stoul(argv[10]), std::stoul(argv[11])
                    };
                }
                if (rank < 2 || parameters.square_alpha.denominator == 0
                    || parameters.lower_gap.denominator == 0
                    || parameters.upper_start.denominator == 0
                    || parameters.r.denominator == 0) {
                    die("invalid full-Q3 square-alpha probe");
                }
                const Real probability = exact_full_q3_event_probability(
                    parameters, rank, factorials
                );
                const FullTailResult tail = full_tail_from_probability(
                    family, rank, parameters.r, parameters.character_cutoff,
                    probability, factorials
                );
                print_tail(family, rank, tail, probability);
                return EXIT_SUCCESS;
            }
            if (argc == 6 && std::string(argv[1]) == "--probe-b-template") {
                const int rank = std::stoi(argv[2]);
                const int template_rank = std::stoi(argv[3]);
                ParameterSet parameters{};
                bool found = false;
                for (const LowBCExactSchedule& schedule : LOW_BC_EXACT_SCHEDULES) {
                    if (schedule.parameters.family == 'B'
                        && schedule.parameters.rank_low == template_rank) {
                        parameters = schedule.parameters;
                        found = true;
                        break;
                    }
                }
                if (!found) die("unknown full-Q3 B template row");
                parameters.r = {std::stoul(argv[4]), std::stoul(argv[5])};
                if (rank < 2 || parameters.r.denominator == 0) {
                    die("invalid full-Q3 B template probe");
                }
                const Real probability = exact_full_q3_event_probability(
                    parameters, rank, factorials
                );
                const FullTailResult tail = full_tail_from_probability(
                    'B', rank, parameters.r, parameters.character_cutoff,
                    probability, factorials
                );
                print_tail('B', rank, tail, probability);
                return EXIT_SUCCESS;
            }
            if (argc == 17 && std::string(argv[1]) == "--probe-d-layer") {
                const int rank = std::stoi(argv[2]);
                ParameterSet parameters = PARAMETER_SETS[4];
                parameters.r = {std::stoul(argv[3]), std::stoul(argv[4])};
                parameters.square_alpha = {
                    std::stoul(argv[5]), std::stoul(argv[6])
                };
                parameters.character_cutoff = {
                    std::stoul(argv[7]), std::stoul(argv[8])
                };
                parameters.lower_gap = {
                    std::stoul(argv[9]), std::stoul(argv[10])
                };
                parameters.upper_start = {
                    std::stoul(argv[11]), std::stoul(argv[12])
                };
                parameters.upper_step = {
                    std::stoul(argv[13]), std::stoul(argv[14])
                };
                parameters.square_lambda_fraction = {
                    std::stoul(argv[15]), std::stoul(argv[16])
                };
                if (rank < 4 || parameters.r.denominator == 0
                    || parameters.square_alpha.denominator == 0
                    || parameters.character_cutoff.denominator == 0
                    || parameters.lower_gap.denominator == 0
                    || parameters.upper_start.denominator == 0
                    || parameters.upper_step.denominator == 0
                    || parameters.square_lambda_fraction.denominator == 0) {
                    die("invalid full-Q3 D layered probe");
                }
                const Real probability = exact_full_q3_event_probability(
                    parameters, rank, factorials
                );
                const FullTailResult tail = full_tail_from_probability(
                    'D', rank, parameters.r, parameters.character_cutoff,
                    probability, factorials
                );
                print_tail('D', rank, tail, probability);
                return EXIT_SUCCESS;
            }
            if (argc == 5 && std::string(argv[1]) == "--probe-d-r") {
                const int rank = std::stoi(argv[2]);
                ParameterSet parameters = PARAMETER_SETS[4];
                parameters.r = {
                    std::stoul(argv[3]), std::stoul(argv[4])
                };
                if (rank < 4 || parameters.r.denominator == 0) {
                    die("invalid full-Q3 D r probe");
                }
                const Real probability = exact_full_q3_event_probability(
                    parameters, rank, factorials
                );
                const FullTailResult tail = full_tail_from_probability(
                    'D', rank, parameters.r, parameters.character_cutoff,
                    probability, factorials
                );
                print_tail('D', rank, tail, probability);
                return EXIT_SUCCESS;
            }
            if (argc == 4 && std::string(argv[1]) == "--probe-d-range") {
                const int rank_low = std::stoi(argv[2]);
                const int rank_high = std::stoi(argv[3]);
                if (rank_low < 4 || rank_high < rank_low) {
                    die("invalid full-Q3 D probe range");
                }
                const ParameterSet& parameters = PARAMETER_SETS[4];
                for (int rank = rank_low; rank <= rank_high; ++rank) {
                    const Real probability = exact_full_q3_event_probability(
                        parameters, rank, factorials
                    );
                    const FullTailResult tail = full_tail_from_probability(
                        'D',
                        rank,
                        parameters.r,
                        parameters.character_cutoff,
                        probability,
                        factorials
                    );
                    print_tail('D', rank, tail, probability);
                }
                return EXIT_SUCCESS;
            }
            if (argc == 12 && std::string(argv[1]) == "--probe-layer") {
                const std::string family_text = argv[2];
                if (family_text != "B" && family_text != "C") {
                    die("full-Q3 layer probe family must be B or C");
                }
                const char family = family_text[0];
                const int rank = std::stoi(argv[3]);
                bool found = false;
                ParameterSet parameters{};
                for (const LowBCExactSchedule& schedule : LOW_BC_EXACT_SCHEDULES) {
                    if (schedule.parameters.family == family
                        && schedule.parameters.rank_low == rank) {
                        parameters = schedule.parameters;
                        found = true;
                        break;
                    }
                }
                if (!found) die("full-Q3 layer probe row has no exact schedule");
                parameters.r = {std::stoul(argv[4]), std::stoul(argv[5])};
                parameters.lower_lambda_fraction = {
                    std::stoul(argv[6]), std::stoul(argv[7])
                };
                parameters.upper_start = {
                    std::stoul(argv[8]), std::stoul(argv[9])
                };
                parameters.upper_lambda_fraction = {
                    std::stoul(argv[10]), std::stoul(argv[11])
                };
                if (parameters.r.denominator == 0
                    || parameters.lower_lambda_fraction.denominator == 0
                    || parameters.upper_start.denominator == 0
                    || parameters.upper_lambda_fraction.denominator == 0) {
                    die("zero full-Q3 layer probe denominator");
                }
                const Real probability = exact_full_q3_event_probability(
                    parameters, rank, factorials
                );
                const FullTailResult tail = full_tail_from_probability(
                    family,
                    rank,
                    parameters.r,
                    parameters.character_cutoff,
                    probability,
                    factorials
                );
                print_tail(family, rank, tail, probability);
                return EXIT_SUCCESS;
            }
            if (argc == 12 && std::string(argv[1]) == "--probe-one-step") {
                const std::string family_text = argv[2];
                if (family_text != "B" && family_text != "C"
                    && family_text != "D") {
                    die("full-Q3 one-step probe family must be B, C, or D");
                }
                const char family = family_text[0];
                const int rank = std::stoi(argv[3]);
                const RationalParameter r{
                    std::stoul(argv[4]), std::stoul(argv[5])
                };
                const RationalParameter q{
                    std::stoul(argv[6]), std::stoul(argv[7])
                };
                const RationalParameter central{
                    std::stoul(argv[8]), std::stoul(argv[9])
                };
                const RationalParameter h{
                    std::stoul(argv[10]), std::stoul(argv[11])
                };
                if (r.denominator == 0 || q.denominator == 0
                    || central.denominator == 0 || h.denominator == 0) {
                    die("zero full-Q3 one-step probe denominator");
                }
                const Real probability =
                    exact_one_step_full_q3_event_probability(
                        family, rank, r, q, central, h, factorials
                    );
                const FullTailResult tail = full_tail_from_probability(
                    family, rank, r, central, probability, factorials
                );
                print_tail(family, rank, tail, probability);
                return EXIT_SUCCESS;
            }
            if (argc != 6 || std::string(argv[1]) != "--probe-r") {
                std::cerr << "usage: " << argv[0]
                          << " [--probe-square-alpha B|C RANK NUM DEN"
                          << " [LOW_GAP_NUM LOW_GAP_DEN"
                          << " UPPER_START_NUM UPPER_START_DEN"
                          << " [R_NUM R_DEN]]]"
                          << " [--probe-b-template RANK TEMPLATE_RANK RNUM RDEN]"
                          << " [--probe-d-layer RANK RNUM RDEN QNUM QDEN"
                          << " DNUM DDEN LGNUM LGDEN USNUM USDEN"
                          << " STEPNUM STEPDEN SFNUM SFDEN]"
                          << " [--probe-d-range RANK_LOW RANK_HIGH]"
                          << " [--probe-d-r RANK NUMERATOR DENOMINATOR]"
                          << " [--probe-r B|C RANK NUMERATOR DENOMINATOR]"
                          << " [--probe-layer B|C RANK"
                          << " RNUM RDEN LFNUM LFDEN USNUM USDEN UFNUM UFDEN]"
                          << " [--probe-one-step B|C|D RANK"
                          << " RNUM RDEN QNUM QDEN DNUM DDEN HNUM HDEN]\n";
                return EXIT_FAILURE;
            }
            const std::string family_text = argv[2];
            if (family_text != "B" && family_text != "C") {
                die("full-Q3 r-probe family must be B or C");
            }
            const char family = family_text[0];
            const int rank = std::stoi(argv[3]);
            const unsigned long numerator = std::stoul(argv[4]);
            const unsigned long denominator = std::stoul(argv[5]);
            if (denominator == 0) die("zero full-Q3 r-probe denominator");
            bool found = false;
            ParameterSet parameters{};
            for (const LowBCExactSchedule& schedule : LOW_BC_EXACT_SCHEDULES) {
                if (schedule.parameters.family == family
                    && schedule.parameters.rank_low == rank) {
                    parameters = schedule.parameters;
                    found = true;
                    break;
                }
            }
            if (!found) die("full-Q3 r-probe row has no exact schedule");
            parameters.r = {numerator, denominator};
            const Real probability = exact_full_q3_event_probability(
                parameters, rank, factorials
            );
            const FullTailResult tail = full_tail_from_probability(
                family,
                rank,
                parameters.r,
                parameters.character_cutoff,
                probability,
                factorials
            );
            print_tail(family, rank, tail, probability);
            return EXIT_SUCCESS;
        }
        std::set<std::pair<char, int>> checked_rows;
        for (const LowBCExactSchedule& schedule : LOW_BC_EXACT_SCHEDULES) {
            ParameterSet parameters = schedule.parameters;
            const int rank = parameters.rank_low;
            if (parameters.rank_high != rank) {
                die("low-B/C exact schedule is not a singleton");
            }
            if (parameters.family == 'B' && rank == 7) {
                parameters.square_alpha = {6, 1};
            } else if (parameters.family == 'B' && rank == 9) {
                parameters.square_alpha = {15, 2};
            } else if (parameters.family == 'B' && rank == 11) {
                parameters.r = {47, 10};
                parameters.square_alpha = {7, 1};
                parameters.lower_gap = {3, 2};
                parameters.upper_start = {7, 2};
            } else if (parameters.family == 'B' && rank == 12) {
                parameters.r = {27, 5};
                parameters.square_alpha = {7, 1};
                parameters.lower_gap = {5, 4};
                parameters.upper_start = {7, 2};
            } else if (parameters.family == 'B' && rank == 13) {
                parameters.r = {5, 1};
                parameters.square_alpha = {7, 1};
                parameters.lower_gap = {1, 1};
                parameters.upper_start = {7, 2};
            } else if (parameters.family == 'C' && rank == 4) {
                parameters.r = {33, 10};
            } else if (parameters.family == 'C' && rank == 5) {
                parameters.square_alpha = {9, 2};
            } else if (parameters.family == 'C' && rank == 6) {
                parameters.square_alpha = {5, 1};
            } else if (parameters.family == 'C' && rank == 7) {
                parameters.square_alpha = {11, 2};
            } else if (parameters.family == 'C' && rank == 8) {
                parameters.square_alpha = {6, 1};
            } else if (parameters.family == 'C' && rank == 9) {
                parameters.square_alpha = {15, 2};
            } else if (parameters.family == 'C' && rank == 13) {
                parameters.r = {495, 100};
            } else if (parameters.family == 'C' && rank == 14) {
                parameters.r = {485, 100};
            } else if (parameters.family == 'C' && rank == 17) {
                parameters.square_alpha = {13, 2};
            }
            const Real probability = exact_full_q3_event_probability(
                parameters, rank, factorials
            );
            const FullTailResult tail = full_tail_from_probability(
                parameters.family,
                rank,
                parameters.r,
                parameters.character_cutoff,
                probability,
                factorials
            );
            record_directed_tail(
                parameters.family,
                rank,
                parameters,
                tail,
                probability,
                checked_rows
            );
        }

        ParameterSet b8_parameters{};
        bool found_b8_template = false;
        for (const LowBCExactSchedule& schedule : LOW_BC_EXACT_SCHEDULES) {
            if (schedule.parameters.family == 'B'
                && schedule.parameters.rank_low == 9) {
                b8_parameters = schedule.parameters;
                found_b8_template = true;
                break;
            }
        }
        if (!found_b8_template) die("missing full-Q3 B8 tail template");
        b8_parameters.label = "B8-full-Q3";
        b8_parameters.rank_low = 8;
        b8_parameters.rank_high = 8;
        b8_parameters.r = {550, 100};
        b8_parameters.square_alpha = {17, 2};
        const Real b8_probability = exact_full_q3_event_probability(
            b8_parameters, 8, factorials
        );
        const FullTailResult b8_tail = full_tail_from_probability(
            'B',
            8,
            b8_parameters.r,
            b8_parameters.character_cutoff,
            b8_probability,
            factorials
        );
        record_directed_tail(
            'B', 8, b8_parameters, b8_tail, b8_probability, checked_rows
        );

        for (int rank : {14, 15}) {
            const ParameterSet parameters{
                "B14-15-full-Q3", 'B', rank, rank,
                {16, 5}, {529, 100}, {213, 50}, {143, 100},
                {109, 100}, {13, 10}, {31, 1000}, {1, 1}, {1, 1},
                {1, 1}, {509, 500}, 5, rank - 6, false
            };
            const Real probability = exact_full_q3_event_probability(
                parameters, rank, factorials
            );
            const FullTailResult tail = full_tail_from_probability(
                'B',
                rank,
                RationalParameter{16, 5},
                RationalParameter{143, 100},
                probability,
                factorials
            );
            record_directed_tail(
                'B', rank, parameters, tail, probability, checked_rows
            );
        }

        for (std::size_t parameter_index : {std::size_t(0), std::size_t(1)}) {
            const ParameterSet& parameters = PARAMETER_SETS[parameter_index];
            const int rank = parameters.rank_low;
            const Real probability = exact_full_q3_event_probability(
                parameters, rank, factorials
            );
            const FullTailResult tail = full_tail_from_probability(
                'B',
                rank,
                parameters.r,
                parameters.character_cutoff,
                probability,
                factorials
            );
            record_directed_tail(
                'B', rank, parameters, tail, probability, checked_rows
            );
        }

        const ParameterSet& c_parameters = PARAMETER_SETS[3];
        for (int rank = 20; rank <= 28; ++rank) {
            const Real probability = exact_full_q3_event_probability(
                c_parameters, rank, factorials
            );
            const FullTailResult tail = full_tail_from_probability(
                'C',
                rank,
                c_parameters.r,
                c_parameters.character_cutoff,
                probability,
                factorials
            );
            record_directed_tail(
                'C', rank, c_parameters, tail, probability, checked_rows
            );
        }

        for (int rank = 4; rank <= 30; ++rank) {
            const full_q3_bcd_remaining::RowCutoff* cutoff =
                full_q3_bcd_remaining::find_row('D', rank);
            if (cutoff == nullptr) die("missing full-Q3 D tail row");
            if (cutoff->tail_method
                != full_q3_bcd_remaining::TailMethod::directed_interval) {
                continue;
            }
            ParameterSet parameters = PARAMETER_SETS[4];
            parameters.label = "D8-30-full-Q3";
            parameters.rank_low = rank;
            parameters.rank_high = rank;
            const unsigned long rank_ui = static_cast<unsigned long>(rank);
            if (rank % 2 == 0) {
                parameters.r = {2000001, 1000000};
            } else {
                parameters.r = {
                    2000001UL * static_cast<unsigned long>(rank - 2),
                    1000000UL * rank_ui
                };
            }
            unsigned long square_alpha_tenths = 0;
            unsigned long layer_numerator = 12;
            if (rank == 8) {
                square_alpha_tenths = 50;
            } else if (rank == 10 || rank == 11) {
                square_alpha_tenths = 50;
                layer_numerator = 14;
            } else if (rank == 12) {
                square_alpha_tenths = 50;
            } else if (rank == 14) {
                square_alpha_tenths = 50;
            } else if (rank <= 22) {
                square_alpha_tenths = 45;
            } else if (rank == 23) {
                square_alpha_tenths = 40;
            } else if (rank == 24) {
                square_alpha_tenths = 45;
            } else if (rank == 25) {
                square_alpha_tenths = 40;
            } else if (rank == 26) {
                square_alpha_tenths = 45;
            } else {
                square_alpha_tenths = 40;
            }
            parameters.square_alpha = {
                square_alpha_tenths, 10
            };
            parameters.character_cutoff = {3, 2};
            parameters.lower_gap = {layer_numerator, 10};
            parameters.upper_start = {layer_numerator, 10};
            parameters.upper_step = {3, 100};
            parameters.square_lambda_fraction = {1, 1};
            const Real probability = exact_full_q3_event_probability(
                parameters, rank, factorials
            );
            const FullTailResult tail = full_tail_from_probability(
                'D',
                rank,
                parameters.r,
                parameters.character_cutoff,
                probability,
                factorials
            );
            record_directed_tail(
                'D', rank, parameters, tail, probability, checked_rows
            );
        }

        for (const full_q3_bcd_remaining::RowCutoff& cutoff :
             full_q3_bcd_remaining::row_cutoffs) {
            const bool should_be_checked = cutoff.tail_method
                == full_q3_bcd_remaining::TailMethod::directed_interval;
            const bool was_checked = checked_rows.count(
                std::make_pair(cutoff.family, cutoff.rank)
            ) == 1U;
            if (should_be_checked != was_checked) {
                die("full-Q3 directed-tail coverage mismatch");
            }
        }
        if (checked_rows.size()
            != full_q3_bcd_remaining::directed_interval_rows) {
            die("full-Q3 low-tail row ledger mismatch");
        }
        std::cout << "FULL_Q3_LOW_TAIL rows_checked="
                  << checked_rows.size() << '\n';
        std::cout << "FULL_Q3_BCD_LOW_TAIL VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FULL_Q3_BCD_LOW_TAIL VERIFICATION: FAIL: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
