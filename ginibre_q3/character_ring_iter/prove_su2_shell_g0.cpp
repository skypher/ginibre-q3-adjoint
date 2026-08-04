// Exact finite certificate for the diagonal G_0 and G_1 shell families.
//
// In the free band d=K-1-2Q>=11, write d=H+11 and K=2Q+H+12.
// The default G_0 margin on the only nonzero crossing domain V=x is
//
//     (2Q+1)d_5(x)-f_5 d_4(x),
//
// where d_a is evaluated by the three Kac--Walton tail intervals.  This
// Passing --g1 instead certifies f_4 d_4(V)-f_5 d_3(V).  Passing --g2
// certifies f_4 d_3(V)-f_5 d_2(V) on the wide pointwise reserve cone
// 3Q>2d.  The source
// enumerates the exact Presburger endpoint/activation fan and certifies the
// selected rational polynomial in each chamber.  It embeds the general
// three-variable Newton/cone checker so all arithmetic is exact.

#define SU2_T4_GROUPS_EMBEDDED
#include "prove_su2_t4_groups.cpp"
#undef SU2_T4_GROUPS_EMBEDDED

namespace {

using G0Mask = boost::multiprecision::uint128_t;

struct BranchSelection {
    bool nonempty;
    bool lower_zero;
    bool lower_raw;
    bool upper_product;
    bool upper_raw;
    bool lower_minus_one_nonnegative;
    std::vector<bool> upper_images;
    std::vector<bool> lower_images;
};

struct G0Selection {
    bool return_active;
    std::array<std::vector<BranchSelection>, 2> powers;
};

struct G0Record {
    G0Mask mask;
    G0Selection selection;
};

Polynomial choose_binomial(
    const Polynomial& endpoint,
    const Polynomial& label,
    int power,
    int image
) {
    return binomial(
        endpoint
            - scale(label, 2 * image)
            - constant(image)
            + constant(power - 1),
        power - 1
    );
}

void append_constraint(
    std::vector<Polynomial>& constraints,
    const Polynomial& slack,
    bool value
) {
    constraints.push_back(
        value ? slack : constant(-1) - slack
    );
}

Polynomial selected_crossing_power(
    int power,
    int first_power,
    const G0Selection& selection,
    std::vector<Polynomial>& constraints
) {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial x = Polynomial::variable(2);
    // K=2Q+H+12, so H=d-11.
    const Polynomial k = scale(q, 2) + h + constant(12);
    Polynomial result;
    const std::vector<BranchSelection>& branches =
        selection.powers[static_cast<std::size_t>(power - first_power)];
    if (branches.size() != 3U) {
        throw std::runtime_error("unexpected G0 branch count");
    }
    for (int branch = 0; branch < 3; ++branch) {
        const BranchSelection& choice =
            branches[static_cast<std::size_t>(branch)];
        Polynomial lower_label;
        Polynomial upper_label;
        if (branch == 0) {
            lower_label = k - q - x;
            upper_label = q + x;
        } else if (branch == 1) {
            lower_label = scale(k, 2) - q - x + constant(1);
            upper_label = k + q + x + constant(1);
        } else {
            lower_label = scale(k, 3) - q - x + constant(2);
            upper_label = scale(k, 2) + q + x + constant(2);
        }
        const Polynomial product = scale(q, power);
        const Polynomial raw_lower = product - upper_label;
        const Polynomial raw_upper = product - lower_label;
        const Polynomial lower = choice.lower_zero
            ? constant(0) : raw_lower;
        const Polynomial upper = choice.upper_product
            ? product : raw_upper;
        append_constraint(
            constraints, upper - lower, choice.nonempty
        );
        append_constraint(
            constraints, constant(0) - raw_lower, choice.lower_zero
        );
        append_constraint(
            constraints, raw_lower, choice.lower_raw
        );
        append_constraint(
            constraints, raw_upper - product, choice.upper_product
        );
        append_constraint(
            constraints, product - raw_upper, choice.upper_raw
        );
        append_constraint(
            constraints,
            lower - constant(1),
            choice.lower_minus_one_nonnegative
        );
        if (
            choice.upper_images.size()
                != static_cast<std::size_t>(power)
            || choice.lower_images.size()
                != static_cast<std::size_t>(power)
        ) {
            throw std::runtime_error("unexpected G0 activation count");
        }
        for (int image = 1; image <= power; ++image) {
            append_constraint(
                constraints,
                upper - scale(q, 2 * image) - constant(image),
                choice.upper_images[static_cast<std::size_t>(image - 1)]
            );
            append_constraint(
                constraints,
                lower - constant(1) - scale(q, 2 * image)
                    - constant(image),
                choice.lower_images[static_cast<std::size_t>(image - 1)]
            );
        }
        if (!choice.nonempty) {
            continue;
        }
        Polynomial contribution = choose_binomial(upper, q, power, 0);
        if (choice.lower_minus_one_nonnegative) {
            contribution -= choose_binomial(
                lower - constant(1), q, power, 0
            );
        }
        for (int image = 1; image <= power; ++image) {
            const long coefficient = binomial_long(power, image);
            if (choice.upper_images[static_cast<std::size_t>(image - 1)]) {
                contribution += scale(
                    choose_binomial(upper, q, power, image),
                    image % 2 == 0 ? coefficient : -coefficient
                );
            }
            if (choice.lower_images[static_cast<std::size_t>(image - 1)]) {
                contribution -= scale(
                    choose_binomial(
                        lower - constant(1), q, power, image
                    ),
                    image % 2 == 0 ? coefficient : -coefficient
                );
            }
        }
        result += scale(contribution, branch % 2 == 0 ? 2 : -2);
    }
    return result;
}

Chamber make_g0_chamber(
    const G0Record& record,
    int first_power,
    std::uint64_t identifier
) {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial x = Polynomial::variable(2);
    // The active crossing domain.  Off it d_4=d_5=0, so the G_0 margin
    // vanishes and needs no certificate record.
    std::vector<Polynomial> constraints{
        q - constant(7),
        h,
        x,
        scale(q, 2) + h + constant(11) - scale(x, 2),
        scale(x, 2) - h - constant(12)
    };
    if (first_power == 2) {
        // d=H+11 and 3Q>2d are integral, hence 3Q-2H-23>=0.
        constraints.push_back(scale(q, 3) - scale(h, 2) - constant(23));
    }
    // f_5 has its unique affine correction precisely when Q-2d-3>=0.
    const Polynomial return_slack = q - scale(h, 2) - constant(25);
    append_constraint(
        constraints, return_slack, record.selection.return_active
    );
    const Polynomial lower = selected_crossing_power(
        first_power, first_power, record.selection, constraints
    );
    const Polynomial upper = selected_crossing_power(
        first_power + 1, first_power, record.selection, constraints
    );
    const Polynomial f4 = scale(q, 2) + constant(1);
    Polynomial f5 = (
        scale(q * q + q, 5) + constant(2)
    ) * Polynomial(Rational(1, 2));
    if (record.selection.return_active) {
        f5 -= binomial(q - scale(h, 2) - constant(23), 2);
    }
    return Chamber{
        identifier,
        std::move(constraints),
        f4 * upper - f5 * lower
    };
}

std::vector<G0Record> enumerate_g0_records(int first_power) {
    z3::context context;
    z3::solver solver(context);
    const z3::expr q = context.int_const("Q");
    const z3::expr h = context.int_const("H");
    const z3::expr x = context.int_const("X");
    const z3::expr k = 2 * q + h + 12;
    solver.add(q >= 7);
    solver.add(h >= 0);
    solver.add(x >= 0);
    solver.add(2 * x < k);
    solver.add(k - q - x <= q + x);
    if (first_power == 2) {
        solver.add(3 * q > 2 * (h + 11));
    }

    std::vector<z3::expr> hinges;
    const z3::expr return_active = q - 2 * h - 25 >= 0;
    hinges.push_back(return_active);
    struct SymbolicBranch {
        z3::expr nonempty;
        z3::expr lower_zero;
        z3::expr lower_raw;
        z3::expr upper_product;
        z3::expr upper_raw;
        z3::expr lower_minus_one_nonnegative;
        std::vector<z3::expr> upper_images;
        std::vector<z3::expr> lower_images;
    };
    std::array<std::vector<SymbolicBranch>, 2> symbolic{};
    for (int power = first_power; power <= first_power + 1; ++power) {
        std::vector<SymbolicBranch>& power_branches =
            symbolic[static_cast<std::size_t>(power - first_power)];
        for (int branch = 0; branch < 3; ++branch) {
            z3::expr lower_label = context.int_val(0);
            z3::expr upper_label = context.int_val(0);
            if (branch == 0) {
                lower_label = k - q - x;
                upper_label = q + x;
            } else if (branch == 1) {
                lower_label = 2 * k - q - x + 1;
                upper_label = k + q + x + 1;
            } else {
                lower_label = 3 * k - q - x + 2;
                upper_label = 2 * k + q + x + 2;
            }
            const z3::expr product = power * q;
            const z3::expr raw_lower = product - upper_label;
            const z3::expr raw_upper = product - lower_label;
            const z3::expr lower = z3::ite(
                context.int_val(0) >= raw_lower,
                context.int_val(0), raw_lower
            );
            const z3::expr upper = z3::ite(
                product <= raw_upper, product, raw_upper
            );
            SymbolicBranch descriptor{
                lower <= upper,
                lower == 0,
                lower == raw_lower,
                upper == product,
                upper == raw_upper,
                lower - 1 >= 0,
                {}, {}
            };
            hinges.push_back(descriptor.nonempty);
            hinges.push_back(descriptor.lower_zero);
            hinges.push_back(descriptor.lower_raw);
            hinges.push_back(descriptor.upper_product);
            hinges.push_back(descriptor.upper_raw);
            hinges.push_back(descriptor.lower_minus_one_nonnegative);
            for (int image = 1; image <= power; ++image) {
                const z3::expr activation = image * (2 * q + 1);
                descriptor.upper_images.push_back(upper >= activation);
                descriptor.lower_images.push_back(lower - 1 >= activation);
                hinges.push_back(descriptor.upper_images.back());
                hinges.push_back(descriptor.lower_images.back());
            }
            power_branches.push_back(std::move(descriptor));
        }
    }
    const std::size_t expected_hinges = static_cast<std::size_t>(
        1 + 3 * (6 + 2 * first_power)
            + 3 * (6 + 2 * (first_power + 1))
    );
    if (hinges.size() != expected_hinges) {
        throw std::runtime_error("unexpected diagonal hinge count");
    }

    std::vector<G0Record> result;
    while (solver.check() == z3::sat) {
        const z3::model model = solver.get_model();
        G0Mask mask = 0;
        std::size_t bit = 0U;
        const auto read = [&model, &context, &mask, &bit](
            const z3::expr& expression
        ) {
            const bool value = z3::eq(
                model.eval(expression, true), context.bool_val(true)
            );
            if (value) {
                mask |= G0Mask(1) << bit;
            }
            ++bit;
            return value;
        };
        G0Selection selection;
        selection.return_active = read(return_active);
        for (int power = first_power;
             power <= first_power + 1;
             ++power) {
            std::vector<BranchSelection>& selected_branches =
                selection.powers[static_cast<std::size_t>(
                    power - first_power
                )];
            for (const SymbolicBranch& branch
                 : symbolic[static_cast<std::size_t>(
                       power - first_power
                   )]) {
                BranchSelection selected{
                    read(branch.nonempty),
                    read(branch.lower_zero),
                    read(branch.lower_raw),
                    read(branch.upper_product),
                    read(branch.upper_raw),
                    read(branch.lower_minus_one_nonnegative),
                    {}, {}
                };
                for (int image = 0; image < power; ++image) {
                    selected.upper_images.push_back(read(
                        branch.upper_images[static_cast<std::size_t>(image)]
                    ));
                    selected.lower_images.push_back(read(
                        branch.lower_images[static_cast<std::size_t>(image)]
                    ));
                }
                selected_branches.push_back(std::move(selected));
            }
        }
        if (bit != hinges.size()) {
            throw std::runtime_error("diagonal hinge read mismatch");
        }
        z3::expr block = context.bool_val(false);
        for (std::size_t index = 0U; index < hinges.size(); ++index) {
            const bool value = (mask & (G0Mask(1) << index)) != 0;
            block = block || (value ? !hinges[index] : hinges[index]);
        }
        solver.add(block);
        result.push_back(G0Record{mask, std::move(selection)});
    }
    std::sort(
        result.begin(), result.end(),
        [](const G0Record& left, const G0Record& right) {
            return left.mask < right.mask;
        }
    );
    return result;
}

bool certify_g0_chamber(const Chamber& chamber) {
    const std::vector<Polynomial> constraints =
        irredundant_constraints(chamber);
    const std::vector<bool> forced_zero =
        forced_zero_constraints(constraints);
    bool passed = bounded_group_integer_certificate(chamber, constraints);
    if (!passed) {
        passed = bounded_qy_h_ray_newton_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = forced_equality_slice_ray_newton_certificate(
            chamber, constraints
        );
    }
    if (!passed) {
        passed = parity_unit_offset_square_cone_certificate(
            chamber, constraints
        );
    }
    if (!passed) {
        passed = single_recession_ray_newton_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = recession_fan_newton_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = direct_facet_certificate(chamber, constraints, forced_zero);
    }
    if (!passed) {
        passed = sum_cone_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = constant_three_sum_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = equal_sum_square_cone_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = double_sum_cone_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = pair_cut_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = pair_cut_tree_certificate(chamber, constraints);
    }
    return passed;
}

Integer integer_binomial(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

Integer tensor_weight(int power, int label, int index) {
    Integer result = 0;
    for (int image = 0; image <= power; ++image) {
        const Integer term = binomial_long(power, image)
            * integer_binomial(
                index - image * (2 * label + 1) + power - 1,
                power - 1
            );
        result += image % 2 == 0 ? term : -term;
    }
    return result;
}

Integer direct_crossing_weight(int k, int q, int power, int x) {
    const int product = power * q;
    const int lower_depth = k - q - x;
    const int upper_depth = q + x;
    if (lower_depth > upper_depth) {
        return 0;
    }
    Integer result = 0;
    for (int branch = 0; branch < 3; ++branch) {
        int lower_label = 0;
        int upper_label = 0;
        int sign = 1;
        if (branch == 0) {
            lower_label = lower_depth;
            upper_label = upper_depth;
        } else if (branch == 1) {
            lower_label = 2 * k - q - x + 1;
            upper_label = k + q + x + 1;
            sign = -1;
        } else {
            lower_label = 3 * k - q - x + 2;
            upper_label = 2 * k + q + x + 2;
        }
        const int lower = std::max(0, product - upper_label);
        const int upper = std::min(product, product - lower_label);
        if (lower <= upper) {
            result += sign * (
                tensor_weight(power, q, upper)
                - tensor_weight(power, q, lower - 1)
            );
        }
    }
    return 2 * result;
}

void audit_g0_formula(
    const std::vector<G0Record>& records,
    int first_power,
    int maximum_k
) {
    std::vector<Chamber> chambers;
    chambers.reserve(records.size());
    for (std::size_t index = 0U; index < records.size(); ++index) {
        chambers.push_back(make_g0_chamber(
            records[index], first_power, static_cast<std::uint64_t>(index)
        ));
    }
    std::uint64_t parameters = 0U;
    std::uint64_t entries = 0U;
    std::uint64_t inactive_entries = 0U;
    for (int k = 3; k <= maximum_k; ++k) {
        for (int q = 7; 2 * q < k; ++q) {
            const int d = k - 1 - 2 * q;
            if (d < 11) {
                continue;
            }
            if (first_power == 2 && 3 * q <= 2 * d) {
                continue;
            }
            ++parameters;
            const Integer f4 = 2 * q + 1;
            const Integer f5 = (5 * Integer(q) * q + 5 * q + 2) / 2
                - integer_binomial(q - 2 * d - 1, 2);
            for (int x = 0; 2 * x < k; ++x) {
                const Integer direct = f4 * direct_crossing_weight(
                        k, q, first_power + 1, x
                    ) - f5 * direct_crossing_weight(k, q, first_power, x);
                const std::array<int, 3> values{q, d - 11, x};
                std::size_t matches = 0U;
                for (const Chamber& chamber : chambers) {
                    if (!std::all_of(
                            chamber.constraints.begin(),
                            chamber.constraints.end(),
                            [&values](const Polynomial& constraint) {
                                return evaluate(constraint, values) >= 0;
                            }
                        )) {
                        continue;
                    }
                    ++matches;
                    const Rational margin = evaluate(chamber.margin, values);
                    if (margin.denominator() != 1 || margin.numerator() != direct) {
                        throw std::runtime_error(
                            "diagonal polynomial/direct mismatch K="
                            + std::to_string(k)
                            + " Q=" + std::to_string(q)
                            + " x=" + std::to_string(x)
                        );
                    }
                    // Tied endpoint selectors can place a coordinate in
                    // several records.  The first exact record already
                    // validates the formula; the chamber certificate later
                    // proves every tied record independently.
                    break;
                }
                if (matches == 0U) {
                    if (direct != 0) {
                        throw std::runtime_error(
                            "uncovered nonzero diagonal coordinate K="
                            + std::to_string(k)
                            + " Q=" + std::to_string(q)
                            + " x=" + std::to_string(x)
                        );
                    }
                    ++inactive_entries;
                }
                ++entries;
            }
        }
    }
    std::cout
        << "SU2_SHELL_DIAGONAL_AUDIT"
        << " first_power=" << first_power
        << " maximum_K=" << maximum_k
        << " parameters=" << parameters
        << " entries=" << entries
        << " inactive_zero_entries=" << inactive_entries
        << " result=PASS_EXACT_REPLAY\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int first_power = 4;
        int argument = 1;
        if (argc > 1 && std::string(argv[1]) == "--g1") {
            first_power = 3;
            ++argument;
        } else if (argc > 1 && std::string(argv[1]) == "--g2") {
            first_power = 2;
            ++argument;
        }
        const std::string target = first_power == 4
            ? "G0"
            : (first_power == 3 ? "G1" : "H2_EVEN_P2_MARGIN");
        std::optional<std::size_t> selected_position;
        std::optional<std::pair<std::size_t, std::size_t>> selected_range;
        std::optional<int> audit_bound;
        bool masks_only = false;
        const int remaining = argc - argument;
        if (
            remaining == 1
            && std::string(argv[argument]) == "--masks-only"
        ) {
            masks_only = true;
        } else if (
            remaining == 2
            && std::string(argv[argument]) == "--position"
        ) {
            selected_position = static_cast<std::size_t>(
                std::stoull(argv[argument + 1])
            );
        } else if (
            remaining == 3
            && std::string(argv[argument]) == "--range"
        ) {
            const std::size_t begin = static_cast<std::size_t>(
                std::stoull(argv[argument + 1])
            );
            const std::size_t end = static_cast<std::size_t>(
                std::stoull(argv[argument + 2])
            );
            if (begin >= end) {
                throw std::runtime_error(
                    "diagonal range must satisfy BEGIN < END"
                );
            }
            selected_range = std::make_pair(begin, end);
        } else if (
            remaining == 2
            && std::string(argv[argument]) == "--audit"
        ) {
            audit_bound = std::stoi(argv[argument + 1]);
            if (*audit_bound < 3) {
                throw std::runtime_error("audit bound must be at least three");
            }
        } else if (remaining != 0) {
            throw std::runtime_error(
                "usage: [--g1|--g2] [--masks-only|--position INDEX|"
                "--range BEGIN END|--audit MAXIMUM_K]"
            );
        }
        const std::vector<G0Record> records =
            enumerate_g0_records(first_power);
        const std::size_t hinge_count = static_cast<std::size_t>(
            1 + 3 * (6 + 2 * first_power)
                + 3 * (6 + 2 * (first_power + 1))
        );
        std::cout
            << "SU2_SHELL_" << target << "_MASKS"
            << " d_min=11"
            << " hinges=" << hinge_count
            << " masks=" << records.size()
            << " result=PASS_EXACT_PRESBURGER_CENSUS\n";
        if (masks_only) {
            return EXIT_SUCCESS;
        }
        if (audit_bound.has_value()) {
            audit_g0_formula(records, first_power, *audit_bound);
            return EXIT_SUCCESS;
        }
        if (
            selected_position.has_value()
            && *selected_position >= records.size()
        ) {
            throw std::runtime_error("diagonal chamber position is out of range");
        }
        if (
            selected_range.has_value()
            && selected_range->second > records.size()
        ) {
            throw std::runtime_error("diagonal chamber range is out of range");
        }
        std::size_t attempted = 0U;
        std::size_t certified = 0U;
        for (std::size_t index = 0U; index < records.size(); ++index) {
            if (selected_position.has_value() && index != *selected_position) {
                continue;
            }
            if (
                selected_range.has_value()
                && (
                    index < selected_range->first
                    || index >= selected_range->second
                )
            ) {
                continue;
            }
            ++attempted;
            const Chamber chamber = make_g0_chamber(
                records[index], first_power, static_cast<std::uint64_t>(index)
            );
            if (certify_g0_chamber(chamber)) {
                ++certified;
            } else {
                std::cout
                    << "SU2_SHELL_" << target << "_UNRESOLVED"
                    << " position=" << index
                    << " mask=" << records[index].mask
                    << " result=NEEDS_STRONGER_CERTIFICATE\n";
            }
            std::cerr
                << "SU2_SHELL_" << target << "_PROGRESS"
                << " position=" << index + 1U
                << '/' << records.size()
                << " certified=" << certified << '\n';
        }
        const bool complete = certified == attempted;
        std::cout
            << "SU2_SHELL_" << target
            << " d_min=11"
            << " attempted=" << attempted
            << " certified=" << certified
            << " scope="
            << (
                selected_position.has_value()
                    ? "position"
                    : (selected_range.has_value() ? "range" : "all")
            )
            << " result="
            << (complete ? "PASS_EXACT_CERTIFICATE" : "INCOMPLETE")
            << '\n';
        return complete ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SHELL_DIAGONAL FAILURE: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
