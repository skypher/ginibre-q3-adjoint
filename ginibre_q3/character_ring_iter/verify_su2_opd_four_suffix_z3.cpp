#include <array>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>

#include <z3++.h>

namespace {

z3::expr absolute(const z3::expr& value) {
    return z3::ite(value >= 0, value, -value);
}

z3::expr maximum(const z3::expr& first, const z3::expr& second) {
    return z3::ite(first >= second, first, second);
}

z3::expr minimum(const z3::expr& first, const z3::expr& second) {
    return z3::ite(first <= second, first, second);
}

z3::expr maximum3(
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& third
) {
    return maximum(maximum(first, second), third);
}

z3::expr indicator(const z3::expr& proposition) {
    return z3::ite(
        proposition,
        proposition.ctx().int_val(1),
        proposition.ctx().int_val(0)
    );
}

z3::expr fusion_contains(
    const z3::expr& left,
    const z3::expr& right,
    const z3::expr& output
) {
    return output >= absolute(left - right)
        && output <= left + right
        && z3::mod(left + right - output, 2) == 0;
}

z3::expr truncated_interval_cardinality(
    const z3::expr& low,
    const z3::expr& high,
    const z3::expr& parity,
    int cap
) {
    z3::context& context = low.ctx();
    const z3::expr exact = z3::ite(
        parity && low <= high,
        (high - low) / 2 + 1,
        context.int_val(0)
    );
    return z3::ite(exact <= cap, exact, context.int_val(cap));
}

z3::expr fusion_interval_intersection_cardinality(
    const z3::expr& first_left,
    const z3::expr& first_right,
    const z3::expr& second_left,
    const z3::expr& second_right
) {
    const z3::expr low = maximum(
        absolute(first_left - first_right),
        absolute(second_left - second_right)
    );
    const z3::expr high = minimum(
        first_left + first_right,
        second_left + second_right
    );
    const z3::expr parity = z3::mod(
        first_left + first_right - second_left - second_right,
        2
    ) == 0;
    return z3::ite(
        parity && low <= high,
        (high - low) / 2 + 1,
        low.ctx().int_val(0)
    );
}

z3::expr triple_support_low(
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& third
) {
    const z3::expr total = first + second + third;
    const z3::expr defect
        = 2 * maximum3(first, second, third) - total;
    return z3::ite(defect > 0, defect, z3::mod(total, 2));
}

z3::expr truncated_threefold_multiplicity(
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& third,
    const z3::expr& target,
    int cap
) {
    const z3::expr low = maximum(
        absolute(first - second), absolute(third - target)
    );
    const z3::expr high = minimum(first + second, third + target);
    const z3::expr parity
        = z3::mod(first + second - third - target, 2) == 0;
    return z3::ite(
        target >= 0,
        truncated_interval_cardinality(low, high, parity, cap),
        first.ctx().int_val(0)
    );
}

z3::expr main_path_candidates(
    int p_value,
    const z3::expr& q,
    const z3::expr& r,
    const z3::expr& s,
    const z3::expr& t,
    const z3::expr& a,
    int second_cap,
    int third_cap
) {
    z3::context& context = q.ctx();
    const z3::expr p = context.int_val(p_value);
    const z3::expr suffix_total = s + t + a;
    const z3::expr suffix_low = triple_support_low(s, t, a);
    z3::expr result = context.int_val(0);
    for (int first_index = 0; first_index < p_value; ++first_index) {
        const z3::expr first = q - p + 1 + 2 * first_index;
        const z3::expr second_low = maximum(
            absolute(first - r), suffix_low
        );
        const z3::expr second_high = minimum(first + r, suffix_total);
        const z3::expr second_parity
            = z3::mod(first + r - suffix_total, 2) == 0;
        for (int second_index = 0; second_index < second_cap;
             ++second_index) {
            const z3::expr second = second_low + 2 * second_index;
            const z3::expr third_low = maximum(
                absolute(second - s), absolute(a - t)
            );
            const z3::expr third_high = minimum(second + s, a + t);
            const z3::expr third_parity
                = z3::mod(second + s - a - t, 2) == 0;
            result = result + z3::ite(
                second_parity && second <= second_high,
                truncated_interval_cardinality(
                    third_low, third_high, third_parity, third_cap
                ),
                context.int_val(0)
            );
        }
    }
    return result;
}

z3::expr fourfold_invariant_candidates(
    const std::array<z3::expr, 4>& labels,
    int cap
) {
    const z3::expr low = maximum(
        absolute(labels[0] - labels[1]),
        absolute(labels[2] - labels[3])
    );
    const z3::expr high = minimum(
        labels[0] + labels[1], labels[2] + labels[3]
    );
    const z3::expr parity = z3::mod(
        labels[0] + labels[1] - labels[2] - labels[3], 2
    ) == 0;
    return truncated_interval_cardinality(low, high, parity, cap);
}

z3::expr fourfold_target_candidates(
    const z3::expr& q,
    const z3::expr& r,
    const z3::expr& s,
    const z3::expr& t,
    int target,
    int cap
) {
    z3::context& context = q.ctx();
    const z3::expr target_expr = context.int_val(target);
    const z3::expr support_total = s + t + target_expr;
    const z3::expr support_low = triple_support_low(s, t, target_expr);
    const z3::expr first_low = maximum(
        absolute(q - r), support_low
    );
    const z3::expr first_high = minimum(q + r, support_total);
    const z3::expr first_parity
        = z3::mod(q + r - support_total, 2) == 0;
    z3::expr result = context.int_val(0);
    for (int first_index = 0; first_index < cap; ++first_index) {
        const z3::expr first = first_low + 2 * first_index;
        const z3::expr second_low = maximum(
            absolute(s - t), absolute(first - target_expr)
        );
        const z3::expr second_high = minimum(s + t, first + target_expr);
        const z3::expr second_parity
            = z3::mod(s + t - first - target_expr, 2) == 0;
        result = result + z3::ite(
            first_parity && first <= first_high,
            truncated_interval_cardinality(
                second_low, second_high, second_parity, cap
            ),
            context.int_val(0)
        );
    }
    return result;
}

z3::expr negative_load(
    const z3::expr& p,
    const std::array<z3::expr, 4>& labels,
    const z3::expr& a
) {
    z3::context& context = a.ctx();
    z3::expr result = context.int_val(0);
    for (int first = 0; first < 4; ++first) {
        for (int second = first + 1; second < 4; ++second) {
            std::array<int, 2> complement{};
            int complement_size = 0;
            for (int index = 0; index < 4; ++index) {
                if (index != first && index != second) {
                    complement[static_cast<std::size_t>(complement_size)]
                        = index;
                    ++complement_size;
                }
            }
            result = result + indicator(
                fusion_contains(labels[static_cast<std::size_t>(first)],
                                labels[static_cast<std::size_t>(second)],
                                p - 1)
                && fusion_contains(
                    labels[static_cast<std::size_t>(complement[0])],
                    labels[static_cast<std::size_t>(complement[1])], a
                )
            );
            const z3::expr pair_target_multiplicity
                = fusion_interval_intersection_cardinality(
                    labels[static_cast<std::size_t>(complement[0])],
                    labels[static_cast<std::size_t>(complement[1])],
                    a,
                    p
                );
            result = result + z3::ite(
                fusion_contains(
                    labels[static_cast<std::size_t>(first)],
                    labels[static_cast<std::size_t>(second)],
                    context.int_val(1)
                ),
                pair_target_multiplicity,
                context.int_val(0)
            );
        }
    }
    for (int omitted = 0; omitted < 4; ++omitted) {
        std::array<int, 3> triple{};
        int triple_size = 0;
        for (int index = 0; index < 4; ++index) {
            if (index != omitted) {
                triple[static_cast<std::size_t>(triple_size)] = index;
                ++triple_size;
            }
        }
        const z3::expr& first
            = labels[static_cast<std::size_t>(triple[0])];
        const z3::expr& second
            = labels[static_cast<std::size_t>(triple[1])];
        const z3::expr& third
            = labels[static_cast<std::size_t>(triple[2])];
        const z3::expr triple_to_one
            = indicator(fusion_contains(first, second, third - 1))
            + indicator(fusion_contains(first, second, third + 1));
        result = result + z3::ite(
            fusion_contains(
                labels[static_cast<std::size_t>(omitted)], p, a
            ),
            triple_to_one,
            context.int_val(0)
        );
    }
    return result;
}

z3::expr boundary_payments(
    int p_value,
    const std::array<z3::expr, 4>& labels,
    const z3::expr& a,
    int cap
) {
    z3::context& context = a.ctx();
    const z3::expr p = context.int_val(p_value);
    z3::expr singleton = context.int_val(0);
    for (int omitted = 0; omitted < 4; ++omitted) {
        std::array<int, 3> complement{};
        int complement_size = 0;
        for (int index = 0; index < 4; ++index) {
            if (index != omitted) {
                complement[static_cast<std::size_t>(complement_size)] = index;
                ++complement_size;
            }
        }
        const z3::expr multiplicity = truncated_threefold_multiplicity(
                labels[static_cast<std::size_t>(complement[0])],
                labels[static_cast<std::size_t>(complement[1])],
                labels[static_cast<std::size_t>(complement[2])],
                a - 1,
                cap
            )
            + truncated_threefold_multiplicity(
                labels[static_cast<std::size_t>(complement[0])],
                labels[static_cast<std::size_t>(complement[1])],
                labels[static_cast<std::size_t>(complement[2])],
                a + 1,
                cap
            );
        singleton = singleton + z3::ite(
            labels[static_cast<std::size_t>(omitted)] == p,
            multiplicity,
            context.int_val(0)
        );
    }
    const z3::expr invariant = fourfold_invariant_candidates(labels, cap);
    const z3::expr target = fourfold_target_candidates(
        labels[0], labels[1], labels[2], labels[3], p_value, cap
    );
    return singleton
        + z3::ite(a == p - 1, invariant, context.int_val(0))
        + z3::ite(a == 1, target, context.int_val(0));
}

void print_model_value(
    const z3::model& model,
    const char* name,
    const z3::expr& value
) {
    std::cout << ' ' << name << '=' << model.eval(value, true);
}

bool verify_negative_bound() {
    z3::context context;
    z3::tactic qflia(context, "qflia");
    z3::solver solver = qflia.mk_solver();
    const z3::expr p = context.int_const("p");
    const z3::expr q = context.int_const("q");
    const z3::expr r = context.int_const("r");
    const z3::expr s = context.int_const("s");
    const z3::expr t = context.int_const("t");
    const z3::expr a = context.int_const("a");
    const std::array<z3::expr, 4> labels{q, r, s, t};
    solver.add(
        p >= 2 && p <= q && q <= r && r <= s && s <= t && a >= 1
    );
    solver.add(z3::mod(p - 1 + q + r + s + t - a, 2) == 0);
    for (const z3::expr& label : labels) {
        solver.add(a != label);
    }
    const z3::expr negative = negative_load(p, labels, a);
    solver.add(negative > 4 * p + 12);
    std::cout << "solving symbolic negative bound N<=4p+12" << std::endl;
    const z3::check_result result = solver.check();
    if (result == z3::sat) {
        const z3::model model = solver.get_model();
        std::cout << "SU2_OPD_FOUR_SUFFIX_NEGATIVE_BOUND_Z3 FAIL";
        print_model_value(model, "p", p);
        print_model_value(model, "q", q);
        print_model_value(model, "r", r);
        print_model_value(model, "s", s);
        print_model_value(model, "t", t);
        print_model_value(model, "a", a);
        print_model_value(model, "negative", negative);
        std::cout << '\n';
        return false;
    }
    if (result == z3::unknown) {
        std::cout << "SU2_OPD_FOUR_SUFFIX_NEGATIVE_BOUND_Z3 UNKNOWN reason="
                  << solver.reason_unknown() << '\n';
        return false;
    }
    std::cout << "SU2_OPD_FOUR_SUFFIX_NEGATIVE_BOUND_Z3 PASS result=unsat\n";
    return true;
}

bool verify_p4_endpoint_exclusion() {
    z3::context context;
    z3::tactic qflia(context, "qflia");
    z3::solver solver = qflia.mk_solver();
    const z3::expr q = context.int_const("q");
    const z3::expr r = context.int_const("r");
    const z3::expr s = context.int_const("s");
    const z3::expr t = context.int_const("t");
    const z3::expr a = context.int_const("a");
    const std::array<z3::expr, 4> labels{q, r, s, t};
    solver.add(4 <= q && q <= r && r <= s && s <= t && a >= 1);
    solver.add(z3::mod(3 + q + r + s + t - a, 2) == 0);
    for (const z3::expr& label : labels) {
        solver.add(a != label);
    }
    const z3::expr negative = negative_load(context.int_val(4), labels, a);
    solver.add(negative > 24 && (a == 1 || a == 3));
    std::cout << "solving p=4 residue endpoint exclusion" << std::endl;
    const z3::check_result result = solver.check();
    if (result == z3::sat) {
        const z3::model model = solver.get_model();
        std::cout << "SU2_OPD_FOUR_SUFFIX_P4_ENDPOINT_Z3 FAIL";
        print_model_value(model, "q", q);
        print_model_value(model, "r", r);
        print_model_value(model, "s", s);
        print_model_value(model, "t", t);
        print_model_value(model, "a", a);
        print_model_value(model, "negative", negative);
        std::cout << '\n';
        return false;
    }
    if (result == z3::unknown) {
        std::cout << "SU2_OPD_FOUR_SUFFIX_P4_ENDPOINT_Z3 UNKNOWN reason="
                  << solver.reason_unknown() << '\n';
        return false;
    }
    std::cout << "SU2_OPD_FOUR_SUFFIX_P4_ENDPOINT_Z3 PASS result=unsat\n";
    return true;
}

bool verify_spin(int p_value, int case_index) {
    z3::context context;
    z3::tactic qflia(context, "qflia");
    z3::solver solver = qflia.mk_solver();
    const z3::expr q = context.int_const("q");
    const z3::expr r = context.int_const("r");
    const z3::expr s = context.int_const("s");
    const z3::expr t = context.int_const("t");
    const z3::expr a = context.int_const("a");
    const std::array<z3::expr, 4> labels{q, r, s, t};
    solver.add(q >= p_value && q <= r && r <= s && s <= t && a >= 1);
    solver.add(
        z3::mod(p_value - 1 + q + r + s + t - a, 2) == 0
    );
    for (const z3::expr& label : labels) {
        solver.add(a != label);
    }
    std::string case_name = "all";
    if (p_value == 2) {
        if (case_index == 0) {
            case_name = "endpoint";
            solver.add(a == 1);
        } else if (case_index == 1) {
            case_name = "strict-interior-low-target";
            solver.add(a > 1 && a < q && q > 2);
        } else if (case_index == 2) {
            case_name = "strict-interior-high-target";
            solver.add(a > q && q > 2);
        } else if (case_index == 3) {
            case_name = "repeated-minimum-interior";
            solver.add(a > 1 && q == 2);
        } else {
            throw std::runtime_error("invalid p=2 case index");
        }
    } else if (case_index == 0 || case_index == 1) {
        const bool endpoint = case_index == 0;
        case_name = endpoint
            ? "repeated-minimum-endpoint"
            : "repeated-minimum-interior";
        solver.add(q == p_value);
        solver.add(endpoint
            ? (a == 1 || a == p_value - 1)
            : (a != 1 && a != p_value - 1));
    } else if (case_index == 2) {
        case_name = "strict-minimum-endpoint";
        solver.add(q > p_value);
        solver.add(a == 1 || a == p_value - 1);
    } else if (case_index == 3) {
        case_name = "strict-minimum-interior-low-target";
        solver.add(q > p_value && a > 1 && a != p_value - 1 && a < q);
    } else if (case_index == 4) {
        case_name = "strict-minimum-interior-high-target";
        solver.add(q > p_value && a > q);
    } else {
        throw std::runtime_error("invalid small-spin case index");
    }
    // This is only the depth of each explicit path family.  The solver checks
    // the candidate inequality directly; no upper bound on the negative load
    // is assumed as an axiom.
    const int negative_cap = 4 * p_value + 12;
    const z3::expr negative = negative_load(
        context.int_val(p_value), labels, a
    );
    const int second_cap
        = (negative_cap + p_value - 1) / p_value;
    const int third_cap = 4;
    const z3::expr candidates = main_path_candidates(
        p_value, q, r, s, t, a, second_cap, third_cap
    );
    const z3::expr payments = p_value == 2 && case_index != 1
        ? boundary_payments(p_value, labels, a, negative_cap)
        : context.int_val(0);
    // Lemma 25 pays every channel-1 case with N<=p(p+2), while an N_A-only
    // load is paid by Proposition 1C.  Only the strict excess remains here.
    solver.add(negative > p_value * (p_value + 2));
    solver.add(candidates + payments < negative);
    std::cout << "solving p=" << p_value << " case=" << case_name
              << " caps=" << second_cap << 'x' << third_cap << std::endl;
    std::mutex progress_mutex;
    std::condition_variable progress_condition;
    bool finished = false;
    int elapsed_minutes = 0;
    std::jthread progress([&]() {
        std::unique_lock<std::mutex> lock(progress_mutex);
        while (!progress_condition.wait_for(
            lock,
            std::chrono::minutes(1),
            [&]() { return finished; }
        )) {
            ++elapsed_minutes;
            std::cout << "progress p=" << p_value
                      << " case=" << case_name
                      << " elapsed_minutes=" << elapsed_minutes
                      << std::endl;
        }
    });
    const z3::check_result result = solver.check();
    {
        std::lock_guard<std::mutex> lock(progress_mutex);
        finished = true;
    }
    progress_condition.notify_one();
    progress.join();
    if (result == z3::sat) {
        const z3::model model = solver.get_model();
        std::cout << "SU2_OPD_FOUR_SUFFIX_Z3 FAIL";
        print_model_value(model, "p", context.int_val(p_value));
        print_model_value(model, "q", q);
        print_model_value(model, "r", r);
        print_model_value(model, "s", s);
        print_model_value(model, "t", t);
        print_model_value(model, "a", a);
        print_model_value(model, "candidates", candidates);
        print_model_value(model, "payments", payments);
        print_model_value(model, "negative", negative);
        std::cout << " case=" << case_name;
        std::cout << '\n';
        return false;
    }
    if (result == z3::unknown) {
        std::cout << "SU2_OPD_FOUR_SUFFIX_Z3 UNKNOWN p=" << p_value
                  << " case=" << case_name
                  << " reason=" << solver.reason_unknown() << '\n';
        return false;
    }
    std::cout << "p=" << p_value << " case=" << case_name
              << " result=unsat caps=" << second_cap << 'x' << third_cap
              << '\n';
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--negative-bound") {
            return verify_negative_bound() ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 2 && std::string(argv[1]) == "--p4-endpoints") {
            return verify_p4_endpoint_exclusion()
                ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        int first_spin = 2;
        int last_spin = 4;
        int requested_case = -1;
        if (argc == 2 || argc == 3) {
            first_spin = std::stoi(argv[1]);
            last_spin = first_spin;
            if (argc == 3) {
                requested_case = std::stoi(argv[2]);
            }
        } else if (argc != 1) {
            throw std::runtime_error(
                "usage: verify_su2_opd_four_suffix_z3 [SPIN [CASE_INDEX]]"
            );
        }
        if (first_spin < 2 || first_spin > 4) {
            throw std::runtime_error("spin must lie in 2,3,4");
        }
        for (int p = first_spin; p <= last_spin; ++p) {
            const int case_count = p == 2 ? 4 : 5;
            const int first_case = requested_case < 0 ? 0 : requested_case;
            const int last_case
                = requested_case < 0 ? case_count - 1 : requested_case;
            if (first_case < 0 || last_case >= case_count) {
                throw std::runtime_error("invalid case index for spin");
            }
            for (int case_index = first_case; case_index <= last_case;
                 ++case_index) {
                if (!verify_spin(p, case_index)) {
                    return EXIT_FAILURE;
                }
            }
        }
        std::cout << "SU2_OPD_FOUR_SUFFIX_Z3 PASS spins=" << first_spin;
        if (last_spin != first_spin) {
            std::cout << ",...," << last_spin;
        }
        std::cout << '\n';
        return EXIT_SUCCESS;
    } catch (const z3::exception& error) {
        std::cerr << error.msg() << '\n';
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
