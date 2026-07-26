#include <array>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <z3++.h>

namespace {

z3::expr zabs(const z3::expr& value) {
    return z3::ite(value >= 0, value, -value);
}

z3::expr zmin(const z3::expr& lhs, const z3::expr& rhs) {
    return z3::ite(lhs <= rhs, lhs, rhs);
}

z3::expr zmax(const z3::expr& lhs, const z3::expr& rhs) {
    return z3::ite(lhs >= rhs, lhs, rhs);
}

z3::expr fusion(
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& target
) {
    return zabs(a - b) <= target
        && target <= a + b
        && target <= 2 * k - a - b
        && (a + b + target) % 2 == 0;
}

z3::expr fourfold_rank(
    z3::context& ctx,
    const z3::expr& k,
    const z3::expr& a,
    const z3::expr& b,
    const z3::expr& c,
    const z3::expr& d
) {
    const z3::expr lo1 = zabs(a - b);
    const z3::expr hi1 = zmin(a + b, 2 * k - a - b);
    const z3::expr lo2 = zabs(c - d);
    const z3::expr hi2 = zmin(c + d, 2 * k - c - d);
    const z3::expr lo = zmax(lo1, lo2);
    const z3::expr hi = zmin(hi1, hi2);
    const z3::expr compatible = lo1 % 2 == lo2 % 2;
    return z3::ite(compatible && lo <= hi, (hi - lo) / 2 + 1,
                   ctx.int_val(0));
}

z3::expr labels_equal(
    z3::context& ctx,
    const std::array<z3::expr, 5>& plus,
    const std::array<z3::expr, 5>& values
) {
    z3::expr result = ctx.bool_val(true);
    for (std::size_t i = 0; i < plus.size(); ++i) {
        result = result && plus[i] == values[i];
    }
    return result;
}

}  // namespace

int main() {
    z3::context ctx;
    z3::solver solver(ctx);

    const z3::expr k = ctx.int_const("k");
    const z3::expr r = ctx.int_const("r");
    const z3::expr s = ctx.int_const("s");
    const std::array<z3::expr, 5> plus{
        ctx.int_const("p0"),
        ctx.int_const("p1"),
        ctx.int_const("p2"),
        ctx.int_const("p3"),
        ctx.int_const("p4")
    };

    solver.add(k >= 6);
    solver.add(r >= 2 && r <= k - 2);
    solver.add(s >= 2 && s <= k - 2);
    for (std::size_t i = 0; i < plus.size(); ++i) {
        solver.add(plus[i] >= 1 && plus[i] <= k - 1);
        solver.add(plus[i] != r && plus[i] != s);
        if (i + 1U < plus.size()) {
            solver.add(plus[i] <= plus[i + 1U]);
        }
    }

    std::vector<z3::expr> active_r;
    std::vector<z3::expr> active_s;
    std::vector<z3::expr> rank_r;
    std::vector<z3::expr> rank_s;
    std::vector<std::array<int, 2>> edges;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            std::array<int, 3> rest{};
            int next = 0;
            for (int h = 0; h < 5; ++h) {
                if (h != i && h != j) {
                    rest[static_cast<std::size_t>(next)] = h;
                    ++next;
                }
            }
            const z3::expr rs = fourfold_rank(
                ctx, k, s, plus[static_cast<std::size_t>(rest[0])],
                plus[static_cast<std::size_t>(rest[1])],
                plus[static_cast<std::size_t>(rest[2])]
            );
            const z3::expr rr = fourfold_rank(
                ctx, k, r, plus[static_cast<std::size_t>(rest[0])],
                plus[static_cast<std::size_t>(rest[1])],
                plus[static_cast<std::size_t>(rest[2])]
            );
            const z3::expr ar = fusion(
                k, plus[static_cast<std::size_t>(i)],
                plus[static_cast<std::size_t>(j)], r
            ) && rs > 0;
            const z3::expr as = fusion(
                k, plus[static_cast<std::size_t>(i)],
                plus[static_cast<std::size_t>(j)], s
            ) && rr > 0;
            solver.add(z3::implies(ar, rs <= 2));
            solver.add(z3::implies(as, rr <= 2));
            active_r.push_back(ar);
            active_s.push_back(as);
            rank_r.push_back(rs);
            rank_s.push_back(rr);
            edges.push_back({i, j});
        }
    }

    z3::expr c = ctx.int_val(0);
    z3::expr h = ctx.int_val(0);
    z3::expr other_c = ctx.int_val(0);
    const z3::expr zero = ctx.int_val(0);
    const z3::expr one_count = ctx.int_val(1);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        c = c + z3::ite(active_r[e], one_count, zero);
        h = h + z3::ite(
            active_r[e] && rank_r[e] == 2, one_count, zero
        );
        other_c = other_c + z3::ite(active_s[e], one_count, zero);
    }
    solver.add(c > 0);
    solver.add(other_c == 0);
    solver.add((h == 0 && c >= 8) || (h > 0 && c + h >= 12));

    std::array<z3::expr, 5> low{
        plus[0] == 1 || plus[0] == k - 1,
        plus[1] == 1 || plus[1] == k - 1,
        plus[2] == 1 || plus[2] == k - 1,
        plus[3] == 1 || plus[3] == k - 1,
        plus[4] == 1 || plus[4] == k - 1
    };
    z3::expr low_count = ctx.int_val(0);
    for (const z3::expr& is_low : low) {
        low_count = low_count + z3::ite(is_low, one_count, zero);
    }
    solver.add(low_count >= 2);

    for (std::size_t e = 0; e < edges.size(); ++e) {
        const int i = edges[e][0];
        const int j = edges[e][1];
        const z3::expr maximal = z3::ite(
            h > 0, active_r[e] && rank_r[e] == 2, active_r[e]
        );
        solver.add(z3::implies(
            low_count == 2 && maximal,
            !low[static_cast<std::size_t>(i)]
                && !low[static_cast<std::size_t>(j)]
        ));
    }

    const z3::expr one = ctx.int_val(1);
    const z3::expr three = ctx.int_val(3);
    const z3::expr jone = k - 1;
    const z3::expr jthree = k - 3;
    const std::array<z3::expr, 5> a0{one, one, one, one, one};
    const std::array<z3::expr, 5> a1{jone, jone, jone, jone, jone};
    const std::array<z3::expr, 5> b0{one, one, one, three, three};
    const std::array<z3::expr, 5> b1{
        jthree, jthree, jone, jone, jone
    };
    const std::array<z3::expr, 5> c0{one, one, jthree, jthree, jone};
    const std::array<z3::expr, 5> c1{one, three, three, jone, jone};

    const std::vector<z3::expr> families{
        r == 2 && s == 3 && labels_equal(ctx, plus, a0),
        r == 2 && s == k - 3 && labels_equal(ctx, plus, a1),
        r == 2 && s == 5 && labels_equal(ctx, plus, b0),
        r == 2 && s == k - 5 && labels_equal(ctx, plus, b1),
        r == k - 2 && s == 3 && labels_equal(ctx, plus, c0),
        r == k - 2 && s == k - 3 && labels_equal(ctx, plus, c1)
    };
    z3::expr classified = ctx.bool_val(false);
    for (const z3::expr& family : families) {
        classified = classified || family;
    }

    if (solver.check() != z3::sat) {
        std::cerr << "base constraints are unexpectedly infeasible\n";
        return EXIT_FAILURE;
    }
    for (std::size_t i = 0; i < families.size(); ++i) {
        solver.push();
        solver.add(families[i]);
        if (solver.check() != z3::sat) {
            std::cerr << "missing family witness index=" << i << '\n';
            return EXIT_FAILURE;
        }
        solver.pop();
    }

    solver.add(!classified);
    const z3::check_result result = solver.check();
    if (result != z3::unsat) {
        std::cerr << "classification failure result=" << result << '\n';
        if (result == z3::sat) {
            std::cerr << solver.get_model() << '\n';
        }
        return EXIT_FAILURE;
    }

    std::cout
        << "SU2_D12_CAP2_Z3 k>=6 base=SAT family_witnesses=6"
        << " counterexample=UNSAT result=PASS\n";
    return EXIT_SUCCESS;
}
