"""Exact finite audit of the q=2 plus row inequality. No files are written."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from math import comb


def row(a, e):
    n = a + e
    return tuple(sum((-1)**t * comb(e, t) * comb(a, k-t)
                     for t in range(e+1) if 0 <= k-t <= a)
                 for k in range(n+1))


def row_value(c, n, k):
    return c[k] if 0 <= k <= n else 0


def slack_and_W(c, n, j):
    C = lambda k: row_value(c, n, k)
    D = lambda k: C(k)**2 - C(k-1)*C(k+1)
    B = lambda k: C(k-1) + C(k+1)
    W = B(j+3)*C(j) - C(j+3)*B(j)
    return D(j)-D(j+3)+W, W


def slack_from_initial(n, d, j, x, y):
    c = {j: Fraction(x), j+1: Fraction(y)}
    for k in (j+1, j+2, j+3):
        c[k+1] = (d*c[k] - (n-k+1)*c[k-1]) / (k+1)
    c[j-1] = (d*c[j] - (j+1)*c[j+1]) / (n-j+1)
    C = lambda k: c.get(k, Fraction(0))
    D = lambda k: C(k)**2 - C(k-1)*C(k+1)
    B = lambda k: C(k-1) + C(k+1)
    return D(j)-D(j+3)+B(j+3)*C(j)-C(j+3)*B(j)


def quadratic(n, d, j):
    A = slack_from_initial(n, d, j, 1, 0)
    C = slack_from_initial(n, d, j, 0, 1)
    B = slack_from_initial(n, d, j, 1, 1) - A - C
    return A, B, C


def multiplier_slack(c, n, j):
    C = lambda k: row_value(c, n, k)
    D = lambda k: C(k)**2 - C(k-1)*C(k+1)
    R = lambda k: C(k) + C(k-2)  # coefficients of (1+z^2)P
    DR = lambda k: R(k)**2 - R(k-1)*R(k+1)
    return DR(j+2)-DR(j+3) + 2*(D(j+1)-D(j+2))


def audit(amax):
    stats = Counter()
    minimum = None
    first_core_negative_W = None
    first_uncovered = None
    mismatch_recurrence = 0
    mismatch_multiplier = 0

    for a in range(amax+1):
        for e in range(amax+1):
            n = a+e
            c = row(a,e)
            for j in range((n+1)//2, n):
                val, W = slack_and_W(c,n,j)
                stats['pairs'] += 1
                stats['negative_slack'] += val < 0
                stats['W_negative'] += W < 0
                stats['W_nonnegative'] += W >= 0
                mismatch_multiplier += multiplier_slack(c,n,j) != val
                x = c[j]
                y = c[j+1] if j+1 <= n else 0
                mismatch_recurrence += (
                    slack_from_initial(n,a-e,j,x,y) != val
                )
                item = (val,a,e,j)
                if minimum is None or item < minimum:
                    minimum = item

                if a >= e+2 and e >= 3 and W < 0:
                    itemW = (n,a,e,j,W,val)
                    if (first_core_negative_W is None
                            or itemW[:4] < first_core_negative_W[:4]):
                        first_core_negative_W = itemW

                if W < 0:
                    stats['W_negative_rows_tested_by_form'] += 1
                    A,B,Cq = quadratic(n,a-e,j)
                    disc = B*B-4*A*Cq
                    if A > 0 and disc < 0:
                        stats['positive_definite_when_W_negative'] += 1
                    else:
                        stats['W_negative_not_PD'] += 1
                        itemU = (n,a,e,j,val,W,A,B,Cq,disc)
                        if (first_uncovered is None
                                or itemU[:4] < first_uncovered[:4]):
                            first_uncovered = itemU

    assert stats['negative_slack'] == 0
    assert mismatch_recurrence == 0
    assert mismatch_multiplier == 0
    print('box: 0 <= a,e <=', amax)
    print('row-index census:', dict(stats))
    print('minimum (E+,a,e,j):', minimum)
    print('first W<0 in folded open chamber (N,a,e,j,W,E+):',
          first_core_negative_W)
    print('first W<0 with non-PD quadratic '
          '(N,a,e,j,E+,W,A,B,C,disc):', first_uncovered)
    print('recurrence mismatches:', mismatch_recurrence)
    print('multiplier identity mismatches:', mismatch_multiplier)

    for a,e,j in ((4,6,5),(5,3,5),(19,3,20)):
        n = a+e
        c = row(a,e)
        val,W = slack_and_W(c,n,j)
        print('example (a,e,j), E+, W =', (a,e,j), val, W)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Exact q=2 plus row audit; writes no files.'
    )
    parser.add_argument('--amax', type=int, default=40,
                        help='scan 0 <= a,e <= amax (default: 40)')
    audit(parser.parse_args().amax)
