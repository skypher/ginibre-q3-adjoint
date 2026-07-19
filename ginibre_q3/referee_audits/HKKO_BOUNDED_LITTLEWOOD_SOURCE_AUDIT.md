# Direct source audit of the bounded-Littlewood determinant identities

Primary source: Huh--Kim--Krattenthaler--Okada, *Bounded Littlewood
identities with fixed number of odd rows or odd columns*, Electron. J.
Combin. 33 (2026), P3.5.  The archived author PDF
`references/bounded_littlewood_odd_rows_columns_2025.pdf` has SHA-256
`ae94e04a70c94131e35408ec0e938d1f38422a04f9a66ab4600f029dd8dd8bb2`.
The checks below were made by direct reading of PDF page 3, not from the
project's transcription or C++ implementation.

The source defines

```text
f_r(x) = sum_{n in Z} e_n(x)e_{n+r}(x),
e(x) = sum_{n>=0} e_n(x),
bar_e(x) = sum_{n>=0} (-1)^n e_n(x).
```

It also defines `c(lambda)` as the number of odd columns.

## Type B

Set `w=b>=1` and `u=0` in Theorem 1.3, equations (1.9) and (1.10).
Equation (1.9) becomes

```text
sum_{lambda_1<=2b+1, c(lambda)=0 or 2b+1} s_lambda
 = e det(f_{i-j}-f_{i+j-1}).
```

Equation (1.10) becomes the difference of those two sectors:

```text
sum_{lambda_1<=2b+1, c(lambda)=0} s_lambda
- sum_{lambda_1<=2b+1, c(lambda)=2b+1} s_lambda
 = bar_e det(f_{i-j}+f_{i+j-1}).
```

Their half-sum is exactly the paper's `mathcal B_b`.  In particular, both
the factor `1/2` and the use of `bar_e` in the plus determinant are required
and are present in `full_q3_extension.tex`.

## Type C

Equation (1.6), with `w=b`, is

```text
sum_{lambda_1<=2b, r(lambda)=0} s_lambda
 = det(f_{i-j}-f_{i+j}).
```

The paper applies the Hall-isometric involution `omega`, which sends `e_r` to
`h_r`, `e_2` to `h_2`, and `s_lambda` to `s_{lambda'}`.  This changes `f^e`
to `f^h` and gives exactly the displayed `mathcal C_b`, with no additional
factor.

## Type D

Set `w=b>=1` and `u=0` in equation (1.7).  The two surviving sectors are
`c(lambda)=0` and `c(lambda)=2b`, while only the `k=0` determinant survives
on the right.  The source formula already contains the prefactor `1/2`, so

```text
sum_{lambda_1<=2b, c(lambda)=0 or 2b} s_lambda
 = (1/2) det(f_{i-j}+f_{i+j-2}),
```

which is exactly the paper's `mathcal D_b`.  No second factor of two is
missing.

Thus all three determinant identities in Lemma
`lem:bounded-littlewood-classical-moments` agree with the cited primary
source, including their width bounds, signs, index shifts, alternating
elementary factor, and normalization.
