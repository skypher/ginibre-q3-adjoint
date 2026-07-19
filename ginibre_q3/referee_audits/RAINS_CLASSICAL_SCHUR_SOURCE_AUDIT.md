# Direct source audit of the classical Schur integral criteria

Primary source: E. M. Rains, *Increasing subsequences and the classical
groups*, Electron. J. Combin. 5 (1998), R12.  The archived PDF
`references/rains_1998_increasing_subsequences.pdf` has SHA-256
`d5242a4d38f52c7af73ed5ce65e83dc3cc8468f9f37636101ffcbc721f80cdc3`.
This audit uses the paragraph immediately preceding Lemma 3.1 on journal
pages 5--6.

Rains states the following two normalized Haar integrals for a partition
`lambda`:

- over `O(N)`, the integral of `s_lambda` is one exactly when `lambda` has
  at most `N` parts and every part is even, and is zero otherwise;
- over `Sp(2b)`, it is one exactly when `lambda` has at most `2b` parts and
  every part size occurs with even multiplicity, equivalently every part of
  `lambda'` is even, and is zero otherwise.

These are precisely the two criteria used in Lemma
`lem:finite-orthogonal-adjoint-moments`.

For type `B_b`, the adjoint character is the character of `wedge^2 V`.
The central map `g -> -g` interchanges the two components of `O(2b+1)` and
acts trivially on `wedge^2 V`; hence its tensor-power moments over `SO(2b+1)`
equal the `O(2b+1)` moments.  Rains' first criterion admits exactly the even
partitions `2 nu` with `ell(nu)<=2b+1`, giving the paper's stable sum minus
the length-wall complement.

For type `C_b`, the adjoint character is the character of `Sym^2 V`.
Rains' second criterion admits precisely

```text
lambda=(nu_1,nu_1,nu_2,nu_2,...),  ell(nu)<=b.
```

Combined with the symmetric-function involution interchanging `h_2` and
`e_2`, this is exactly the paper's type-C formula.

For type `D_b`, normalized component averages satisfy

```text
integral_SO(2b) f = integral_O(2b) f + integral_O(2b) det(g) f(g).
```

The first term is Rains' even-part sector.  In the second term,
`det*s_lambda` is the Schur character obtained by adding one full column of
height `2b`.  Applying Rains' even-part criterion after that determinant
twist shows that the original `lambda` must have exactly `2b` positive odd
parts.  Removing one full column therefore gives uniquely

```text
lambda=(1+2nu_1,...,1+2nu_(2b)),  |nu|=j-b,
```

which is the correction term in the paper.  Thus the source criterion and
the elementary component decomposition give all three finite classical
moment formulas with the stated signs and without an omitted factor of two.
