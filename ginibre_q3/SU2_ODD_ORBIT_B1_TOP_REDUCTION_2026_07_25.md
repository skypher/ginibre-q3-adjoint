# A top-ray reduction for the first two-label odd-orbit chamber

Date: 2026-07-25

## 1. Universal identity

Let `O_(2m-1)` be the simple-current orbit ring at odd level, with basis
`B_0,...,B_(m-1)`, and put

```text
T=B_(m-1).
```

The truncated `SU(2)` fusion rule gives

```text
T^2=B_0+B_1.                                           (1.1)
```

In the doubled ring write

```text
D_a=B_a tensor B_0-B_0 tensor B_a,
S_a=B_a tensor B_0+B_0 tensor B_a.
```

If `x=T tensor B_0` and `y=B_0 tensor T`, then (1.1) gives the exact identity

```text
D_1=x^2-y^2=(x-y)(x+y)=D_T S_T.                       (1.2)
```

No spectral approximation or positivity estimate is involved.

## 2. Reduction to the top-orbit ray

For arbitrary nonnegative integers `r,a,b`,

```text
D_1^r D_T^a S_T^b = D_T^(r+a) S_T^(r+b).              (2.1)
```

If the original signed word has even minus parity, then `r+a` is even.  The
right side of (2.1) is therefore an even-minus word supported entirely on the
top orbit label.  The existing all-length odd top-orbit-ray theorem proves
its scalar corner nonnegative.

Consequently:

> **Lemma 2.1.**  In every odd simple-current orbit ring, every signed word
> supported on `D_1,D_T,S_T`, with an even total number of minus factors, has
> nonnegative scalar corner.

This is an unbounded statement in all three exponents.

## 3. The former first `O_11` obstruction

For `O_11`, `T=B_5`.  The chamber previously isolated by the transport fan is

```text
D_1^(2+2p) S_5^(1+2q),              p,q>=0.
```

Equation (2.1) converts it identically to

```text
D_5^(2+2p) S_5^(3+2p+2q),           p,q>=0.            (3.1)
```

The minus exponent in (3.1) is even, so the top-ray theorem proves the chamber
immediately.  The five-region weighted-AM-GM certificate remains a useful
independent exact replay, but it is not needed for this theorem.

More generally, every `O_11` chamber of the form

```text
D_1^r D_5^a S_5^b,                  r+a even,
```

is now removed from the unresolved frontier at once.

## 4. Significance

The identity explains why the first rank-six transport obstruction was
positive despite resisting coordinatewise Hall transport: it was not a new
rank-six phenomenon.  It was a disguised one-label top-ray word.  Future
frontier censuses should apply (1.2) before spectral transport or AM-GM
allocation, both to reduce the search and to avoid certifying the same
positivity twice.
