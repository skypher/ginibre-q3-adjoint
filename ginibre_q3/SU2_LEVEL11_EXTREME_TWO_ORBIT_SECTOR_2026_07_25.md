# The complete extreme two-orbit sector of `SU(2)_11`

Date: 2026-07-25

## 1. Full finite theorem

At level eleven the odd simple-current lift has orbit fibres

```text
{V_2,V_9}   -> B_1,
{V_1,V_10}  -> B_5.
```

The rank-six orbit theorem in
`O11_B1_TOP_TWO_LABEL_THEOREM_2026_07_25.md` proves scalar `GKS2*` for every
signed word supported on `B_1,B_5`, with arbitrary exponents and signs.

Consider now any even-minus word in `SU(2)_11` supported on

```text
{V_1,V_2,V_9,V_10}.
```

If the number of odd-labelled occurrences (`V_1` and `V_9`) is odd, its
scalar corner is zero by the parity selection rule.  Otherwise the exact
odd-level lift is the average of two even-minus orbit words supported on
`B_1,B_5`; both are nonnegative by the complete two-label orbit theorem.
Thus:

> **Theorem 1.1.** Every even-minus signed word of arbitrary length in
> `SU(2)_11`, supported on `{V_1,V_2,V_9,V_10}`, has nonnegative scalar
> corner.

Adjoining one plus factor with any label in the same four-element set stays
inside the theorem.  The plus-factor identity therefore proves the partial
columns

```text
c in {1,2,9,10}
```

nonnegative for the same words.

## 2. Ordinary stable consequence

For an ordinary `SU(2)` word supported on `V_1,V_2`, with arbitrary signs,
the level-eleven theorem applies whenever the sharp stability threshold

```text
K(c)=max(2,c,ceil((T+c)/2))
```

is at most eleven.  Hence every scalar corner in this unrestricted-sign
sector is nonnegative when

```text
ceil(T/2)<=11,
```

and targets `c=1,2` are nonnegative whenever `K(c)<=11`.

This finite-window result complements the all-length theorem in
`SU2_V1_V2_MINUS_SECTOR_2026_07_25.md`, which drops the level-eleven bound
when every `V_2` factor has minus sign.
