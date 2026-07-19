# F_4 primary-invariant structural moments

Computed exactly via Wick contraction on the $e^{-|w|^2}\prod(\alpha\cdot w)^2 dw$
measure in $\mathbb R^4$ with the 24 F_4 positive roots.

## Single-primary (verified in c3_rigorous_F4.py)

- $\langle I_6(F_4)\rangle_\mu = 0$ exactly
- $\langle I_6(F_4) Q\rangle_\mu = 0$ exactly
- $\langle I_6(F_4) Q^2\rangle_\mu = 0$ exactly
- $\langle I_6(F_4) Q^3\rangle_\mu = 0$ exactly
- $\langle I_8(F_4)\rangle_\mu = -485537/(4860 \cdot 20160/52^3) = -485537/(4860) / \text{coef}$... just $\langle S_{\text{primary}}\rangle_\mu = -485537/(4 \cdot $scale$)$
  with $I_8 := \sum(\alpha\cdot w)^8 - (357/32)|w|^8$

## Two-primary (new this pass)

- $\langle I_6(F_4) \cdot I_8(F_4)\rangle_\mu = 4144529025/2$ — **NONZERO**

This nonzero product-primary expectation enters $c_4(F_4)$ via terms like $G_1 G_3$ and $G_2^2/2$
in the Laplace-expansion formula, preventing a simple "all primary expectations vanish"
reduction for $c_4$.  Hence $c_4(F_4)$ requires explicit primary-invariant tracking.
