# E8 arXiv m100 finite bridge certificate

This note records the exact finite-bridge replay for the remaining E8 odd
window below the rectangular tail chain.

## Lemma 1 (moment source)

Let `m_k` be the E8 adjoint character moment.  Use `m_0,...,m_70` from
`character_ring_iter/logs/e8_70.log` and `m_71,...,m_100` from the
`adjoint_tensor_ranks.tex` ancillary file to Bourjaily-Plesser-Vergu,
arXiv:2412.21189, recorded locally in
`references/arxiv_2412_21189_e8_m71_m100.txt`.

Proof.  The checker
`character_ring_iter/direct_chain_e8_arxiv_m100_finite_bridge.py` reuses the
source audit from `direct_chain_rect_e8_n81_certificate.py`.  It checks
`m_0,...,m_70` against the local character-ring log, `m_0,...,m_30` against
`references/oeis_A179663.txt`, and `m_71,...,m_100` against the arXiv
ancillary transcription.  It returns

```text
moment_log_0_100_arxiv_extension: OK
```

## Lemma 2 (finite Q3 formula)

For odd `n <= 98`,

```text
Q_3^{E8}(n)
= 2 sum_{k=0}^n binom(n,k)
    (m_{k+2}m_{n-k} - m_{k+1}m_{n-k+1}).
```

Proof.  This is the same exact character-ring formula used by
`character_ring_iter/verify_chain.py`.  Since Lemma 1 supplies the moments up
to `m_100`, every term needed for `Q_3(n)` with odd `n <= 98` is present.

## Proposition 3 (E8 finite bridge)

The finite bridge values

```text
Q_3^{E8}(n),  n = 69, 71, 73, 75, 77,
```

are positive.  The direct Chain differences

```text
D_E8(n) = Q_3^{E8}(n+2) - 4 Q_3^{E8}(n),
  n = 67, 69, 71, 73, 75, 77,
```

are also positive.

Proof.  The checker evaluates the integers exactly and returns

```text
log10_min_bridge_q3=100.94
log10_min_chain_diff=100.94
bridge_q3_69_77_positive: OK
chain_diff_67_77_positive: OK
```

It also verifies the longer available exact Chain window

```text
D_E8(n) >= 0 for n = 67,69,...,95,
```

returning

```text
extended_chain_diff_67_95_positive: OK
```

## Corollary 4 (E8 finite bridge discharged)

The condition `E8-FiniteBridge(77)`, namely

```text
Q_3^{E8}(n) >= 0 for every odd 67 < n <= 77,
```

holds.

Proof.  This is the `bridge_q3_69_77_positive` conclusion of Proposition 3.
Together with `DCT_RECT_TRIG_E8_N77.md`, the direct Chain propagates from this
finite bridge to every odd `n >= 77`.
