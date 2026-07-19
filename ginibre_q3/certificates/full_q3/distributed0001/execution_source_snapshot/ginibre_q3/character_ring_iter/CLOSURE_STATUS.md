# Closure status: Q_3 non-abelian, pre-asymptotic residues

## Per-family summary

| Family | $c_1$ | $c_2$ | $c_3$ | Diagnostic threshold $(8|c_3|)^{1/3}$ | $n_{\max}(G)$ | Live closure |
|---|---|---|---|---|---|---|
| $G_2$ | $-63/4$ (rigorous) | $165851/1152$ (Wick rigorous) | $-3080875/3072$ (rigorous) | $20.03$ | $197$ | DCT-RectTrig tail certificate for odd `n >= 21` |
| $F_4$ | $-182$ (rigorous) | $8401705/486$ (Wick rigorous) | $-44918139383/39366$ (rigorous) | $208.08$ | $217$ | DCT-RectTrig tail certificate for odd `n >= 65` |
| $E_6$ | $-1599/4$ (rigorous) | $10524215/128$ (Dunkl exact) | $-11872316417/1024$ (Dunkl exact) | $452.66$ | $77$ | DCT-RectTrig tail certificate for odd `n >= 75` |
| $E_7$ | $-18221/16$ (rigorous) | $82042395295/124416$ (Dunkl exact) | $-41730554417849633/161243136$ (Dunkl exact) | $1274.54$ | $67$ | DCT-RectTrig tail certificate for odd `n >= 65` |
| $E_8$ | $-3906$ (rigorous) | $346377539/45$ (rank-$6$ isotropy rigorous) | $-11478332145838/1125$ (Dunkl exact) | $4337.83$ | $67$ | DCT-RectTrig tail certificate for odd `n >= 133`, wide-cap bridge steps `n=131,129,127,125,123,121,119,117,115,113,111,109,107,105,103,101,99,97,95,93,91,89,87,85,83,81,79,77`, and arXiv m100 finite-bridge replay |

## What's done this pass

### Rigorous
0. `direct_chain_rect_g2_f4_certificate.cpp` closes the G2 and F4 direct
   Chain tails by exact GMP rectangular certificates from odd `n >= 21` and
   `n >= 65`.  It verifies the non-simply-laced short-root normalization,
   quartic root identities, base domination, and all-later odd-step
   propagation.  The numbered proof note is `../DCT_RECT_TRIG_G2_F4.md`.
1. `c_2(G_2), c_2(F_4), c_2(E_8)` exact rational (Wick / rank-6 isotropy).
2. `c_3(G_2) = -3080875/3072` via polynomial identity $\prod(\alpha \cdot w)^2 = (27/1024)(|w|^{12} - I_6^{(0)2})$ + angular vanishing.
3. `c_3(F_4) = -44918139383/39366` via Wick verification that `<I_6(F_4) Q^k>_mu = 0` for $k=0,1,2,3$, plus exact computation of rank-$8$ primary contribution.
4. Character-ring verification that empirical `c_3(G_2)` matches the closed form to 6 decimals at $n=200$.
5. Direct integer-arithmetic verification of shift-2 log-convexity across 554 cases (both parities, all 5 families).
6. Exact Chain replay from `logs/f4_220c.log`: F_4 verified through
   Chain `m=107`, hence odd `n=217`, with OEIS overlap `m_0..m_20`
   matching and last difference
   `1.36 × 10^371`; this exceeds the reduced threshold `208.08`.
7. Exact Chain replay from `logs/e6_80.log`, `logs/e7_70.log`,
   `logs/e8_70.log`: E_6 verified through Chain `m=37`, E_7 through
   `m=32`, and E_8 through `m=32`, each with OEIS overlap matching.
8. The legacy independent `dunkl_exceptional_coefficients.py` computes the remaining exceptional
   coefficient inputs by the Dunkl heat-kernel identity. It first reproduces
   the existing F_4 Wick normalizations and exact coefficients:
   `c_2(F_4) = 8401705/486` and
   `c_3(F_4) = -44918139383/39366`.
9. The same exact computation gives
   `c_2(E_6) = 10524215/128`,
   `c_3(E_6) = -11872316417/1024`,
   `c_2(E_7) = 82042395295/124416`, and
   `c_3(E_7) = -41730554417849633/161243136`.
10. The E_8 degree-$8$ slot gives
    `<S_{E_8}(w)> = 3667734000`, `<S_{E_8}(u)> = 67921/15`,
    `<S_{\mathrm{primary}}(u)> = -651/10`, and
    `c_3(E_8) = -11478332145838/1125`.
11. `direct_chain_tail_constants.py` records the exact constants for the
    direct Chain-tail reduction
    `D_G(n) = Q_3^G(n+2) - 4 Q_3^G(n)`.
12. `direct_chain_a0_dominance.py` computes the E_6/E_7/E_8 leading
    positive-saddle coefficient `A_0(G)` from the Macdonald-Mehta Gaussian
    integral and checks the conditional one-term half-remainder dominance
    margin.
13. `direct_chain_halflead_audit.py` compares the half-leading lower bound
    with `Q_3(n) + 2(2C_G)^n` throughout the exact logs and rules out the
    half-leading estimate on every currently available odd E_6/E_7/E_8
    Chain-window value.
14. `direct_chain_asymptotic_ratio_audit.py` compares the exact logged
    `Q_3` and direct Chain difference scales against the Macdonald-Mehta
    leading coefficient `A_0(G)`.
15. `direct_chain_truncated_mass.py` computes the finite-n truncated
    Gamma-Wronskian radial mass available under root-angle caps
    `|alpha.v| <= delta`, for both `Q_3` and the direct Chain difference
    `D_G(n)=Q_3(n+2)-4Q_3(n)`.
16. `direct_chain_cap_thresholds.py` inverts that truncated-mass comparison
    at the current E_6/E_7/E_8 prefixes and solves the cap size needed by
    the identity-saddle model to match the exact logged `Q_3` and `D_G`
    coefficients.
17. `direct_chain_rect_lower_model.py` tests a proof-shaped rectangular
    radial wedge for the wide-region trigonometric lower-bound route.  The
    E_6 prefix has positive margin in this conservative model; E_7 and E_8
    need a stronger positive-region estimate.
18. `direct_chain_rect_e6_e7_rank_certificate.cpp` replaces both former
    single-rectangle scripts.  Those scripts used the false character
    endpoints `C_E6=2` and `C_E7=3` and are deleted.  The new OpenMP/GMP
    checker uses the source-audited universal bound `C_G=rank(G)`, exact
    root-generated quartic identities, wide-cap scalar Bernstein checks,
    and rational rectangle unions.  It closes E6 from odd `n>=75` with
    `23088` cells and normalized base log-margin `8.904`.  For E7 it combines `66925`
    cells with a degree-26 Chebyshev negative-region majorant from the exact
    moments `m_0,...,m_70`, giving normalized base log-margin `0.426` and a cellwise
    odd-step ratio above `101/100`.  The numbered notes are
    `../DCT_RECT_TRIG_E6.md` and `../DCT_RECT_TRIG_E7.md`.
19. Garibaldi--Guralnick--Rains, Theorem 2.1, is the active source for
    `Tr Ad(g)>=-rank(G)`.  Their Theorem 2.3/Table 1 gives the sharp E6
    minimum `-3`, explicitly exposing the former `-2` input as invalid.
20. `direct_chain_quartic_rect_model.py` applies the same quartic-improved
    rectangle-union diagnostic to E_7 and E_8.  With an `80 x 80` grid it
    gives margin `+6.88` for E_7 and `-68.51` for E_8.
21. `direct_chain_rect_e8_tail_certificate.py` replaces the negative current
    prefix diagnostic by an exact-rational higher-start E_8 tail certificate.
    The positive region is a `1/40 x 1/40` rectangular union with `13558`
    retained cells, and the numbered proof note is
    `../DCT_RECT_TRIG_E8_TAIL.md`.  It proves the direct Chain tail
    comparison for every odd `n >= 133`; the finite E_8 bridge below that
    start is handled by the bridge certificates and the arXiv m100
    finite-bridge replay below.
22. `direct_chain_rect_e8_n131_certificate.py` closes the next E_8 bridge
    step `D_E8(131) = Q_3^{E8}(133)-4Q_3^{E8}(131)` by an exact-rational
    wide-cap rectangular-union certificate.  It uses the E_8 quartic and
    sextic identities, the cap `|alpha.v| <= 19/6`, and the numbered proof
    note `../DCT_RECT_TRIG_E8_N131.md`.
23. `direct_chain_rect_e8_n129_certificate.py` closes the next E_8 bridge
    step `D_E8(129) = Q_3^{E8}(131)-4Q_3^{E8}(129)` by an exact-rational
    wide-cap rectangular-union certificate.  It uses the same E_8 quartic
    and sextic identities, the cap `|alpha.v| <= 19/6`, the sharper
    rational exponential lower bound from `e < 1457/536`, and the numbered
    proof note `../DCT_RECT_TRIG_E8_N129.md`.
24. `direct_chain_rect_e8_n127_certificate.py` closes the next E_8 bridge
    step `D_E8(127) = Q_3^{E8}(129)-4Q_3^{E8}(127)` by an exact-rational
    wide-cap rectangular-union certificate.  It widens the cap to
    `|alpha.v| <= 10/3`, weakens the gap coefficient to `1/450`, keeps the
    sine denominator `21`, and uses the numbered proof note
    `../DCT_RECT_TRIG_E8_N127.md`.
25. `direct_chain_rect_e8_n125_certificate.py` closes the next E_8 bridge
    step `D_E8(125) = Q_3^{E8}(127)-4Q_3^{E8}(125)` by an exact-rational
    wide-cap rectangular-union certificate.  It widens the cap to
    `|alpha.v| <= 7/2`, keeps the gap coefficient `1/450` and sine
    denominator `21`, and uses the numbered proof note
    `../DCT_RECT_TRIG_E8_N125.md`.
26. `direct_chain_rect_e8_n123_certificate.py` closes the next E_8 bridge
    step `D_E8(123) = Q_3^{E8}(125)-4Q_3^{E8}(123)` by an exact-rational
    wide-cap rectangular-union certificate.  It widens the cap to
    `|alpha.v| <= 4`, uses the gap coefficient `1/480` and sine denominator
    `20`, and uses the numbered proof note `../DCT_RECT_TRIG_E8_N123.md`.
27. `direct_chain_rect_e8_n121_certificate.py` closes the next E_8 bridge
    step `D_E8(121) = Q_3^{E8}(123)-4Q_3^{E8}(121)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the cap
    `|alpha.v| <= 4`, uses the gap coefficient `1/475`, the eighth-moment
    Cauchy lower bound `sum (alpha.u)^8 >= Q(u)^4/300000`, and the numbered
    proof note `../DCT_RECT_TRIG_E8_N121.md`.
28. `direct_chain_rect_e8_n119_certificate.py` closes the next E_8 bridge
    step `D_E8(119) = Q_3^{E8}(121)-4Q_3^{E8}(119)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the cap
    `|alpha.v| <= 4`, adds the matching upper bound
    `2(1-cos z) <= z^2 - z^4/12 + z^6/360` to the gap estimate, and uses
    the numbered proof note `../DCT_RECT_TRIG_E8_N119.md`.
29. `direct_chain_rect_e8_n117_certificate.py` closes the next E_8 bridge
    step `D_E8(117) = Q_3^{E8}(119)-4Q_3^{E8}(117)` by an exact-rational
    wide-cap rectangular-union certificate.  It replaces the Taylor
    coefficients in the gap and average-character estimates by rational
    cubic scalar polynomials certified by exact Sturm root counts, and uses
    the numbered proof note `../DCT_RECT_TRIG_E8_N117.md`.
30. `direct_chain_rect_e8_n115_certificate.py` closes the next E_8 bridge
    step `D_E8(115) = Q_3^{E8}(117)-4Q_3^{E8}(115)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the cubic scalar
    polynomials and improves the Weyl-sine density loss to
    `log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2100`, certified by an exact
    Sturm root count; the numbered proof note is
    `../DCT_RECT_TRIG_E8_N115.md`.
31. `direct_chain_rect_e8_n113_certificate.py` closes the next E_8 bridge
    step `D_E8(113) = Q_3^{E8}(115)-4Q_3^{E8}(113)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the cubic scalar
    polynomials and quartic Weyl-sine loss, and replaces the universal
    negative-region bound by the exact low-moment bound
    `176 * 16^(n-2)`; the numbered proof note is
    `../DCT_RECT_TRIG_E8_N113.md`.
32. `direct_chain_rect_e8_n111_certificate.py` closes the next E_8 bridge
    step `D_E8(111) = Q_3^{E8}(113)-4Q_3^{E8}(111)` by an exact-rational
    wide-cap rectangular-union certificate.  It reuses the `n=113`
    scalar, sine, and low-moment negative machinery at `N=111`; the
    numbered proof note is `../DCT_RECT_TRIG_E8_N111.md`.
33. `direct_chain_rect_e8_n109_certificate.py` closes the next E_8 bridge
    step `D_E8(109) = Q_3^{E8}(111)-4Q_3^{E8}(109)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the same positive
    rectangular machinery and uses the `p=12` exact moment negative bound
    `1092455198188 * 16^(n-12)`; the numbered proof note is
    `../DCT_RECT_TRIG_E8_N109.md`.
34. `direct_chain_rect_e8_n107_certificate.py` closes the next E_8 bridge
    step `D_E8(107) = Q_3^{E8}(109)-4Q_3^{E8}(107)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the same positive
    rectangular machinery and replaces the negative-region estimate by the
    asymmetric degree-28 moment majorant
    `(63/65)*16^95*S^12*(S^2+4)*((496-S)/512)^14`; the numbered proof note
    is `../DCT_RECT_TRIG_E8_N107.md`.
35. `direct_chain_rect_e8_n105_certificate.py` closes the next E_8 bridge
    step `D_E8(105) = Q_3^{E8}(107)-4Q_3^{E8}(105)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the same positive
    rectangular machinery and uses the Chebyshev degree-9 negative-region
    majorant `(63/65)*16^97*S^8*(S^2+4)*C_9(S)^2`, where
    `C_9(S)=T_9(S/248-1)/T_9(-33/31)`; the numbered proof note is
    `../DCT_RECT_TRIG_E8_N105.md`.
36. `direct_chain_rect_e8_n103_certificate.py` closes the next E_8 bridge
    step `D_E8(103) = Q_3^{E8}(105)-4Q_3^{E8}(103)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the same positive
    rectangular machinery and uses the shifted Chebyshev degree-14
    negative-region majorant `(63/65)*16^95*S^8*(S^2+4)*C_14(S)^2`,
    where `C_14(S)=T_14((S-249)/247)/T_14(-265/247)`; the numbered proof
    note is `../DCT_RECT_TRIG_E8_N103.md`.
37. `direct_chain_rect_e8_n101_certificate.py` closes the next E_8 bridge
    step `D_E8(101) = Q_3^{E8}(103)-4Q_3^{E8}(101)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the same positive
    rectangular machinery and uses the shifted Chebyshev degree-18
    negative-region majorant `(63/65)*16^93*S^8*(S^2+4)*C_18(S)^2`,
    where `C_18(S)=T_18((S-249)/247)/T_18(-265/247)`.  The interval
    residuals are certified by exact Bernstein-coefficient positivity; the
    numbered proof note is `../DCT_RECT_TRIG_E8_N101.md`.
38. `direct_chain_rect_e8_n99_n97_certificate.py` closes the next two E_8
    bridge steps `D_E8(99) = Q_3^{E8}(101)-4Q_3^{E8}(99)` and
    `D_E8(97) = Q_3^{E8}(99)-4Q_3^{E8}(97)` by one exact-rational wide-cap
    rectangular-union certificate.  It keeps the same positive rectangular
    machinery and uses the shared shifted Chebyshev degree-26
    negative-region majorant `(63/65)*16^(n-8)*S^8*(S^2+4)*C_26(S)^2`,
    where `C_26(S)=T_26((S-249)/247)/T_26(-265/247)`.  The interval
    residuals are certified by exact Bernstein-coefficient positivity; the
    numbered proof note is `../DCT_RECT_TRIG_E8_N99_N97.md`.
39. `direct_chain_rect_e8_n95_certificate.py` closes the next E_8 bridge
    step `D_E8(95) = Q_3^{E8}(97)-4Q_3^{E8}(95)` by an exact-rational
    wide-cap rectangular-union certificate.  It keeps the same positive
    rectangular machinery and uses the shifted Chebyshev degree-27
    negative-region majorant `(63/65)*16^83*S^12*(S^2+4)*C_27(S)^2`,
    where `C_27(S)=T_27((2S-503)/489)/T_27(-535/489)`.  The interval
    residuals are certified by exact Bernstein-coefficient positivity; the
    numbered proof note is `../DCT_RECT_TRIG_E8_N95.md`.
40. `direct_chain_rect_e8_n93_certificate.py` closes the next E_8 bridge
    step `D_E8(93) = Q_3^{E8}(95)-4Q_3^{E8}(93)` by an exact-rational
    widened-cap rectangular-union certificate with root-angle cap `21/5`.
    It uses a widened scalar check, the sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2000`, and the same shifted
    Chebyshev degree-27 negative-region majorant
    `(63/65)*16^81*S^12*(S^2+4)*C_27(S)^2`; the numbered proof note is
    `../DCT_RECT_TRIG_E8_N93.md`.
41. `direct_chain_rect_e8_n91_certificate.py` closes the next E_8 bridge
    step `D_E8(91) = Q_3^{E8}(93)-4Q_3^{E8}(91)` by an exact-rational
    rectangular-union certificate with root-angle cap `13/3`.  It uses the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/115000` and an optimized
    degree-27 square negative-region majorant
    `16^79*S^12*(S^2+4)*P_27(S)^2`; the interval residuals are certified by
    exact Bernstein-coefficient positivity, and the numbered proof note is
    `../DCT_RECT_TRIG_E8_N91.md`.
42. `direct_chain_rect_e8_n89_certificate.py` closes the next E_8 bridge
    step `D_E8(89) = Q_3^{E8}(91)-4Q_3^{E8}(89)` by an exact-rational
    rectangular certificate.  It uses root-angle cap `14/3`, grid step
    `1/20`, the cap-adjusted gap coefficient `169078/100000000`, the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and a degree-28
    square negative-region majorant.  The two new moments `m_71,m_72` are
    sourced from the arXiv 2412.21189 ancillary `adjoint_tensor_ranks.tex`
    `e[8]` block in `../references/arxiv_2412_21189_e8_m71_m72.txt`, with
    local overlap through `m_70` checked against `logs/e8_70.log`.  The
    numbered proof note is `../DCT_RECT_TRIG_E8_N89.md`.
43. `direct_chain_rect_e8_n87_certificate.py` closes the next E_8 bridge
    step `D_E8(87) = Q_3^{E8}(89)-4Q_3^{E8}(87)` by an exact-rational
    rectangular certificate.  It keeps root-angle cap `14/3`, grid step
    `1/20`, the cap-adjusted gap coefficient `169078/100000000`, and the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, while raising the
    square negative-region majorant to degree `30`.  The new moments
    `m_71,...,m_76` are sourced from the arXiv 2412.21189 ancillary
    `adjoint_tensor_ranks.tex` `e[8]` block in
    `../references/arxiv_2412_21189_e8_m71_m76.txt`, with local overlap
    through `m_70` checked against `logs/e8_70.log`.  The numbered proof note
    is `../DCT_RECT_TRIG_E8_N87.md`.
44. `direct_chain_rect_e8_n85_certificate.py` closes the next E_8 bridge
    step `D_E8(85) = Q_3^{E8}(87)-4Q_3^{E8}(85)` by an exact-rational
    rectangular certificate.  It keeps root-angle cap `14/3`, grid step
    `1/20`, the cap-adjusted gap coefficient `169078/100000000`, and the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, while raising the
    square negative-region majorant to degree `34`.  The new moments
    `m_71,...,m_84` are sourced from the arXiv 2412.21189 ancillary
    `adjoint_tensor_ranks.tex` `e[8]` block in
    `../references/arxiv_2412_21189_e8_m71_m84.txt`, with local overlap
    through `m_70` checked against `logs/e8_70.log`.  The numbered proof note
    is `../DCT_RECT_TRIG_E8_N85.md`.
45. `direct_chain_rect_e8_n83_certificate.py` closes the next E_8 bridge
    step `D_E8(83) = Q_3^{E8}(85)-4Q_3^{E8}(83)` by an exact-rational
    rectangular certificate.  It keeps root-angle cap `14/3`, grid step
    `1/20`, the cap-adjusted gap coefficient `169078/100000000`, and the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, while raising the
    square negative-region majorant to degree `38`.  The new moments
    `m_71,...,m_92` are sourced from the arXiv 2412.21189 ancillary
    `adjoint_tensor_ranks.tex` `e[8]` block in
    `../references/arxiv_2412_21189_e8_m71_m92.txt`, with local overlap
    through `m_70` checked against `logs/e8_70.log`.  The numbered proof note
    is `../DCT_RECT_TRIG_E8_N83.md`.
46. `direct_chain_rect_e8_n81_certificate.py` closes the next E_8 bridge
    step `D_E8(81) = Q_3^{E8}(83)-4Q_3^{E8}(81)` by an exact-rational
    rectangular certificate.  It keeps root-angle cap `14/3`, grid step
    `1/20`, the cap-adjusted gap coefficient `169078/100000000`, and the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, while raising the
    square negative-region majorant to degree `42`.  The new moments
    `m_71,...,m_100` are sourced from the arXiv 2412.21189 ancillary
    `adjoint_tensor_ranks.tex` `e[8]` block in
    `../references/arxiv_2412_21189_e8_m71_m100.txt`, with local overlap
    through `m_70` checked against `logs/e8_70.log`.  The numbered proof note
    is `../DCT_RECT_TRIG_E8_N81.md`.
47. `direct_chain_rect_e8_n79_certificate.py` closes the next E_8 bridge
    step `D_E8(79) = Q_3^{E8}(81)-4Q_3^{E8}(79)` by an exact-rational
    rectangular certificate.  It keeps root-angle cap `14/3`, grid step
    `1/20`, the cap-adjusted gap coefficient `169078/100000000`, and the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, while using a
    degree-37 square negative-region majorant at moment weight `p=22`.
    The new moments `m_71,...,m_100` are sourced from the arXiv 2412.21189
    ancillary `adjoint_tensor_ranks.tex` `e[8]` block in
    `../references/arxiv_2412_21189_e8_m71_m100.txt`, with local overlap
    through `m_70` checked against `logs/e8_70.log`.  The numbered proof note
    is `../DCT_RECT_TRIG_E8_N79.md`.
48. `direct_chain_rect_e8_n77_certificate.py` closes the next E_8 bridge
    step `D_E8(77) = Q_3^{E8}(79)-4Q_3^{E8}(77)` by an exact-rational
    rectangular certificate.  It keeps root-angle cap `14/3`, grid step
    `1/20`, the cap-adjusted gap coefficient `169078/100000000`, and the
    sextic sine bound
    `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, while using a
    degree-37 square negative-region majorant at moment weight `p=22`.
    The new moments `m_71,...,m_100` are sourced from the arXiv 2412.21189
    ancillary `adjoint_tensor_ranks.tex` `e[8]` block in
    `../references/arxiv_2412_21189_e8_m71_m100.txt`, with local overlap
    through `m_70` checked against `logs/e8_70.log`.  The numbered proof note
    is `../DCT_RECT_TRIG_E8_N77.md`.
49. `direct_chain_e8_arxiv_m100_finite_bridge.py` discharges the remaining
    E_8 finite bridge by exact integer evaluation of the source-checked
    `m_0,...,m_100` moment table.  It verifies
    `Q_3^{E8}(69),...,Q_3^{E8}(77) > 0`, the bridge differences
    `D_E8(67),D_E8(69),...,D_E8(77) > 0`, and the longer available exact
    Chain window `D_E8(67),D_E8(69),...,D_E8(95) > 0`.  The numbered proof
    note is `../E8_ARXIV_M100_FINITE_BRIDGE.md`.
50. `verify_e8_negative_bounds_gmp.cpp` closes the replay-input gap in the unified
    E8 rectangle checker.  The previously unmanifested
    `e8_rect_negative_bounds.tsv` is now a manifest-pinned certificate output.
    Before the OpenMP/GMP rectangle sum consumes it, the audit checks its exact
    29-row schema, reconstructs every nontrivial pointwise majorant from a
    manifest-pinned rational factor witness, verifies the full `m_0..m_100`
    source chain, integrates the witness with GMP, rederives the ten crude
    rank-bound rows directly, and requires exact rational equality row by row.
    Exact Bernstein subdivision handles the only two rows that need interval
    splitting (`n=81` at depth 3 and `n=85` at depth 5).  No historical Python
    row module is active in this audit.

### Historical diagnostics, now superseded
The earlier MC/fit estimates for `c_2(E_6)` and `c_2(E_7)` are no longer
used in the proof ledger; the exact Dunkl values above replace them.

## Bypassed or archived routes

The source-audited theorem route below no longer depends on these alternative
routes after the G_2/F_4/E_6/E_7 rectangular certificates, the E_8
tail-and-bridge certificates, and the arXiv m100 finite-bridge replay.  They
are retained as diagnostic routes and cross-checks.  In particular, the
displayed coefficient thresholds alone never constituted a remainder theorem.

### Wong `K_{univ}(G)` / Borel-radius bound
Bottleneck: explicit bound on nearest complex singularity of
`-log(chi_adj/d)` in `T_C`. Symbolic work, no compute needed. This route
would convert the empirical `|R_3(n)| <= |c_3|/n^3` estimate into an explicit
bound `<= 2|c_3|/n^3` for `n >= 2 C_Borel(G)`.

### Character-ring extension thresholds
The exact-coefficient reduced thresholds now read:
- E_6: `452.66`, so the current reduced-threshold theorem asks for odd
  `n >= 453`, equivalently Chain `m >= 225`; this is now bypassed by the
  `DCT-RectTrig(E_6)` direct-tail certificate from odd `n >= 75`.
- E_7: `1274.54`, so it asks for odd `n >= 1275`, equivalently
  Chain `m >= 636`; this is now bypassed by the `DCT-RectTrig(E_7)`
  direct-tail certificate from odd `n >= 65`.
- E_8: `4337.83`, so it asks for odd `n >= 4339`, equivalently
  Chain `m >= 2168`.

The E_8 finite bridge between the exact prefix and the sharpened
tail-and-bridge start is discharged by
`direct_chain_e8_arxiv_m100_finite_bridge.py`.

### Direct Chain-tail reduction diagnostics

For odd `n`, set `D_G(n) = Q_3^G(n+2) - 4 Q_3^G(n)`.  If the positive
Cartan saddle has expansion

```text
Q_{3,+}^G(n) = (2d)^n n^{-alpha} (A_0(G) + A_1(G)/n + O_G(n^{-2})),
```

then direct expansion gives

```text
D_{G,+}(n) = (2d)^n n^{-alpha}
             (B_0(G) + B_1(G)/n + O_G(n^{-2})),
B_0(G) = 4(d^2 - 1) A_0(G),
B_1(G) = -8 alpha d^2 A_0(G) + 4(d^2 - 1) A_1(G).
```

The negative-region input already in `NONABELIAN.md` gives
`|I_-(n)| <= 2(2C_G)^n`, hence

```text
|D_{G,-}(n)| <= 2((2C_G)^2 + 4)(2C_G)^n.
```

Exact constants:

| Group | `alpha` | `4(d^2-1)` | `-8 alpha d^2` | negative factor | `d/C_G` | current Chain `n_max` |
|---|---:|---:|---:|---:|---:|---:|
| `G_2` | 16 | 780 | -25088 | 40 | 7 | 197 |
| `F_4` | 54 | 10812 | -1168128 | 136 | 13 | 217 |
| `E_6` | 80 | 24332 | -3893760 | 40 | 39 | 77 |
| `E_7` | 135 | 70752 | -19104120 | 80 | 133/3 | 67 |
| `E_8` | 250 | 246012 | -123008000 | 520 | 31 | 67 |

Intermediate target `DCT-Rem(G)`: compute the positive-saddle constants and
give an explicit positive-saddle remainder bound for `D_{G,+}(n)` which,
together with the displayed negative-region bound, implies `D_G(n) >= 0`
for every odd `n` above the current exact Chain prefix.

Conditional subleaf `DCT-A0(G)` for `G in {E_6,E_7,E_8}`: with the
simply-laced normalization `|alpha|^2 = 2`, put `kappa_G = 2|R_+|/r`.
The checker uses the Macdonald-Mehta integral at exponent `1`,

```text
int exp(-|x|^2/2) prod_{alpha>0}(alpha.x)^2 dx
  = (2 pi)^(r/2) prod_i d_i!,
```

and `|W| = prod_i d_i` to compute

```text
A_0(G) =
  8 a_G d^2 (d/kappa_G)^d prod_i((d_i-1)!)^2 / (2 pi)^r,
where a_G = |R_+| + r/2 = d/2.
```

Under the one-term half-leading input at both `n` and `n+2`,
`|E_G(n)| <= A_0(G)/2` and `|E_G(n+2)| <= A_0(G)/2`, direct subtraction gives

```text
D_{G,+}(n) >= (2d)^n n^(-alpha) A_0(G)
             * (2 d^2 (n/(n+2))^alpha - 6).
```

Using the rational lower bound `pi < 4`, the dominance margins over the
negative-region bound at the current exact prefixes are:

| Group | tail start `n` | `log10(A_0 lower)` | `log10` dominance margin |
|---|---:|---:|---:|
| `E_6` | 77 | 103.00 | 76.19 |
| `E_7` | 67 | 202.75 | 67.48 |
| `E_8` | 67 | 440.24 | 82.79 |

This proves only the conditional dominance step: if a half-leading estimate
is available at a given `n`, the bad-region contribution is dominated.

Exact-log audit.  For odd `n`,
`Q_{3,+}(n) = Q_3(n) + |I_-(n)| <= Q_3(n) + 2(2C_G)^n`.  The current tail
starts give:

| Group | prefix `n` | `log10(Q_3(n)+2(2C_G)^n)` | `log10` half-leading lower | excess |
|---|---:|---:|---:|---:|
| `E_6` | 77 | 113.99 | 120.65 | 6.66 |
| `E_7` | 67 | 97.13 | 118.39 | 21.26 |
| `E_8` | 67 | 97.22 | 164.02 | 66.80 |

Scanning every odd `n` with exact logged moments gives:

| Group | tested odd `n` | last impossible `n` | first not ruled out | minimum excess |
|---|---:|---:|---:|---:|
| `E_6` | `1..77` | 77 | none | 6.66 |
| `E_7` | `1..67` | 67 | none | 21.26 |
| `E_8` | `1..67` | 67 | none | 66.80 |

Therefore `DCT-HalfLead(G)` cannot start anywhere inside the current exact
E_6/E_7/E_8 logs.

Asymptotic-scale audit.  Using `pi > 3` to display an upper estimate for
`A_0(G)`, the exact logs give:

| Group | `n` for `Q_3` | `log10 A_0 upper` | `log10 Q_3` coefficient | gap | `n` for `D_G` | `log10 D_G` coefficient | gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| `E_6` | 77 | 103.75 | 96.04 | 7.71 | 75 | 99.51 | 7.41 |
| `E_7` | 67 | 203.62 | 181.19 | 22.44 | 65 | 184.26 | 22.13 |
| `E_8` | 67 | 441.24 | 373.14 | 68.10 | 65 | 375.24 | 67.76 |

Thus the obstruction is late effective onset of the identity saddle at the
logged `n`, not merely the single prefix comparison.  The next useful
finite-n theorem should bound the actual local saddle mass at these `n`
instead of using a fixed fraction of the full Macdonald-Mehta integral.

Truncated local-mass diagnostic.  For a simply-laced exceptional group,
`Q(v)=kappa_G |v|^2`.  Under the rescaling `v=u/sqrt(n)`, the cap
`|alpha.v| <= delta` is guaranteed by

```text
s = Q(u)/(2d) <= delta^2 kappa_G n/(4d).
```

The radial part of the Gamma-Wronskian leading coefficient inside this cap
has fraction

```text
F_a(S) = (a+1) P(a,S) P(a+2,S) - a P(a+1,S)^2,
```

where `a=d/2`, `S=delta^2 kappa_G n/(4d)`, and `P(a,S)` is the regularized
lower incomplete gamma function.  For the direct Chain difference the
truncated model uses

```text
A_0(G) * [4d^2 (n/(n+2))^alpha F_a(S_{n+2}) - 4F_a(S_n)].
```

Diagnostic values at the current prefixes:

| Group | `n_Q` | `n_D` | `delta` | `S_Q` | `log10 F_Q` | Q model gap | `log10 F_D` | D model gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `E_6` | 77 | 75 | 2.000 | 11.85 | -21.28 | -13.68 | -17.80 | -13.68 |
| `E_6` | 77 | 75 | pi | 29.23 | -3.74 | 3.85 | -0.27 | 3.85 |
| `E_7` | 67 | 65 | 2.000 | 9.07 | -71.16 | -48.87 | -68.09 | -48.87 |
| `E_7` | 67 | 65 | pi | 22.37 | -29.35 | -7.05 | -26.27 | -7.05 |
| `E_8` | 67 | 65 | 2.000 | 8.10 | -200.39 | -132.45 | -198.29 | -132.45 |
| `E_8` | 67 | 65 | pi | 20.00 | -112.48 | -44.53 | -110.38 | -44.53 |

Inverting the same comparison gives the smallest cap size needed by this
truncated identity-saddle model:

| Group | `n_Q` | `n_D` | target `log10 F_Q` | required `delta_Q` | `delta_Q/pi` | target `log10 F_D` | required `delta_D` | `delta_D/pi` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `E_6` | 77 | 75 | -7.59 | 2.7532 | 0.876 | -4.12 | 2.7531 | 0.876 |
| `E_7` | 67 | 65 | -22.30 | 3.4452 | 1.097 | -19.22 | 3.4452 | 1.097 |
| `E_8` | 67 | 65 | -67.94 | 4.0665 | 1.294 | -65.84 | 4.0665 | 1.294 |

Thus the next finite-n theorem cannot be another principal-cap estimate with
`delta <= pi` for E_7/E_8.  The next named hard target is
`DCT-WideRegionTrig(G)`: prove a finite-n lower bound for `D_G` using
explicit trigonometric factors on an affine/wide identity region, or replace
that region by a different positive-region estimate whose constants are
valid from the current prefixes.

Rectangular radial-wedge diagnostic.  On the cap `delta=pi`, the checker
uses elementary inequalities for `2(1-cos t)`, a degree-six lower bound for
`2cos t`, and an exponential lower model for the Weyl sine factor.  It then
searches rectangles `s in [s0,s1]`, `t in [0,t1]` with
`t1 < (4/pi^2)s0`, plus the swapped rectangle.  The displayed margin is
`log10(lower D_{G,+} bound / negative-region bound)`:

The E6/E7 entries below used the now-rejected endpoints `C=2,3`; they are
retained only as historical search diagnostics and are not valid certificate
margins.  The active rank-corrected values are the C++ results stated above.

| Group | `n_D` | cap `S` | best rectangle margin | `s0` | `s1` | `t1` | exact `D_G` margin |
|---|---:|---:|---:|---:|---:|---:|---:|
| `E_6` | 75 | 28.47 | 37.60 | 21.71 | 22.60 | 7.70 | 67.23 |
| `E_7` | 65 | 21.71 | -7.32 | 20.08 | 21.57 | 7.12 | 44.65 |
| `E_8` | 65 | 19.40 | -102.05 | 17.95 | 19.40 | 6.36 | 16.24 |

This isolates the next split: `DCT-RectTrig(E_6)` can be handled by fixing a
clean rational rectangle and writing the elementary inequalities below;
E_7/E_8 require `DCT-WideRegionTrig` to use substantially more than a single
radial rectangle.

E6/E7 rank-corrected rectangle certificate.  The active C++ checker uses
the cap `delta=14/3` and proves the scalar cosine and Weyl-sine bounds by
exact Bernstein coefficients.  It generates the E6/E7 roots by simple
reflections and verifies coefficientwise

```text
16 sum_E6 (alpha.u)^4 = Q_E6(u)^2,
27 sum_E7 (alpha.u)^4 = Q_E7(u)^2.
```

For E6 it uses the valid crude negative bound `296*12^n`.  For E7 it proves
a degree-26 shifted-Chebyshev pointwise majorant on `[-14,-2] union [0,2]`
and integrates that majorant from `m_0,...,m_70`.  The replay returns:

```text
E6 retained_cells=23088  base_log_margin=8.904  step_log_margin=0.573
E7 retained_cells=66925  base_log_margin=0.426  odd_step_ratio>101/100
RESULT: ALL PASS
```

The exact prefixes overlap the tail starts at odd `n=75` and `n=65`.

E_8 rectangle diagnostic at the current prefix.  The same quartic-improved
rectangle-union method
uses

```text
sum_{alpha>0}(alpha.u)^4 = Q(u)^2/50
```

for E_8.  On an `80 x 80` grid, summing every admissible rectangle gives:

```text
E8 nD=65 cells=2583 margin=-68.51
```

So the tested E_7 rectangle-union mechanism does not close E_8.  The next
hard target is an E_8 positive-region estimate that captures more than the
local rectangular radial wedge.

E_8 higher-start tail certificate.  The same quartic identity does give a
higher-start rectangular-union certificate:

```text
Delta = 1/40,
s0 in {37,1481/40,...,793/20},
t0 in {26,1041/40,...,1239/40},
retained cells = 13558.
```

At `n=133`, the root cap is `1995 pi^2/496`, and
`pi > 333/106` gives the rational lower cap
`221223555/5573056 > 1587/40`.  The checker uses
`sin(z/2)/(z/2) >= exp(-z^2/21)` and `e < 11/4`.  It returns:

```text
retained_cells=13558
log10_union_fraction_lower=-66.75
h133_min=176378749/318402000
log10_ratio133_lower=0.02
log10_best_cell_ratio133=-3.25
log10_odd_step_ratio_lower=0.85
h_monotone_from_n133: OK
```

Thus `DCT-RectTrig(E_8-tail)` supplies a direct-tail dominance certificate
from odd `n >= 133`.

E_8 wide-cap bridge step.  The script
`direct_chain_rect_e8_n131_certificate.py` closes the subleaf
`E8-WideCapExp(131)`.  It uses the wider root-angle cap `delta=19/6`, the
degree-six E_8 identity `sum(alpha.u)^6 = Q(u)^3/1800`, and the
trigonometric bounds

```text
2(1-cos z) >= z^2 - z^4/12 + z^6/432,
2cos z >= 2 - z^2 + z^4/12 - z^6/360,
sin(z/2)/(z/2) >= exp(-z^2/21)
```

on `|z| <= 19/6`.  The exact-rational replay uses `e < 11/4` and returns:

```text
retained_cells=17797
log10_union_fraction_lower=-65.31
gap131_min=20929073/123559200000
h131_min=25824732993427/45523842750000
log10_ratio131_lower=1.59
log10_best_cell_ratio131=-1.41
ratio131_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n129_certificate.py`.  It uses the same wide-cap
trigonometric bounds, replacing the coarse exponential ceiling by the
fractional-part lower bound from `e < 1457/536` and
`exp(-r) >= 1-r+r^2/2-r^3/6` for `0 <= r < 1`.  The exact-rational replay
returns:

```text
retained_cells=12719
log10_union_fraction_lower=-66.02
gap129_min=26004883/647002080000
h129_min=4932872459749/8694090450000
log10_ratio129_lower=0.35
log10_best_cell_ratio129=-2.90
ratio129_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n127_certificate.py`.  It widens the cap to
`delta=10/3`, uses the valid scalar gap coefficient `1/450`, and keeps the
same fractional-part exponential lower bound.  The exact-rational replay
returns:

```text
retained_cells=10528
log10_union_fraction_lower=-60.43
gap127_min=604560121/3266122500000
h127_min=11143031103389/20739877875000
log10_ratio127_lower=1.04
log10_best_cell_ratio127=-2.32
ratio127_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n125_certificate.py`.  It widens the cap to
`delta=7/2`, keeps the valid scalar gap coefficient `1/450` and sine
denominator `21`, and trims the rectangular union to
`s0 >= 44`, `t0 in [31,67/2]`.  The exact-rational replay returns:

```text
retained_cells=8497
log10_union_fraction_lower=-55.42
gap125_min=25411339/117187500000
h125_min=2230567373321/4394531250000
log10_ratio125_lower=1.11
log10_best_cell_ratio125=-2.36
ratio125_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n123_certificate.py`.  It widens the cap to
`delta=4`, uses the valid scalar gap coefficient `1/480`, switches the sine
denominator to `20`, and trims the rectangular union to
`s0 in [53,59]`, `t0 in [69/2,75/2]`.  The exact-rational replay returns:

```text
retained_cells=28121
log10_union_fraction_lower=-42.55
gap123_min=170892133/408483000000
h123_min=647489491811/1507302270000
log10_ratio123_lower=1.10
log10_best_cell_ratio123=-3.19
ratio123_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n121_certificate.py`.  It keeps `delta=4`, uses the
valid scalar gap coefficient `1/475`, adds the Cauchy eighth-moment lower
bound, switches the sine denominator to `81/4`, and trims the rectangular
union to `s0 in [45,1815/31]`, `t0 in [30,42]`.  The exact-rational replay
returns:

```text
retained_cells=74053
log10_union_fraction_lower=-42.58
gap121_min=89528719/1173567656250
h121_min=239264176093148353904/572269131771240234375
log10_ratio121_lower=0.10
log10_best_cell_ratio121=-4.29
ratio121_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n119_certificate.py`.  It keeps `delta=4`, uses the
valid scalar gap coefficient `1/475`, adds the upper bound
`2(1-cos z) <= z^2 - z^4/12 + z^6/360` on the smaller-radius side of the
gap, and uses a coarser `1/10 x 1/10` rectangular union.  The exact-rational
replay returns:

```text
retained_cells=15000
log10_union_fraction_lower=-38.16
gap119_min=5527391/5935125
h119_min=7933413921950307734/19828183106689453125
log10_ratio119_lower=0.66
log10_best_cell_ratio119=-2.99
ratio119_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n117_certificate.py`.  It keeps `delta=4` and replaces
the Taylor coefficients in the gap and average-character estimates by
rational cubic scalar polynomials certified by exact Sturm root counts on
`0 <= z^2 <= 16`.  The exact-rational replay returns:

```text
retained_cells=111427
log10_union_fraction_lower=-33.40
gap117_min=183481292003/3422250000000000
h117_min=1317904655485969331/3580529062500000000
log10_ratio117_lower=0.05
log10_best_cell_ratio117=-4.38
gap_lower_sturm: OK
gap_upper_sturm: OK
cos_lower_sturm: OK
ratio117_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n115_certificate.py`.  It keeps the cubic scalar
polynomials and improves the Weyl-sine density loss to
`log(sin(z/2)/(z/2)) >= -z^2/24 - z^4/2100`, certified by an exact Sturm
root count on `0 <= (z/2)^2 <= 4`.  The exact-rational replay returns:

```text
retained_cells=98455
log10_union_fraction_lower=-35.13
gap115_min=850372679773/29756250000000000
h115_min=8242300592449041529/22100214843750000000
log10_ratio115_lower=0.52
log10_best_cell_ratio115=-3.85
sine_quartic_sturm: OK
ratio115_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n113_certificate.py`.  It keeps the same positive
rectangular machinery and improves the negative-region estimate using the
exact low moments

```text
m_0 = 1, m_1 = 0, m_2 = 1, m_3 = 1,
m_4 = 5, m_5 = 16, m_6 = 79.
```

The binomial expansion gives

```text
int int (X-Y)^2 (X+Y)^2 ((X+Y)^2+4) dX dY = 176,
```

so the negative-region bound for odd `n >= 3` is `176 * 16^(n-2)`.
The exact-rational replay returns:

```text
retained_cells=91581
negative_bound=176*16^111
log10_union_fraction_lower=-35.78
gap113_min=30644144709323/143651250000000000
h113_min=61950654956354875633/167736776250000000000
log10_ratio113_lower=1.89
log10_best_cell_ratio113=-2.44
low_moment_sources: OK
low_moment_negative_176: OK
ratio113_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n111_certificate.py`.  It reuses the same cubic scalar
polynomials, quartic Weyl-sine bound, and exact low-moment negative-region
bound at `N=111`.  The exact-rational replay returns:

```text
retained_cells=84702
negative_bound=176*16^109
log10_union_fraction_lower=-36.47
gap111_min=11131456399439/138611250000000000
h111_min=358832265265913163/981401875000000000
log10_ratio111_lower=0.39
log10_best_cell_ratio111=-3.88
low_moment_sources: OK
low_moment_negative_176: OK
ratio111_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n109_certificate.py`.  It uses the `p=12`
negative-region estimate

```text
J_12 = 5428655636,
J_14 = 1070740575644,
J_14 + 4J_12 = 1092455198188,
|D_-(n)| <= 1092455198188 * 16^(n-12).
```

The exact-rational replay returns:

```text
retained_cells=77454
negative_bound=1092455198188*16^97
log10_union_fraction_lower=-37.24
gap109_min=1262084876833/13366125000000000
h109_min=250022064633883337/690583125000000000
log10_ratio109_lower=1.15
log10_best_cell_ratio109=-3.06
moment_sources_0_16: OK
moment_negative_p12: OK
ratio109_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 109`.

The next bridge step is handled by
`direct_chain_rect_e8_n107_certificate.py`.  It uses the asymmetric
degree-28 negative-region estimate

```text
K_28 = sum_{i=0}^{14} binom(14,i) 496^{14-i}(-1)^i 512^{-14}
       * (J_{14+i}+4J_{12+i})
     = 9594616716049666093860048619540906615318349630451
       / 21267647932558653966460912964485513216,
|D_-(107)| <= (63/65) * 16^95 * K_28.
```

The exact-rational replay returns:

```text
retained_cells=70562
negative_bound=(63/65)*16^95*asym_integral
log10_asym_negative_vs_p12=-0.40
log10_union_fraction_lower=-38.04
gap107_min=1428012189147/14311250000000000
h107_min=3400562485123005593/9494083250000000000
log10_ratio107_lower=0.08
log10_best_cell_ratio107=-4.07
moment_sources_0_30: OK
asymmetric_negative_interval_sturm: OK
asymmetric_positive_interval_sturm: OK
asymmetric_moment_integral: OK
ratio107_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 107`.

The next bridge step is handled by
`direct_chain_rect_e8_n105_certificate.py`.  It uses the Chebyshev
degree-9 negative-region estimate

```text
C_9(S) = T_9(S/248 - 1) / T_9(-33/31),
K_28 = int int (X-Y)^2 S^8(S^2+4) C_9(S)^2 dX dY
     = 1374322971742349851054687958226748195749248467
       / 7483467728620557379409154785380638130176,
|D_-(105)| <= (63/65) * 16^97 * K_28.
```

The exact-rational replay returns:

```text
retained_cells=63660
negative_bound=(63/65)*16^97*cheb_integral
log10_cheb_negative_vs_p12=-1.97
log10_union_fraction_lower=-38.91
gap105_min=337595139377/24806250000000000
h105_min=5958587947540894453/16821738281250000000
log10_ratio105_lower=0.19
log10_best_cell_ratio105=-3.88
moment_sources_0_30: OK
chebyshev_negative_interval_sturm: OK
chebyshev_positive_interval_sturm: OK
chebyshev_moment_integral: OK
ratio105_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 105`.

The next bridge step is handled by
`direct_chain_rect_e8_n103_certificate.py`.  It uses the shifted Chebyshev
degree-14 negative-region estimate

```text
C_14(S) = T_14((S-249)/247) / T_14(-265/247),
K_38 = int int (X-Y)^2 S^8(S^2+4) C_14(S)^2 dX dY
     = 211367009391307597006963601980640081834787585716146581122825466841199954628
       / 101897318230030856294954778256196368953775863550819861054942528948894561,
|D_-(103)| <= (63/65) * 16^95 * K_38.
```

The exact-rational replay returns:

```text
retained_cells=56387
negative_bound=(63/65)*16^95*cheb_integral
log10_cheb_negative_vs_p12=-3.92
log10_union_fraction_lower=-39.97
gap103_min=2702191337/19096200000000
h103_min=3709027066736679977/10585792812500000000
log10_ratio103_lower=0.68
log10_best_cell_ratio103=-3.31
moment_log_0_40_oeis_overlap_0_30: OK
shifted_chebyshev_negative_interval_sturm: OK
shifted_chebyshev_positive_interval_sturm: OK
shifted_chebyshev_moment_integral: OK
ratio103_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 103`.

The next bridge step is handled by
`direct_chain_rect_e8_n101_certificate.py`.  It uses the shifted Chebyshev
degree-18 negative-region estimate

```text
C_18(S) = T_18((S-249)/247) / T_18(-265/247),
K_46 = int int (X-Y)^2 S^8(S^2+4) C_18(S)^2 dX dY
     = 128995043741757277312714846478298476789204294652455658518488602039975970471033893288196646788
       / 1728826333011022977416309538283938404536057433264315876432894452017451444007274700356544913,
|D_-(101)| <= (63/65) * 16^93 * K_46.
```

The exact-rational replay returns:

```text
retained_cells=49478
negative_bound=(63/65)*16^93*cheb_integral
log10_cheb_negative_vs_p12=-5.36
log10_union_fraction_lower=-41.19
gap101_min=38039712937/1136250000000000
h101_min=41469541186626839237/119772491250000000000
log10_ratio101_lower=0.68
log10_best_cell_ratio101=-3.24
moment_log_0_48_oeis_overlap_0_30: OK
shifted_chebyshev_negative_bernstein: OK
shifted_chebyshev_positive_bernstein: OK
shifted_chebyshev_moment_integral: OK
ratio101_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 101`.

The next two bridge steps are handled by
`direct_chain_rect_e8_n99_n97_certificate.py`.  They use the shared shifted
Chebyshev degree-26 negative-region estimate

```text
C_26(S) = T_26((S-249)/247) / T_26(-265/247),
K_62 = int int (X-Y)^2 S^8(S^2+4) C_26(S)^2 dX dY
     = 652708718049108501279891206610792628257515262318580914182134886342747101135881782357617401218658855060065796860047212614264650354628
       / 2445200692443756455056633566892545977519407487004403208271763856863104974462205280301950124320130386149081237313785855952197768855841,
|D_-(n)| <= (63/65) * 16^(n-8) * K_62,  n in {99,97}.
```

The exact-rational replay returns:

```text
case n=99 s_end=479/10 t_end=479/10
retained_cells=42662 cap=1485/31
negative_bound=(63/65)*16^91*cheb_integral
log10_cheb_negative_vs_p12=-7.81
log10_union_fraction_lower=-42.52
gap99_min=1856998249/13750000000000
h99_min=3540441037296109973/10254296250000000000
log10_ratio99_lower=1.70
log10_best_cell_ratio99=-2.15
case n=97 s_end=469/10 t_end=469/10
retained_cells=35869 cap=1455/31
negative_bound=(63/65)*16^89*cheb_integral
log10_cheb_negative_vs_p12=-7.81
log10_union_fraction_lower=-43.97
gap97_min=11851190925583/105851250000000000
h97_min=47220675071042327/136724531250000000
log10_ratio97_lower=0.27
log10_best_cell_ratio97=-3.50
moment_log_0_64_oeis_overlap_0_30: OK
shifted_chebyshev_negative_bernstein_99: OK
shifted_chebyshev_positive_bernstein_99: OK
ratio99_gt_1: OK
shifted_chebyshev_negative_bernstein_97: OK
shifted_chebyshev_positive_bernstein_97: OK
ratio97_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 97`.

The next bridge step is handled by
`direct_chain_rect_e8_n95_certificate.py`.  It uses the shifted Chebyshev
degree-27 negative-region estimate on `[7,496]`

```text
C_27(S) = T_27((2S-503)/489) / T_27(-535/489),
K_70 = int int (X-Y)^2 S^12(S^2+4) C_27(S)^2 dX dY
     = 13210904703183651885179876881244444531729278568775924088504240739843154251964895373471546007671055851963458908535843318949535431873701818608573195108157158828
       / 51853866904436819647490281676508755811733798973504689872181530854918857454734960543821244670502133670706573447042449367379208198446013146633495281225231025,
|D_-(95)| <= (63/65) * 16^83 * K_70.
```

The exact-rational replay returns:

```text
retained_cells=29784 cap=1425/31
negative_bound=(63/65)*16^83*cheb_integral
log10_cheb_negative_vs_p12=-9.65
log10_union_fraction_lower=-45.39
gap95_min=9482346092687/101531250000000000
h95_min=17205807832136089069/49834921875000000000
log10_ratio95_lower=0.70
log10_best_cell_ratio95=-3.01
moment_log_0_70_oeis_overlap_0_30: OK
shifted_chebyshev_negative_bernstein: OK
shifted_chebyshev_positive_bernstein: OK
shifted_chebyshev_moment_integral: OK
ratio95_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 95`.

The next bridge step is handled by
`direct_chain_rect_e8_n93_certificate.py`.  It uses the widened root-angle
cap `delta=21/5`, the adjusted scalar lower coefficient `7/4000`, the
quartic sine bound `log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2000`, and the same
shifted Chebyshev degree-27 negative-region estimate

```text
C_27(S) = T_27((2S-503)/489) / T_27(-535/489),
K_70 = int int (X-Y)^2 S^12(S^2+4) C_27(S)^2 dX dY
     = 13210904703183651885179876881244444531729278568775924088504240739843154251964895373471546007671055851963458908535843318949535431873701818608573195108157158828
       / 51853866904436819647490281676508755811733798973504689872181530854918857454734960543821244670502133670706573447042449367379208198446013146633495281225231025,
|D_-(93)| <= (63/65) * 16^81 * K_70.
```

The exact-rational replay returns:

```text
retained_cells=53822 cap=3969/80
wide_gap_lower_coeffs=99119225/100000000,-7802858/100000000,7/4000
sine_bound=exp(-z^2/24-z^4/2000)
negative_bound=(63/65)*16^81*cheb_integral
log10_cheb_negative_vs_p12=-9.65
log10_union_fraction_lower=-40.58
gap93_min=7225097/40500000000
h93_min=962182419750323/3138750000000000
log10_ratio93_lower=0.04
log10_best_cell_ratio93=-4.04
wide_gap_lower_sturm: OK
wide_gap_upper_sturm: OK
wide_cos_lower_sturm: OK
sine_quartic_21_5_sturm: OK
moment_log_0_70_oeis_overlap_0_30: OK
shifted_chebyshev_negative_bernstein: OK
shifted_chebyshev_positive_bernstein: OK
shifted_chebyshev_moment_integral: OK
ratio93_gt_1: OK
```

The next bridge step is handled by
`direct_chain_rect_e8_n91_certificate.py`.  It uses root-angle cap
`delta=13/3`, the adjusted scalar lower coefficient `1562/890661`, the
sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/115000`, and an optimized
degree-27 square negative-region majorant

```text
|D_-(91)| <= 16^79 * int int (X-Y)^2 S^12(S^2+4)P_27(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=69502 cap=76895/1488
gap_lower_coeffs=99119225/100000000,-7802858/100000000,1562/890661
sine_bound=exp(-z^2/24-z^4/2880-z^6/115000)
negative_bound=16^79*optimized_degree27_integral
log10_union_fraction_lower=-38.10
gap91_min=29512424195754583/2676615873750000000000
h91_min=1173568851885107429/4171553750000000000
log10_positive91_lower=96.86
log10_negative91_upper=96.73
log10_ratio91_lower=0.13
n91_gap_lower_sturm: OK
n91_gap_upper_sturm: OK
n91_cos_lower_sturm: OK
sine_sextic_13_3_sturm: OK
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
optimized_negative_moment_integral: OK
ratio91_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 91`.

The next bridge step is handled by
`direct_chain_rect_e8_n89_certificate.py`.  It uses root-angle cap
`delta=14/3`, the adjusted scalar lower coefficient `169078/100000000`, the
sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and an optimized
degree-28 square negative-region majorant

```text
|D_-(89)| <= 16^77 * int int (X-Y)^2 S^12(S^2+4)P_28(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=130740 cap=21805/372
gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000
sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)
negative_bound=16^77*optimized_degree28_integral
log10_union_fraction_lower=-31.17
gap89_min=4441364733419/44555625000000000
h89_min=9192686056706299193/40976323125000000000
log10_positive89_lower=93.44
log10_negative89_upper=92.97
log10_ratio89_lower=0.47
n89_gap_lower_sturm: OK
n89_gap_upper_sturm: OK
n89_cos_lower_sturm: OK
sine_sextic_14_3_sturm: OK
moment_log_0_72_arxiv_extension: OK
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
optimized_negative_moment_integral: OK
ratio89_gt_1: OK
```

#### E8 bridge n=87

`direct_chain_rect_e8_n87_certificate.py` closes

```text
D_E8(87) = Q_3^{E8}(89)-4Q_3^{E8}(87).
```

It keeps the `n=89` cap `14/3`, grid step `1/20`, cap-adjusted scalar gap,
and sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and raises the
negative-region square majorant to degree `30`:

```text
|D_-(87)| <= 16^75 * int int (X-Y)^2 S^12(S^2+4)P_30(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=116526 cap=7105/124
gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000
sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)
negative_bound=16^75*optimized_degree30_integral
log10_union_fraction_lower=-32.59
gap87_min=1828049159293/85151250000000000
h87_min=954625549096143173/4252831875000000000
log10_positive87_lower=89.87
log10_negative87_upper=89.43
log10_ratio87_lower=0.44
n87_gap_lower_sturm: OK
n87_gap_upper_sturm: OK
n87_cos_lower_sturm: OK
sine_sextic_14_3_sturm: OK
moment_log_0_76_arxiv_extension: OK
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
optimized_negative_moment_integral: OK
ratio87_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 87`.

#### E8 bridge n=85

`direct_chain_rect_e8_n85_certificate.py` closes

```text
D_E8(85) = Q_3^{E8}(87)-4Q_3^{E8}(85).
```

It keeps the `n=89` cap `14/3`, grid step `1/20`, cap-adjusted scalar gap,
and sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and raises the
negative-region square majorant to degree `34`:

```text
|D_-(85)| <= 16^73 * int int (X-Y)^2 S^12(S^2+4)P_34(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=103495 cap=20825/372
gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000
sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)
negative_bound=16^73*optimized_degree34_integral
log10_union_fraction_lower=-34.01
gap85_min=3855342619/37353515625000
h85_min=213613375085470649/951893750000000000
log10_positive85_lower=86.34
log10_negative85_upper=86.05
log10_ratio85_lower=0.29
n85_gap_lower_sturm: OK
n85_gap_upper_sturm: OK
n85_cos_lower_sturm: OK
sine_sextic_14_3_sturm: OK
moment_log_0_84_arxiv_extension: OK
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein_32: OK
optimized_negative_moment_integral: OK
ratio85_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 85`.

#### E8 bridge n=83

`direct_chain_rect_e8_n83_certificate.py` closes

```text
D_E8(83) = Q_3^{E8}(85)-4Q_3^{E8}(83).
```

It keeps the `n=89` cap `14/3`, grid step `1/20`, cap-adjusted scalar gap,
and sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and raises the
negative-region square majorant to degree `38`:

```text
|D_-(83)| <= 16^71 * int int (X-Y)^2 S^12(S^2+4)P_38(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=91138 cap=20335/372
gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000
sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)
negative_bound=16^71*optimized_degree38_integral
log10_union_fraction_lower=-35.50
gap83_min=1339382274583/19375312500000000
h83_min=119188825916243597/531761910000000000
log10_positive83_lower=82.83
log10_negative83_upper=82.55
log10_ratio83_lower=0.28
n83_gap_lower_sturm: OK
n83_gap_upper_sturm: OK
n83_cos_lower_sturm: OK
sine_sextic_14_3_sturm: OK
moment_log_0_92_arxiv_extension: OK
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
optimized_negative_moment_integral: OK
ratio83_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 83`.

#### E8 bridge n=81

`direct_chain_rect_e8_n81_certificate.py` closes

```text
D_E8(81) = Q_3^{E8}(83)-4Q_3^{E8}(81).
```

It keeps the `n=89` cap `14/3`, grid step `1/20`, cap-adjusted scalar gap,
and sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and raises the
negative-region square majorant to degree `42`:

```text
|D_-(81)| <= 16^69 * int int (X-Y)^2 S^12(S^2+4)P_42(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=79060 cap=6615/124
gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000
sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)
negative_bound=16^69*optimized_degree42_integral
log10_union_fraction_lower=-37.10
gap81_min=2237512207/164025000000000
h81_min=1733570508555555541/7722502031250000000
log10_positive81_lower=79.34
log10_negative81_upper=79.08
log10_ratio81_lower=0.27
n81_gap_lower_sturm: OK
n81_gap_upper_sturm: OK
n81_cos_lower_sturm: OK
sine_sextic_14_3_sturm: OK
moment_log_0_100_arxiv_extension: OK
optimized_negative_interval_bernstein_8: OK
optimized_positive_interval_bernstein: OK
optimized_negative_moment_integral: OK
ratio81_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 81`.

#### E8 bridge n=79

`direct_chain_rect_e8_n79_certificate.py` closes

```text
D_E8(79) = Q_3^{E8}(81)-4Q_3^{E8}(79).
```

It keeps the `n=89` cap `14/3`, grid step `1/20`, cap-adjusted scalar gap,
and sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and uses a
degree-37 square negative-region majorant at moment weight `p=22`:

```text
|D_-(79)| <= 16^57 * int int (X-Y)^2 S^22(S^2+4)P_37(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=68095 cap=19355/372
gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000
sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)
negative_bound=16^57*optimized_degree37_integral
log10_union_fraction_lower=-38.71
gap79_min=2975877287857/35105625000000000
h79_min=2141716707210133747/9552630625000000000
log10_positive79_lower=75.89
log10_negative79_upper=75.10
log10_ratio79_lower=0.79
n79_gap_lower_sturm: OK
n79_gap_upper_sturm: OK
n79_cos_lower_sturm: OK
sine_sextic_14_3_sturm: OK
moment_log_0_100_arxiv_extension: OK
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
optimized_negative_moment_integral: OK
ratio79_gt_1: OK
```

Therefore the remaining E_8 direct-tail work is reduced to the finite
bridge `67 < n <= 79`.

#### E8 bridge n=77

`direct_chain_rect_e8_n77_certificate.py` closes

```text
D_E8(77) = Q_3^{E8}(79)-4Q_3^{E8}(77).
```

It keeps the `n=89` cap `14/3`, grid step `1/20`, cap-adjusted scalar gap,
and sextic sine bound
`log(sin(z/2)/(z/2)) >= -z^2/24-z^4/2880-z^6/103800`, and uses a
degree-37 square negative-region majorant at moment weight `p=22`:

```text
|D_-(77)| <= 16^55 * int int (X-Y)^2 S^22(S^2+4)P_37(S)^2 dX dY.
```

The exact-rational replay returns:

```text
retained_cells=57817 cap=18865/372
gap_lower_coeffs=99119225/100000000,-7802858/100000000,169078/100000000
sine_bound=exp(-z^2/24-z^4/2880-z^6/103800)
negative_bound=16^55*optimized_degree37_integral
log10_union_fraction_lower=-40.40
gap77_min=1815834031/275625000000000
h77_min=3964993922572267201/17690653750000000000
log10_positive77_lower=72.46
log10_negative77_upper=72.26
log10_ratio77_lower=0.19
n77_gap_lower_sturm: OK
n77_gap_upper_sturm: OK
n77_cos_lower_sturm: OK
sine_sextic_14_3_sturm: OK
moment_log_0_100_arxiv_extension: OK
optimized_negative_interval_bernstein: OK
optimized_positive_interval_bernstein: OK
optimized_negative_moment_integral: OK
ratio77_gt_1: OK
```

Therefore the rectangular direct-tail work is reduced to the finite bridge
`67 < n <= 77`.

#### E8 arXiv m100 finite bridge
`direct_chain_e8_arxiv_m100_finite_bridge.py` evaluates the exact finite
formula for `Q_3^{E8}(n)` from the source-checked `m_0,...,m_100` moment
table and verifies the remaining bridge directly.

The exact integer replay returns:

```text
group=E8 source=m_0..m_100
bridge_odd_n=(69, 71, 73, 75, 77)
chain_odd_n=(67, 69, 71, 73, 75, 77)
extended_chain_odd_n=(67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95)
log10_min_bridge_q3=100.94
log10_min_chain_diff=100.94
log10_min_extended_chain_diff=100.94
moment_log_0_100_arxiv_extension: OK
bridge_q3_69_77_positive: OK
chain_diff_67_77_positive: OK
extended_chain_diff_67_95_positive: OK
```

Thus `E8-FiniteBridge(77)` is discharged by the exact arXiv m100 finite
bridge replay.

## Route to full closure

### Path A: source-audited theorem route
1. Use `DCT_RECT_TRIG_G2_F4.md` and its C++/GMP replay for the G2/F4 tails.
2. Use the source-checked arXiv m100 finite-bridge replay for
   `E8-FiniteBridge(77)`, and use the E_8 tail-and-bridge chain from odd
   `n >= 77`.  The E_6 and E_7 slots are handled by
   `DCT_RECT_TRIG_E6.md` and `DCT_RECT_TRIG_E7.md`; the E_8 tail and the
   `n=131,129,127,125,123,121,119,117,115,113,111,109,107,105,103,101,99,97,95,93,91,89,87,85,83,81,79,77` bridge steps are handled by
   `DCT_RECT_TRIG_E8_TAIL.md`, `DCT_RECT_TRIG_E8_N131.md`,
   `DCT_RECT_TRIG_E8_N129.md`, `DCT_RECT_TRIG_E8_N127.md`, and
   `DCT_RECT_TRIG_E8_N125.md`, `DCT_RECT_TRIG_E8_N123.md`, and
   `DCT_RECT_TRIG_E8_N121.md`, `DCT_RECT_TRIG_E8_N119.md`,
   `DCT_RECT_TRIG_E8_N117.md`, `DCT_RECT_TRIG_E8_N115.md`, and
   `DCT_RECT_TRIG_E8_N113.md`, `DCT_RECT_TRIG_E8_N111.md`,
   `DCT_RECT_TRIG_E8_N109.md`, `DCT_RECT_TRIG_E8_N107.md`, and
   `DCT_RECT_TRIG_E8_N105.md`, `DCT_RECT_TRIG_E8_N103.md`,
   `DCT_RECT_TRIG_E8_N101.md`, `DCT_RECT_TRIG_E8_N99_N97.md`, and
   `DCT_RECT_TRIG_E8_N95.md`, `DCT_RECT_TRIG_E8_N93.md`,
   `DCT_RECT_TRIG_E8_N91.md`, `DCT_RECT_TRIG_E8_N89.md`, and
   `DCT_RECT_TRIG_E8_N87.md`, `DCT_RECT_TRIG_E8_N85.md`, and
   `DCT_RECT_TRIG_E8_N83.md`, `DCT_RECT_TRIG_E8_N81.md`,
   `DCT_RECT_TRIG_E8_N79.md`, `DCT_RECT_TRIG_E8_N77.md`, and
   `E8_ARXIV_M100_FINITE_BRIDGE.md`.
3. Use the resulting tail bridge with the exact Chain prefix or extended
   exact Chain window.

### Path B: compute-heavy character-ring
1. Extend exact moments only after a sharper tail bridge gives a finite
   target close to the current range; the current exact-c₃ thresholds alone
   ask for Chain `m >= 225, 636, 2168` for E_6, E_7, E_8 respectively.
2. Direct integer-arithmetic verification of shift-2 log-convexity to the
   sharpened target.
3. Combine that finite bridge with the exact tail theorem.

### Path C: Hybrid
1. Extend character-ring where cheap (E_7 needs ~3 more iterations past current m=70).
2. Rigorous Wong `K_univ` for F_4, E_6 where `c_3` is tractable.
3. GPU hafnian for the rest.
