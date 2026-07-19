# Classical m3=1 source map

This note isolates the current source map for the non-SU classical families
in the adjoint Q3 problem:

```text
A1 = SU(2) = SO(3),  B_n = SO(2n+1),  C_n = Sp(n),  D_n = SO(2n).
```

It is an integration note.  The final proof text should replace each
obligation below by a numbered theorem, proposition, or cited computation.

## Proposition 1 (classification slot)

The active route splits the classical compact simple adjoint families as
follows:

```text
A_n for n >= 2: handled by SU_N_REPAIR.md.
A1: handled by CLASSICAL_M3_ONE_PROOF.md, Propositions 1-2.
B_n, C_n, D_n: remaining classical m3=1 tail slot.
```

Source location: `NONABELIAN.md`, section "Consolidated closure across all
compact simple gauge groups", especially the invariant-degree classification.

Status: the classical symmetric-cubic degree table is closed in
`CLASSICAL_M3_ONE_PROOF.md`, Proposition 3.  The exceptional families are
handled separately by the direct-chain certificates recorded in
`character_ring_iter/CLOSURE_STATUS.md`.

## Proposition 2 (A1/SU2 slot)

The `A1 = SU(2) = SO(3) = Sp(1)` adjoint case is covered by
`CLASSICAL_M3_ONE_PROOF.md`, Propositions 1-2, using `SU_N_REPAIR.md`,
Theorem 71, at `N = 2`.

Status: this leaf is closed in the standalone extraction note.

## Proposition 3 (classical stable-family slot)

The historical diary records a Wong/Chain route for the m3=1 classical
families.  `CLASSICAL_M3_ONE_PROOF.md`, Propositions 5-7, now extract the
SO/Sp adjoint trace formulas and the stable `Sp(infinity)=SO(infinity)`
moment equality.  Propositions 8-9 extract the exact Q3 moment formula and
the `m_0,...,m_{2r+5}` range needed for Chain step `r`.  Propositions 10-13
close the stable SO/Sp Chain leaf, and Propositions 15-16 close the
stable-rank moment-stabilization reduction in explicit rank ranges.  Lemma
16A.1 and Lemma 16A isolate the determinant-degree `SO(2b)` boundary correction
used for the first `D` middle-slice row, and Lemma 16B extends that correction
through depth four.  Lemma 16C gives the finite Pieri recurrence used for
deeper determinant-boundary rows.
Proposition 17 gives the Weyl torus formula for Q3, and Proposition 18
records the Wong local Laplace theorem with the vanishing-amplitude coefficient
extraction.  Propositions 19-20 compute the identity-saddle expansion and the
positive leading Gaussian coefficient, Proposition 21 bounds the negative
region, Proposition 22 isolates the central saddle set and bounds the
positive secondary region away from it, and Proposition 23 makes the
negative-region base strictly subleading.  Proposition 24 supplies the
explicit positive-saddle tail threshold; Lemmas 24.1-24.2 and Corollary 24.3
insert the classical root-data constants and the exact negative-base constants
`C_G` for `B_b`, `C_b`, and `D_b`; Lemma 24.4 shows that this sufficient
cutoff does not meet the stable Chain window in high rank.  Proposition 24.5
computes the classical Macdonald-Mehta leading coefficient, and Lemma 24.6
shows that the leading coefficient is asymptotically large enough throughout
the post-stable range if a rank-growing tail theorem can make it available.
Lemma 24.7 rules out the fixed local-cap half-leading route at the stable
edge.  Proposition 24.8 and Corollary 24.9 reduce the resulting nonlocal
tail option to a two-sided comparison between the one-dimensional pushforward
tails `P_G(c)` and `N_G(a)`.  Proposition 24.10 reduces that comparison to
defining-trace deviations, and Corollary 24.11A fixes the quadratic-tau
high-rank target `a_G=eta C_G`, `c_G=r C_G`, with explicit rate inequalities
for `(r,eta,A,C_0^quad)`.  Proposition 24.12 proves the `tau_G` trace tail
from the Courteaut-Johansson-Lambert/Rains `m=2` concentration estimate while
retaining the Rains correction, Corollary 24.13 pins the high-rank
positive-trace target to
`Pr(|T_1|>=sqrt((22991901/5000000)C_G))>=exp(-(461/200)C_G)`,
Proposition 24.14
proves that positive-trace target by Courteaut-Johansson-Lambert one-trace
total variation, and Corollary 24.15 closes the high-rank SO/Sp tail.
Proposition 24.16B then replaces the first `945` cutoff by the Rains-square
cutoff `296`, using Courteaut-Johansson-Lambert Theorems 1.1, 1.4, and 4.1
and the sharpened negative-tail mass bound in Proposition 24.10.
Proposition 25 supplies the exact classical Chain bridge through odd `n=23`,
Proposition 26 extends the exact B/C bridge through odd `n=25`, and
Proposition 27 extends the exact D bridge through odd `n=25` for ranks
`4 <= b <= 27`.  Proposition 28 extends the exact B/C bridge through odd
`n=27`, Proposition 29 extends the exact D bridge through odd `n=27`, and
Proposition 30 extends the exact classical bridge through odd `n=29`.
Proposition 31 extends the exact classical bridge through odd `n=31`.
The diagnostic replay
`character_ring_iter/classical_tail_constants.py --trace-cutoff-certificate`
checks the displayed monotone trace-bridge inequalities with
`r=20000001/10000000`, `eta=9973/50000`, `A=461/200`, and the retained tau
correction from `C_G>=945`.  It gives the concrete high-rank cutoff
`C_classical_high=945`.  The sharper diagnostic replay
`character_ring_iter/classical_tail_constants.py --trace-square-cutoff-certificate`
checks the Rains-square trace bridge with `r=2000001/1000000`,
`eta=359/2000`, `A=4571/2000`, and gives the concrete high-rank cutoff
`C_classical_high=296`.  Lemma 24.17B identifies the resulting finite residue
as `881` rows:

```text
B_b: 2 <= b <= 295,
C_b: 2 <= b <= 295,
D_b: 4 <= b <= 295 and b=297.
```

It also defines the row-gated middle bridge target

```text
S_mid^296(G) =
  { m >= 0 : stable_odd_reach(G) < 2m+3 < N_loc(G) },

M_mid^296 =
  max { m : m in S_mid^296(G) for some finite row G with C_G<296 }.
```

The diagnostic replay
`character_ring_iter/classical_tail_constants.py --finite-leading-certificate`
checks that the leading-dominance inequality is already positive from the
stable edge throughout the finite residue; the worst certified row is `B_2`
at odd `n=7`, with log10 margin `0.429493534104`.

Current SO/Sp hard-target meter:

```text
CMO-ClassicalTail active meter: conditional, not 100/100

closed  high-rank trace tail, C_G >= 296, via Proposition 24.16B
closed  row-gated exact bridge through m=29, via Propositions 31A--31FRCD
partial finite-row post-bridge tail remainder:
       57 top finite rows closed at first hit by Proposition 24.17D,
       diagnostic slices m=26..29 closed by Propositions 31FRD--31FRCD,
       the older
       CMO-LocalOnsetBound29 route is obstructed by Proposition 31FRCG;
       post-m29 direct tail rows closed so far: B_2/C_2 by Proposition
       31FRCI, B_3/C_3 by Proposition 31FRCJ, D_4 by Proposition 31FRCK,
       B_4/C_4 by Proposition 31FRCL, D_5 by Proposition 31FRCM, and
       C_5 by Proposition 31FRCN, B_5 by Proposition 31FRCO, C_6 by
       Proposition 31FRCP, D_6/B_6 by Proposition 31FRCR, and
       D_7/B_7/C_7 by Proposition 31FRCS, C_8 by Proposition 31FRCT, and
       D_8 by Proposition 31FRCV, B_8/C_9/D_9 by Proposition 31FRCX,
       C_10 by Proposition 31FRCY, B_9 by Proposition 31FRDB, C_11
       by Proposition 31FRDE, C_12 by Proposition 31FRDI, B_10 by
       Proposition 31FRDN, B_11 by Proposition 31FRDN3, B_12 by
       Proposition 31FRDN6, C_13 by Proposition 31FRDN9, C_14 by
       Proposition 31FRDN12, D_10 by Proposition 31FRDR, D_11 by
       Proposition 31FRDV, C_15 by Proposition 31FRDW, C_16 by
       Proposition 31FRDW2, C_17 by Proposition 31FRDW3, C_18 by
       Proposition 31FRDW4, C_19 by Proposition 31FRDW5, B_13 by
       Proposition 31FRDW6, and D_12..D_97 by Propositions
       31FRDW7--31FRDW92; Proposition
       24.17E closes the 118 high-edge rows B_b/C_b with 218<=b<=275
       or b=277 using the source-supported no-TV Chernoff
       trace supplier; paper Proposition `prop:post29-d-interval` closes
       D_98..D_275,D_277,D_279, and paper Proposition
       `prop:post29-b-unitary-square` closes B_124..B_217 by the
       source-audited Rains/CJL unitary square-trace supplier, and paper
       Proposition `prop:post29-b-tilted-mgf` plus the proved correction
       prefix closes B_62..B_123, while paper Proposition
       `prop:post29-c-tilted-mgf` plus the same prefix closes C_62..C_217;
       Proposition
       24.17F gives only the
       first-hit Chain step for B ranks 11..275 and 277, C ranks 13..275 and
       277, and D ranks 12..275,277,279; after Propositions 31FRDN3,
       31FRDN6, 31FRDN9, 31FRDN12, 31FRDW, 31FRDW2, 31FRDW3,
       31FRDW4, 31FRDW5, 31FRDW6, 31FRDW7, 31FRDW8, 31FRDW9,
       31FRDW10, 31FRDW11, 31FRDW12, 31FRDW13, 31FRDW14,
       31FRDW15, 31FRDW16, 31FRDW17, 31FRDW18, 31FRDW19, 31FRDW20,
       31FRDW21, 31FRDW22, 31FRDW23, 31FRDW24, 31FRDW25, 31FRDW26, 31FRDW27, 31FRDW28, 31FRDW29, 31FRDW30, 31FRDW31, 31FRDW32, 31FRDW33, 31FRDW34, 31FRDW35, 31FRDW36, 31FRDW37, 31FRDW38, 31FRDW39, 31FRDW40, 31FRDW41, 31FRDW42, 31FRDW43, 31FRDW44, 31FRDW45, 31FRDW46, 31FRDW47, 31FRDW48, 31FRDW49, 31FRDW50, 31FRDW51, 31FRDW52, 31FRDW53, 31FRDW54, 31FRDW55, 31FRDW56, 31FRDW57, 31FRDW58, 31FRDW59, 31FRDW60, 31FRDW61, 31FRDW62, 31FRDW63, 31FRDW64, 31FRDW65, 31FRDW66, 31FRDW67, 31FRDW68, 31FRDW69, 31FRDW70, 31FRDW71, 31FRDW72, 31FRDW73, 31FRDW74, 31FRDW75, 31FRDW76, 31FRDW77, 31FRDW78, 31FRDW79, 31FRDW80, 31FRDW81, 31FRDW82, 31FRDW83, 31FRDW84, 31FRDW85, 31FRDW86, 31FRDW87, 31FRDW88, 31FRDW89, 31FRDW90, 31FRDW91, 31FRDW92, and 24.17E remove B_11,
       B_12, B_13, C_13, C_14, C_15, C_16, C_17, C_18, C_19,
       D_12, D_13, D_14, D_15, D_16, D_17, D_18, D_19, D_20, D_21, D_22, D_23, D_24, D_25, D_26, D_27, D_28, D_29, D_30, D_31, D_32, D_33, D_34, D_35, D_36, D_37, D_38, D_39, D_40, D_41, D_42, D_43, D_44, D_45, D_46, D_47, D_48, D_49, D_50, D_51, D_52, D_53, D_54, D_55, D_56, D_57, D_58, D_59, D_60, D_61, D_62, D_63, D_64, D_65, D_66, D_67, D_68, D_69, D_70, D_71, D_72, D_73, D_74, D_75, D_76, D_77, D_78, D_79, D_80, D_81, D_82, D_83, D_84, D_85, D_86, D_87, D_88, D_89, D_90, D_91, D_92, D_93, D_94, D_95, D_96, D_97, and the 118
       high-edge rows; after the two paper propositions also remove the D
       interval, B_62..B_217, and C_62..C_217, and paper Proposition
       `prop:post29-bc-layered-mgf` removes B_14..B_61,C_20..C_61, including
       the exact Toeplitz--Hankel determinant closure of B_14,B_15; no row
       remains in CMO-PostM29AllLaterTail
       diagnostic m=21 bridge submeter: 86/86 rows closed
       (46 status-backed rows by Proposition 31DR, D_27..D_33 by
       Proposition 31DS, D_20..D_26 by Proposition 31DT, D_12..D_18 by
       Proposition 31DU, D_9..D_11 by Proposition 31DV, and B/C ranks
       9..16 by Proposition 31DW; this submeter is not the global
       denominator)
       next diagnostic m=22 bridge submeter: 90/90 rows closed
       (B/C ranks 9..23 by Proposition 31DX and D_29..D_49 by
       Proposition 31DY, D_18..D_28 by Proposition 31DZ, and D_12..D_17
       by Proposition 31EA, B/C ranks 6..8 by Proposition 31EB, and B_5
       by Proposition 31EC, D_4..D_11 by Proposition 31ED, and C_5 by
       Proposition 31EE, B_4 by Proposition 31EF, and C_4 by Proposition
       31EG, B_3/C_3 by Proposition 31EH, and B_2/C_2 by Proposition
       31EI; this submeter is not the global denominator)
       next diagnostic m=23 bridge submeter: 94/94 rows closed
       (B_2/C_2 by Proposition 31EJ, B/C ranks 16..24 by Proposition
       31EK, D_4..D_51 by Propositions 31EL--31EO, B/C ranks 9..15 by
       Proposition 31EP, B/C rank 8 by Proposition 31EQ, and B/C rank 7
       by Proposition 31ER, B/C rank 6 by Proposition 31ES, and B/C rank
       5 by Proposition 31ET, B/C rank 4 by Proposition 31EU, and B/C
       rank 3 by Proposition 31EV; this submeter is not the global
       denominator)
       next diagnostic m=24 bridge submeter: 98/98 rows closed
       (B_2/C_2 by Proposition 31EW, B/C ranks 16..25 by Proposition
       31EX, D_32..D_53 by Proposition 31EY, D_21..D_31 by Proposition
       31EZ, D_13..D_20 by Proposition 31FA, D_4 by Proposition 31FB,
       B/C rank 15 by Proposition 31FC, D_5..D_8 by Propositions
       31FD--31FG, D_9..D_12 by Proposition 31FH, B/C ranks 4..14 by
       Proposition 31FI, and B/C rank 3 by Proposition 31FJ; this
       submeter is not the global denominator)
       next diagnostic m=25 bridge submeter: 102/102 rows closed
       (B_2/C_2 by Proposition 31FL, B_3/C_3 by Proposition 31FM,
       B_4/C_4 by Proposition 31FN, B/C ranks 5..16 by Proposition
       31FO, B/C ranks 17..26 by Proposition 31FP, D_4..D_13 by
       Proposition 31FQ, and D_14..D_55 by Proposition 31FR; this
       submeter is not the global denominator)
       next diagnostic m=26 bridge submeter: 106/106 rows closed
       (D_5..D_11 by Propositions 31FRD--31FRI and 31FRO, B/C ranks
       17..27 by Proposition 31FRJ, B_2/C_2 by Proposition 31FRK,
       B_3/C_3 by Proposition 31FRL, B/C ranks 5..16 by Proposition
       31FRM, D_14..D_57 by Proposition 31FRN, D_4 by Proposition
       31FRP, D_12/D_13 by Proposition 31FRQ, B_4 by Proposition
       31FRR, and C_4 by Proposition 31FRS; this submeter is not the
       global denominator)
       next diagnostic m=27 bridge submeter: 110/110 rows closed
       (B/C ranks 18..28 by Proposition 31FRT, D_14..D_59 by
       Proposition 31FRU, B_2/C_2, B_3/C_3, and B/C ranks 6..16 by
       Proposition 31FRV, B_5/C_5 by Proposition 31FRW, D_13 by
       Proposition 31FRX, D_12 by Proposition 31FRY, D_11..D_5 by
       Propositions 31FRZ--31FRAF, D_4 by Proposition 31FRAG,
       B_17/C_17 by Proposition 31FRAH, and B_4/C_4 by Proposition
       31FRAI; this submeter is not the global denominator)
       next diagnostic m=28 bridge submeter: 114/114 rows closed
       (B/C ranks 19..29 by Proposition 31FRAJ and D_15..D_61 by
       Proposition 31FRAK, B_2/C_2 by Proposition 31FRAL, B_3/C_3 by
       Proposition 31FRAM, B_18/C_18 by Proposition 31FRAN, B/C ranks
       7..14 by Proposition 31FRAO, B_16/C_16 by Proposition 31FRAP,
       B_6/C_6 by Proposition 31FRAQ, B_17/C_17 by Proposition
       31FRAR, B_15/C_15 by Proposition 31FRAS, D_14 by Proposition
       31FRAT, D_13 by Proposition 31FRAU, D_12 by Proposition
       31FRAV, D_11 by Proposition 31FRAW, D_7 by Proposition 31FRAX,
       D_8 by Proposition 31FRAY, D_9 by Proposition 31FRAZ, D_4
       by Proposition 31FRBA, D_5 by Proposition 31FRBB, D_6 by
       Proposition 31FRBC, D_10 by Proposition 31FRBD, B/C rank
       5 by Proposition 31FRBE, and B/C rank 4 by Proposition 31FRBF;
       all m=28 diagnostic rows are closed; this
       submeter is not the global denominator)
       active diagnostic m=29 bridge submeter: 118/118 rows closed
       (B/C ranks 20..30 by Proposition 31FRBG and D_16..D_63 by
       Proposition 31FRBH, B/C rank 19 by Proposition 31FRBI, B/C ranks
       8..16 by Proposition 31FRBJ, and B/C rank 18 by Proposition
       31FRBK, B/C rank 17 by Proposition 31FRBL, and D_15 by
       Proposition 31FRBM, D_14 by Proposition 31FRBN, and D_13 by
       Proposition 31FRBO, D_12 by Proposition 31FRBP, and D_11 by
       Proposition 31FRBQ, D_10 by Proposition 31FRBR, and D_9 by
       Proposition 31FRBS, D_8 by Proposition 31FRBT, D_7 by Proposition
       31FRBU, D_6 by Proposition 31FRBV, D_5 by Proposition 31FRBW,
       D_4 by Proposition 31FRBX, B/C rank 7 by Proposition 31FRBY, and
       B/C rank 6 by Proposition 31FRBZ, B/C rank 5 by Proposition
       31FRCA, B/C rank 2 by Proposition 31FRCB, B/C rank 3 by
       Proposition 31FRCC, and B/C rank 4 by Proposition 31FRCD; all active
       diagnostic m=29 bridge rows are closed; this submeter is not the
       global denominator)
 4/5   conditional integration of CMO-ClassicalTail after the first-hit replay,
       via Corollary 24.17G and Proposition 31FRCH
```

First-hit finite trace closure.  Proposition 24.17D applies Proposition 24.8
at `N_hit(G)=max(45,stable_odd_reach(G)+2)` with

```text
r=20001/10000,   eta=1861/10000.
```

The replay command

```text
python3 -u character_ring_iter/classical_tail_constants.py \
  --first-hit-trace-diagnostic --finite-cutoff 296 --current-max-chain 20
```

returns:

```text
rows tested: 881, trace-covered: 57
  B: ranks 276, 278..295 (19 rows)
  C: ranks 276, 278..295 (19 rows)
  D: ranks 276, 278, 280..295, 297 (19 rows)
  minimum OK margin: D_281 C=279 N_hit=279
    margin=0.66444449761076900114840811256041268244011282766005...
```

These rows no longer require `N_loc` input after `m=20`; Proposition 24.8
propagates their first-hit trace inequality to all later odd exponents.

Finite trace-tail diagnostic.  The Rains-square trace route used for
`C_G>=296` does not by itself discharge the finite-row local onset.  The
replay command

```text
python3 -u character_ring_iter/classical_tail_constants.py \
  --finite-trace-tail-diagnostic --finite-cutoff 296 \
  --current-max-chain 20 --slice-steps 3
```

tests Proposition 24.8 plus Proposition 24.10 with the same
Courteaut-Johansson-Lambert/Rains-square tail inputs used in Proposition
24.16B.  It returns zero trace-covered rows in the first three post-`m=20`
slices:

```text
m=21, odd n=45: rows 86, trace-covered 0
  CJL-source-range: 86
m=22, odd n=47: rows 90, trace-covered 0
  CJL-source-range: 90
m=23, odd n=49: rows 94, trace-covered 0
  CJL-source-range: 94
```

A later-window replay

```text
python3 -u character_ring_iter/classical_tail_constants.py \
  --finite-trace-tail-diagnostic --finite-cutoff 296 \
  --current-max-chain 140 --slice-steps 5
```

also returns zero trace-covered rows at `m=141,...,145`; the rows outside
the source range are blocked by the CJL size hypotheses, and the remaining
source-admissible rows are blocked by the tau upper bound exceeding the
positive-trace lower bound.  Thus the remaining local-onset bucket
must be closed by an effective positive-saddle/local-onset estimate or by a
strictly stronger finite-row tail estimate, not by replaying the high-rank
trace constants unchanged.

Post-`m=29` uniform-supplier scan.  The reproducible diagnostic

```text
python3 -u character_ring_iter/post_m29_uniform_supplier_scan.py
```

starts from the post-`m=29` finite-tail residue after the direct rows and the
top finite trace rows are removed.  After Propositions 31FRDN3, 31FRDN6,
31FRDN9, 31FRDN12, 31FRDW, 31FRDW2, 31FRDW3, and 31FRDW4 remove `B_11`, `B_12`, `C_13`,
`C_14`, `C_15`, `C_16`, `C_17`, and `C_18`, the pre-Chernoff all-later residue has 788 rows.

```text
current post-m29 open finite rows
  count: 788
  B: 13..275, 277 (264); C: 19..275, 277 (258);
  D: 12..275, 277, 279 (266)
```

The corrected fixed-root-angle cap check uses the radial ball inscribed in
the cap, not the radial envelope containing it.  With that correction the
single fixed cap is not a supplier:

```text
delta=2.0: bad=791, worst=C_277 slack=-430472.21
delta=2.3: bad=791, worst=C_277 slack=-419598.41
delta=2.5: bad=791, worst=C_277 slack=-414391.21
delta=2.8: bad=791, worst=C_277 slack=-409338.08
delta=3.0: bad=791, worst=C_277 slack=-407702.12
```

The same scan tested the trace-source supplier.  Replacing the Rains-square
`tau_G` input in the finite-row source range by the source-supported
Courteaut-Johansson-Lambert/Rains Chernoff concentration term without the
additive total-variation floor closes 50 post-m29 high-edge rows with
`r=2001/1000` and `eta=9/20` for B/C:

```text
B: ranks 218, 220, ..., 266 (25 rows)
C: ranks 218, 220, ..., 266 (25 rows)
D: ranks none (0 rows)
```

The source-audited replay

```text
python3 -u character_ring_iter/classical_tail_constants.py \
  --post-m29-chernoff-trace-certificate --finite-cutoff 296
```

checks the one-trace source range, the Chernoff admissibility window, and the
finite-row trace-pushforward comparison.  It closes `0` of the original
`796` post-m29 open rows under the earlier parameter choice; with the checked
high-edge parameters it closes 118 rows of the pre-Chernoff 788-row residue.  After
Propositions 31FRDN3, 31FRDN6, 31FRDN9, 31FRDN12, 31FRDW, 31FRDW2, 31FRDW3, 31FRDW4, 31FRDW5, and 24.17E
remove `B_11`, `B_12`, `C_13`, `C_14`, `C_15`, `C_16`, `C_17`, `C_18`, `C_19`, and the 118 high-edge rows,
the current all-later
residue is

```text
B: 13..217 (205 rows)
C: 20..217 (198 rows)
D: 12..275 plus 277 and 279 (266 rows)
total remaining: 669 rows
```

The proof file uses this checked supplier as Proposition 24.17E.

Current hard submeter for the first post-`m=20` slice:

```text
m=21 bridge closure: 86/86 rows

closed by status-backed exact rows in Proposition 31DR:
  B_2..B_8, B_17..B_22,
  C_2..C_8, C_17..C_22,
  D_4..D_8, D_19, D_34..D_47

closed by the linear D-support theorem in Proposition 31DS:
  D_27..D_33

closed by the D middle lower-bound theorem in Proposition 31DT:
  D_20..D_26

closed by the D low-middle lower-bound theorem in Proposition 31DU:
  D_12..D_18

closed by the D bottom lower-bound theorem in Proposition 31DV:
  D_9..D_11

closed by the B/C boundary lower-bound theorem in Proposition 31DW:
  B_9..B_16,
  C_9..C_16

remaining m=21 rows:
  none
```

Current hard submeter for the next post-`m=21` slice:

```text
m=22 bridge closure: 90/90 rows

candidate rows before imposing 47<N_loc(G):
  B_2..B_23,
  C_2..C_23,
  D_4..D_49

closed by the m=22 B/C boundary lower-bound theorem in Proposition 31DX:
  B_9..B_23,
  C_9..C_23

closed by the m=22 high D linear theorem in Proposition 31DY:
  D_29..D_49

closed by the m=22 middle D lower-bound theorem in Proposition 31DZ:
  D_18..D_28

closed by the m=22 low-middle D lower-bound theorem in Proposition 31EA:
  D_12..D_17

closed by the m=22 low B/C boundary lower-bound theorem in Proposition 31EB:
  B_6..B_8,
  C_6..C_8

closed by the m=22 B5 lower-bound theorem in Proposition 31EC:
  B_5

closed by the m=22 low D lower-bound theorem in Proposition 31ED:
  D_4..D_11

closed by the m=22 C5 lower-bound theorem in Proposition 31EE:
  C_5

closed by the m=22 B4 positive-quadratic lower-bound theorem in Proposition 31EF:
  B_4

closed by the m=22 C4 positive-quadratic lower-bound theorem in Proposition 31EG:
  C_4

closed by the m=22 rank-three B/C positive-quadratic theorem in Proposition 31EH:
  B_3,
  C_3

closed by the m=22 rank-two B/C exact Weyl theorem in Proposition 31EI:
  B_2,
  C_2

remaining m=22 rows:
  none
```

Current hard submeter for the next post-`m=22` slice:

```text
m=23 bridge closure: 94/94 rows

candidate rows before imposing 49<N_loc(G):
  B_2..B_24,
  C_2..C_24,
  D_4..D_51

closed by the m=23 rank-two B/C exact Weyl theorem in Proposition 31EJ:
  B_2,
  C_2

closed by the m=23 high B/C lower-bound theorem in Proposition 31EK:
  B_16..B_24,
  C_16..C_24

closed by the m=23 high D linear theorem in Proposition 31EL:
  D_30..D_51

closed by the m=23 middle D lower-bound theorem in Proposition 31EM:
  D_20..D_29

closed by the m=23 low-middle D lower-bound theorem in Proposition 31EN:
  D_12..D_19

closed by the m=23 bottom D lower-bound theorem in Proposition 31EO:
  D_4..D_11

closed by the m=23 B/C rank-nine-through-fifteen shared exact-window theorem
in Proposition 31EP:
  B_9..B_15,
  C_9..C_15

closed by the m=23 B/C rank-eight shared exact-window theorem in Proposition
31EQ:
  B_8,
  C_8

closed by the m=23 B/C rank-seven shared exact-window theorem in Proposition
31ER:
  B_7,
  C_7

closed by the m=23 B/C rank-six shared exact-window theorem in Proposition
31ES:
  B_6,
  C_6

closed by the m=23 B/C rank-five shared exact-window theorem in Proposition
31ET:
  B_5,
  C_5

closed by the m=23 B/C rank-four shared exact-window theorem in Proposition
31EU:
  B_4,
  C_4

closed by the m=23 B/C rank-three shared exact-window theorem in Proposition
31EV:
  B_3,
  C_3

remaining m=23 rows:
  none
```

Current hard submeter for the next post-`m=23` slice:

```text
m=24 bridge closure: 98/98 rows

candidate rows before imposing 51<N_loc(G):
  B_2..B_25,
  C_2..C_25,
  D_4..D_53

closed by the m=24 rank-two B/C exact Weyl theorem in Proposition 31EW:
  B_2,
  C_2

closed by the m=24 high B/C lower-bound theorem in Proposition 31EX:
  B_16..B_25,
  C_16..C_25

closed by the m=24 high D linear theorem in Proposition 31EY:
  D_32..D_53

closed by the m=24 middle D lower-bound theorem in Proposition 31EZ:
  D_21..D_31

closed by the m=24 low-middle D lower-bound theorem in Proposition 31FA:
  D_13..D_20

closed by the m=24 bottom D lower-bound theorems in Propositions 31FB and
31FD--31FH:
  D_4..D_12

closed by the m=24 B/C rank-fifteen shared exact-window theorem in
Proposition 31FC:
  B_15,
  C_15

closed by the m=24 B/C rank-four-through-fourteen shared exact-window theorem
in Proposition 31FI:
  B_4..B_14,
  C_4..C_14

closed by the m=24 rank-three B/C exact Weyl theorem in Proposition 31FJ:
  B_3,
  C_3

remaining m=24 rows:
  none
```

Current hard submeter for the next post-`m=24` slice:

```text
m=25 bridge closure: 102/102 rows

candidate rows before imposing 53<N_loc(G):
  B_2..B_26,
  C_2..C_26,
  D_4..D_55

closed by the m=25 rank-two B/C exact Weyl theorem in Proposition 31FL:
  B_2,
  C_2

closed by the m=25 rank-three B/C exact Weyl theorem in Proposition 31FM:
  B_3,
  C_3

closed by the m=25 rank-four B/C exact Weyl theorem in Proposition 31FN:
  B_4,
  C_4

closed by the m=25 B/C reused exact-window theorem in Proposition 31FO:
  B_5..B_16,
  C_5..C_16

closed by the m=25 high B/C lower-bound theorem in Proposition 31FP:
  B_17..B_26,
  C_17..C_26

closed by the m=25 bottom D reused determinant-delta theorem in Proposition
31FQ:
  D_4..D_13

closed by the m=25 D lower-bound theorem in Proposition 31FR:
  D_14..D_55

remaining m=25 rows:
  none
```

Historical row-gate record for the earlier `m=15` exact Chain slice:
Corollary 24.17A and
`character_ring_iter/classical_tail_constants.py --middle-slice-plan
--current-max-chain 14 --slice-steps 4` give `62` rows at `m=15`, namely
`B_2..B_16`, `C_2..C_16`, and `D_4..D_35`, before imposing the additional
row condition `33<N_loc(G)`.  Proposition 31A closes the boundary row `D_35`
at this slice by the Pfaffian correction
`m_35(D_35)=m_35(stable)+1`, and Proposition 31B closes the adjacent `D_34`,
`D_33`, `D_32`, and `D_31` rows by the depth-four determinant corrections of Lemma
16B.  Proposition 31C closes `D_30` through `D_27` by the finite Pieri
recurrence of Lemma 16C, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth8-certificate`.
Proposition 31D closes `D_26` through `D_24` by the same recurrence, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth11-certificate`.
Proposition 31E closes `D_23` through `D_20`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth15-certificate`.
Proposition 31F closes `D_19`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth16-certificate`.
Proposition 31G closes `D_18`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth17-certificate`.
Proposition 31H closes `D_17`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth18-certificate`.
Proposition 31I closes `D_16`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth19-certificate`.
Proposition 31J closes `D_15`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth20-certificate`.
Proposition 31K closes `D_14`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth21-certificate`.
Proposition 31L closes `D_13`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth22-certificate`.
Proposition 31M closes `D_12`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth23-certificate`.
Proposition 31N closes `D_11`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth24-certificate`.
Proposition 31O closes `D_10`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth25-certificate`.
Proposition 31P closes `D_9`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth26-certificate`.
Proposition 31Q closes `D_8`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth27-certificate`.
Proposition 31R closes `D_7`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth28-certificate`.
Proposition 31S closes `D_6`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth29-certificate`.
Proposition 31T closes `D_5`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth30-certificate`.
Proposition 31U closes `D_4`, replayed by
`character_ring_iter/classical_tail_constants.py --d-boundary-m15-depth31-certificate`.
Proposition 31V closes `B_16`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank16-certificate`.
Proposition 31W closes `C_16`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank16-certificate`.
Proposition 31X closes `B_15`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 15`.
Proposition 31Y closes `C_15`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 15`.
Proposition 31Z closes `B_14`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 14`.
Proposition 31AA closes `C_14`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 14`.
Proposition 31AB closes `B_13`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 13`.
Proposition 31AC closes `C_13`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 13`.
Proposition 31AD closes `B_12`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 12`.
Proposition 31AE closes `C_12`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 12`.
Proposition 31AF closes `B_11`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 11`.
Proposition 31AG closes `C_11`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 11`.
Proposition 31AH closes `B_10`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 10`.
Proposition 31AI closes `C_10`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 10`.
Proposition 31AJ closes `B_9`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 9`.
Proposition 31AK closes `C_9`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 9`.
Proposition 31AL closes `B_8`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 8`.
Proposition 31AM closes `C_8`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 8`.
Proposition 31AN closes `B_7`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 7`.
Proposition 31AO closes `C_7`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 7`.
Proposition 31AP closes `B_6`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 6`.
Proposition 31AQ closes `C_6`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 6`.
Proposition 31AR closes `B_5`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 5`.
Proposition 31AS closes `C_5`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 5`.
Proposition 31AT closes `B_4`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 4`.
Proposition 31AU closes `C_4`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 4`.
Proposition 31AV closes `B_3`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 3`.
Proposition 31AW closes `C_3`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 3`.
Proposition 31AX closes `B_2`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--b-boundary-m15-rank 2`.
Proposition 31AY closes `C_2`, replayed by the C++ certificate compiled from
`character_ring_iter/classical_boundary_certificate.cpp` and run with
`--c-boundary-m15-rank 2`.
The first pre-Wong slice has no remaining B/C/D row after Proposition 31AY.
Proposition 31AZ closes the next boundary rows `D_37`, `B_17`, and `C_17`
at `m=16`, replayed by the same C++ certificate with
`--d-boundary-rank 37 --chain-m 16`,
`--b-boundary-rank 17 --chain-m 16`, and
`--c-boundary-rank 17 --chain-m 16`.
Proposition 31BA closes the adjacent rows `D_36`, `B_16`, and `C_16`,
replayed by the same C++ certificate with
`--d-boundary-rank 36 --chain-m 16`,
`--b-boundary-rank 16 --chain-m 16`, and
`--c-boundary-rank 16 --chain-m 16`.
Proposition 31BB closes the next adjacent rows `D_35`, `B_15`, and `C_15`,
replayed by the same C++ certificate with
`--d-boundary-rank 35 --chain-m 16`,
`--b-boundary-rank 15 --chain-m 16`, and
`--c-boundary-rank 15 --chain-m 16`.
Proposition 31BC closes the next adjacent rows `D_34`, `B_14`, and `C_14`,
replayed by the same C++ certificate with
`--d-boundary-rank 34 --chain-m 16`,
`--b-boundary-rank 14 --chain-m 16`, and
`--c-boundary-rank 14 --chain-m 16`.
Proposition 31BD closes the next adjacent rows `D_33`, `B_13`, and `C_13`,
replayed by the same C++ certificate with
`--d-boundary-rank 33 --chain-m 16`,
`--b-boundary-rank 13 --chain-m 16`, and
`--c-boundary-rank 13 --chain-m 16`.
Proposition 31BE closes the next adjacent rows `D_32`, `B_12`, and `C_12`,
replayed by the same C++ certificate with
`--d-boundary-rank 32 --chain-m 16`,
`--b-boundary-rank 12 --chain-m 16`, and
`--c-boundary-rank 12 --chain-m 16`.
Proposition 31BF closes the next adjacent rows `D_31`, `B_11`, and `C_11`,
replayed by the same C++ certificate with
`--d-boundary-rank 31 --chain-m 16`,
`--b-boundary-rank 11 --chain-m 16`, and
`--c-boundary-rank 11 --chain-m 16`.
Proposition 31BG closes the next adjacent rows `D_30`, `B_10`, and `C_10`,
replayed by the same C++ certificate with
`--d-boundary-rank 30 --chain-m 16`,
`--b-boundary-rank 10 --chain-m 16`, and
`--c-boundary-rank 10 --chain-m 16`.
Proposition 31BH closes the next adjacent rows `D_29`, `B_9`, and `C_9`,
replayed by the same C++ certificate with
`--d-boundary-rank 29 --chain-m 16`,
`--b-boundary-rank 9 --chain-m 16`, and
`--c-boundary-rank 9 --chain-m 16`.
Proposition 31BI closes the next adjacent rows `D_28`, `B_8`, and `C_8`,
replayed by the same C++ certificate with
`--d-boundary-rank 28 --chain-m 16`,
`--b-boundary-rank 8 --chain-m 16`, and
`--c-boundary-rank 8 --chain-m 16`.
Proposition 31BJ closes the next adjacent rows `D_27`, `B_7`, and `C_7`,
replayed by the same C++ certificate with
`--d-boundary-rank 27 --chain-m 16`,
`--b-boundary-rank 7 --chain-m 16`, and
`--c-boundary-rank 7 --chain-m 16`.
Proposition 31BK closes the next adjacent rows `D_26`, `B_6`, and `C_6`,
replayed by the same C++ certificate with
`--d-boundary-rank 26 --chain-m 16`,
`--b-boundary-rank 6 --chain-m 16`, and
`--c-boundary-rank 6 --chain-m 16`.
Proposition 31BL closes the next adjacent rows `D_25`, `B_5`, and `C_5`,
replayed by the same C++ certificate with
`--d-boundary-rank 25 --chain-m 16`,
`--b-boundary-rank 5 --chain-m 16`, and
`--c-boundary-rank 5 --chain-m 16`.
Proposition 31BM closes the next adjacent rows `D_24`, `B_4`, and `C_4`,
replayed by the same C++ certificate with
`--d-boundary-rank 24 --chain-m 16`,
`--b-boundary-rank 4 --chain-m 16`, and
`--c-boundary-rank 4 --chain-m 16`.
Proposition 31BN closes the next adjacent rows `D_23`, `B_3`, and `C_3`,
replayed by the same C++ certificate with
`--d-boundary-rank 23 --chain-m 16`,
`--b-boundary-rank 3 --chain-m 16`, and
`--c-boundary-rank 3 --chain-m 16`.
Proposition 31BO closes the next adjacent rows `D_22`, `B_2`, and `C_2`,
replayed by the same C++ certificate with
`--d-boundary-rank 22 --chain-m 16`,
`--b-boundary-rank 2 --chain-m 16`, and
`--c-boundary-rank 2 --chain-m 16`.
Proposition 31BP closes the next adjacent rows `D_21` through `D_18`,
replayed by the same C++ certificate with
`--d-boundary-rank 21 --chain-m 16`,
`--d-boundary-rank 20 --chain-m 16`,
`--d-boundary-rank 19 --chain-m 16`, and
`--d-boundary-rank 18 --chain-m 16`.
Proposition 31BQ closes the next adjacent rows `D_17` through `D_14`,
replayed by the same C++ certificate with
`--d-boundary-rank 17 --chain-m 16`,
`--d-boundary-rank 16 --chain-m 16`,
`--d-boundary-rank 15 --chain-m 16`, and
`--d-boundary-rank 14 --chain-m 16`.
Proposition 31BR closes the next adjacent rows `D_13` through `D_10`,
replayed by the same C++ certificate with
`--d-boundary-rank 13 --chain-m 16`,
`--d-boundary-rank 12 --chain-m 16`,
`--d-boundary-rank 11 --chain-m 16`, and
`--d-boundary-rank 10 --chain-m 16`.
Proposition 31BS closes the next adjacent rows `D_9` through `D_6`,
replayed by the same C++ certificate with
`--d-boundary-rank 9 --chain-m 16`,
`--d-boundary-rank 8 --chain-m 16`,
`--d-boundary-rank 7 --chain-m 16`, and
`--d-boundary-rank 6 --chain-m 16`.
Proposition 31BT closes the final rows `D_5` and `D_4`, replayed by the same
C++ certificate with
`--d-boundary-rank 5 --chain-m 16` and
`--d-boundary-rank 4 --chain-m 16`.  Thus the second pre-Wong slice is
reduced from `66` rows to no remaining row before imposing the additional
condition `35<N_loc(G)`.
For the next Chain step, Corollary 24.17A gives `70` pre-Wong rows at
`m=17`, namely `B_2..B_18`, `C_2..C_18`, and `D_4..D_39`, before imposing
the additional row condition `37<N_loc(G)`.  Proposition 31BU closes the
boundary rows `D_39`, `B_18`, and `C_18`, replayed by the same C++
certificate with
`--d-boundary-rank 39 --chain-m 17`,
`--b-boundary-rank 18 --chain-m 17`, and
`--c-boundary-rank 18 --chain-m 17`.
Proposition 31BV closes the adjacent rows `D_38`, `B_17`, and `C_17`,
replayed by the same C++ certificate with
`--d-boundary-rank 38 --chain-m 17`,
`--b-boundary-rank 17 --chain-m 17`, and
`--c-boundary-rank 17 --chain-m 17`.
Proposition 31BW closes the adjacent rows `D_37`, `B_16`, and `C_16`,
replayed by the same C++ certificate with
`--d-boundary-rank 37 --chain-m 17`,
`--b-boundary-rank 16 --chain-m 17`, and
`--c-boundary-rank 16 --chain-m 17`.
Proposition 31BX closes the adjacent rows `D_36`, `B_15`, and `C_15`,
replayed by the same C++ certificate with
`--d-boundary-rank 36 --chain-m 17`,
`--b-boundary-rank 15 --chain-m 17`, and
`--c-boundary-rank 15 --chain-m 17`.
Proposition 31BY closes the adjacent rows `D_35`, `B_14`, and `C_14`,
replayed by the same C++ certificate with
`--d-boundary-rank 35 --chain-m 17`,
`--b-boundary-rank 14 --chain-m 17`, and
`--c-boundary-rank 14 --chain-m 17`.
Proposition 31BZ closes the adjacent rows `D_34`, `B_13`, and `C_13`,
replayed by the same C++ certificate with
`--d-boundary-rank 34 --chain-m 17`,
`--b-boundary-rank 13 --chain-m 17`, and
`--c-boundary-rank 13 --chain-m 17`.
Proposition 31CA closes the adjacent rows `D_33`, `B_12`, and `C_12`,
replayed by the same C++ certificate with
`--d-boundary-rank 33 --chain-m 17`,
`--b-boundary-rank 12 --chain-m 17`, and
`--c-boundary-rank 12 --chain-m 17`.
Proposition 31CB closes the adjacent rows `D_32`, `B_11`, and `C_11`,
replayed by the same C++ certificate with
`--d-boundary-rank 32 --chain-m 17`,
`--b-boundary-rank 11 --chain-m 17`, and
`--c-boundary-rank 11 --chain-m 17`.
Proposition 31CC closes the adjacent rows `D_31`, `B_10`, and `C_10`,
replayed by the same C++ certificate with
`--d-boundary-rank 31 --chain-m 17`,
`--b-boundary-rank 10 --chain-m 17`, and
`--c-boundary-rank 10 --chain-m 17`.
Proposition 31CD closes the adjacent rows `D_30`, `B_9`, and `C_9`,
replayed by the same C++ certificate with
`--d-boundary-rank 30 --chain-m 17`,
`--b-boundary-rank 9 --chain-m 17`, and
`--c-boundary-rank 9 --chain-m 17`.
Proposition 31CE closes the adjacent rows `D_29`, `B_8`, and `C_8`,
replayed by the same C++ certificate with
`--d-boundary-rank 29 --chain-m 17`,
`--b-boundary-rank 8 --chain-m 17`, and
`--c-boundary-rank 8 --chain-m 17`.
Proposition 31CF closes the adjacent rows `D_28`, `B_7`, and `C_7`,
replayed by the same C++ certificate with
`--d-boundary-rank 28 --chain-m 17`,
`--b-boundary-rank 7 --chain-m 17`, and
`--c-boundary-rank 7 --chain-m 17`.
Proposition 31CG closes the adjacent rows `D_27`, `B_6`, and `C_6`,
replayed by the same C++ certificate with
`--d-boundary-rank 27 --chain-m 17`,
`--b-boundary-rank 6 --chain-m 17`, and
`--c-boundary-rank 6 --chain-m 17`.
Proposition 31CH closes the adjacent rows `D_26`, `B_5`, and `C_5`,
replayed by the same C++ certificate with
`--d-boundary-rank 26 --chain-m 17`,
`--b-boundary-rank 5 --chain-m 17`, and
`--c-boundary-rank 5 --chain-m 17`.
Proposition 31CI closes the adjacent rows `D_25`, `B_4`, and `C_4`,
replayed by the same C++ certificate with
`--d-boundary-rank 25 --chain-m 17`,
`--b-boundary-rank 4 --chain-m 17`, and
`--c-boundary-rank 4 --chain-m 17`.
Proposition 31CJ closes the adjacent rows `D_24`, `B_3`, and `C_3`,
replayed by the same C++ certificate with
`--d-boundary-rank 24 --chain-m 17`,
`--b-boundary-rank 3 --chain-m 17`, and
`--c-boundary-rank 3 --chain-m 17`.
Proposition 31CK closes the adjacent rows `D_23`, `B_2`, and `C_2`,
replayed by the same C++ certificate with
`--d-boundary-rank 23 --chain-m 17`,
`--b-boundary-rank 2 --chain-m 17`, and
`--c-boundary-rank 2 --chain-m 17`.
Proposition 31CL closes the adjacent rows `D_22`, `D_21`, and `D_20`,
replayed by the same C++ certificate with
`--d-boundary-rank 22 --chain-m 17`,
`--d-boundary-rank 21 --chain-m 17`, and
`--d-boundary-rank 20 --chain-m 17`.
Proposition 31CM closes the adjacent rows `D_19`, `D_18`, and `D_17`,
replayed by the same C++ certificate with
`--d-boundary-rank 19 --chain-m 17`,
`--d-boundary-rank 18 --chain-m 17`, and
`--d-boundary-rank 17 --chain-m 17`.
Proposition 31CN closes the adjacent rows `D_16`, `D_15`, and `D_14`,
replayed by the same C++ certificate with
`--d-boundary-rank 16 --chain-m 17`,
`--d-boundary-rank 15 --chain-m 17`, and
`--d-boundary-rank 14 --chain-m 17`.
Proposition 31CO closes the adjacent rows `D_13`, `D_12`, `D_11`, and
`D_10`, replayed by the same C++ certificate with
`--d-boundary-rank 13 --chain-m 17`,
`--d-boundary-rank 12 --chain-m 17`,
`--d-boundary-rank 11 --chain-m 17`, and
`--d-boundary-rank 10 --chain-m 17`.
Proposition 31CP closes the adjacent rows `D_9`, `D_8`, `D_7`, and `D_6`,
replayed by the same C++ certificate with
`--d-boundary-rank 9 --chain-m 17`,
`--d-boundary-rank 8 --chain-m 17`,
`--d-boundary-rank 7 --chain-m 17`, and
`--d-boundary-rank 6 --chain-m 17`.
Proposition 31CQ closes the final rows `D_5` and `D_4`, replayed by the same
C++ certificate with `--d-boundary-rank 5 --chain-m 17` and
`--d-boundary-rank 4 --chain-m 17`.  Thus the third pre-Wong slice has no
remaining row before imposing `37<N_loc(G)`.
For the next Chain step, Corollary 24.17A gives `74` pre-Wong rows at
`m=18`, namely `B_2..B_19`, `C_2..C_19`, and `D_4..D_41`, before imposing
the additional row condition `39<N_loc(G)`.  Proposition 31CR closes the
boundary rows `D_41`, `B_19`, and `C_19`, replayed by the same C++
certificate with
`--d-boundary-rank 41 --chain-m 18`,
`--b-boundary-rank 19 --chain-m 18`, and
`--c-boundary-rank 19 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `71` rows `B_2..B_18`,
`C_2..C_18`, and `D_4..D_40`, before imposing `39<N_loc(G)`.
Proposition 31CS closes the adjacent rows `D_40`, `B_18`, and `C_18`,
replayed by the same C++ certificate with
`--d-boundary-rank 40 --chain-m 18`,
`--b-boundary-rank 18 --chain-m 18`, and
`--c-boundary-rank 18 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `68` rows `B_2..B_17`,
`C_2..C_17`, and `D_4..D_39`, before imposing `39<N_loc(G)`.
Proposition 31CT closes the adjacent rows `D_39`, `B_17`, and `C_17`,
replayed by the same C++ certificate with
`--d-boundary-rank 39 --chain-m 18`,
`--b-boundary-rank 17 --chain-m 18`, and
`--c-boundary-rank 17 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `65` rows `B_2..B_16`,
`C_2..C_16`, and `D_4..D_38`, before imposing `39<N_loc(G)`.
Proposition 31CU closes the adjacent rows `D_38`, `B_16`, and `C_16`,
replayed by the same C++ certificate with
`--d-boundary-rank 38 --chain-m 18`,
`--b-boundary-rank 16 --chain-m 18`, and
`--c-boundary-rank 16 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `62` rows `B_2..B_15`,
`C_2..C_15`, and `D_4..D_37`, before imposing `39<N_loc(G)`.
Proposition 31CV closes the adjacent rows `D_37`, `B_15`, and `C_15`,
replayed by the same C++ certificate with
`--d-boundary-rank 37 --chain-m 18`,
`--b-boundary-rank 15 --chain-m 18`, and
`--c-boundary-rank 15 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `59` rows `B_2..B_14`,
`C_2..C_14`, and `D_4..D_36`, before imposing `39<N_loc(G)`.
Proposition 31CW closes the adjacent rows `D_36`, `B_14`, and `C_14`,
replayed by the same C++ certificate with
`--d-boundary-rank 36 --chain-m 18`,
`--b-boundary-rank 14 --chain-m 18`, and
`--c-boundary-rank 14 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `56` rows `B_2..B_13`,
`C_2..C_13`, and `D_4..D_35`, before imposing `39<N_loc(G)`.
Proposition 31CX closes the adjacent rows `D_35`, `B_13`, and `C_13`,
replayed by the same C++ certificate with
`--d-boundary-rank 35 --chain-m 18`,
`--b-boundary-rank 13 --chain-m 18`, and
`--c-boundary-rank 13 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `53` rows `B_2..B_12`,
`C_2..C_12`, and `D_4..D_34`, before imposing `39<N_loc(G)`.
Proposition 31CY closes the adjacent rows `D_34`, `B_12`, and `C_12`,
replayed by the same C++ certificate with
`--d-boundary-rank 34 --chain-m 18`,
`--b-boundary-rank 12 --chain-m 18`, and
`--c-boundary-rank 12 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `50` rows `B_2..B_11`,
`C_2..C_11`, and `D_4..D_33`, before imposing `39<N_loc(G)`.
Proposition 31CZ closes the adjacent rows `D_33`, `B_11`, and `C_11`,
replayed by the same C++ certificate with
`--d-boundary-rank 33 --chain-m 18`,
`--b-boundary-rank 11 --chain-m 18`, and
`--c-boundary-rank 11 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `47` rows `B_2..B_10`,
`C_2..C_10`, and `D_4..D_32`, before imposing `39<N_loc(G)`.
Proposition 31DA closes the adjacent rows `D_32`, `B_10`, and `C_10`,
replayed by the same C++ certificate with
`--d-boundary-rank 32 --chain-m 18`,
`--b-boundary-rank 10 --chain-m 18`, and
`--c-boundary-rank 10 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `44` rows `B_2..B_9`,
`C_2..C_9`, and `D_4..D_31`, before imposing `39<N_loc(G)`.
Proposition 31DB closes the adjacent rows `D_31`, `B_9`, and `C_9`,
replayed by the same C++ certificate with
`--d-boundary-rank 31 --chain-m 18`,
`--b-boundary-rank 9 --chain-m 18`, and
`--c-boundary-rank 9 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `41` rows `B_2..B_8`,
`C_2..C_8`, and `D_4..D_30`, before imposing `39<N_loc(G)`.
Proposition 31DC closes the adjacent rows `D_30`, `B_8`, and `C_8`,
replayed by the same C++ certificate with
`--d-boundary-rank 30 --chain-m 18`,
`--b-boundary-rank 8 --chain-m 18`, and
`--c-boundary-rank 8 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `38` rows `B_2..B_7`,
`C_2..C_7`, and `D_4..D_29`, before imposing `39<N_loc(G)`.
Proposition 31DD closes the adjacent rows `D_29`, `B_7`, and `C_7`,
replayed by the same C++ certificate with
`--d-boundary-rank 29 --chain-m 18`,
`--b-boundary-rank 7 --chain-m 18`, and
`--c-boundary-rank 7 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `35` rows `B_2..B_6`,
`C_2..C_6`, and `D_4..D_28`, before imposing `39<N_loc(G)`.
Proposition 31DE closes the adjacent rows `D_28`, `B_6`, and `C_6`,
replayed by the same C++ certificate with
`--d-boundary-rank 28 --chain-m 18`,
`--b-boundary-rank 6 --chain-m 18`, and
`--c-boundary-rank 6 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `32` rows `B_2..B_5`,
`C_2..C_5`, and `D_4..D_27`, before imposing `39<N_loc(G)`.
Proposition 31DF closes the adjacent rows `D_27`, `B_5`, and `C_5`,
replayed by the same C++ certificate with
`--d-boundary-rank 27 --chain-m 18`,
`--b-boundary-rank 5 --chain-m 18`, and
`--c-boundary-rank 5 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `29` rows `B_2..B_4`,
`C_2..C_4`, and `D_4..D_26`, before imposing `39<N_loc(G)`.
Proposition 31DG closes the adjacent rows `D_26`, `B_4`, and `C_4`,
replayed by the same C++ certificate with
`--d-boundary-rank 26 --chain-m 18`,
`--b-boundary-rank 4 --chain-m 18`, and
`--c-boundary-rank 4 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `26` rows `B_2..B_3`,
`C_2..C_3`, and `D_4..D_25`, before imposing `39<N_loc(G)`.
Proposition 31DH closes the adjacent rows `D_25`, `B_3`, and `C_3`,
replayed by the same C++ certificate with
`--d-boundary-rank 25 --chain-m 18`,
`--b-boundary-rank 3 --chain-m 18`, and
`--c-boundary-rank 3 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `23` rows `B_2`, `C_2`, and
`D_4..D_24`, before imposing `39<N_loc(G)`.
Proposition 31DI closes the adjacent rows `D_24`, `B_2`, and `C_2`,
replayed by the same C++ certificate with
`--d-boundary-rank 24 --chain-m 18`,
`--b-boundary-rank 2 --chain-m 18`, and
`--c-boundary-rank 2 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `20` rows `D_4..D_23`, before
imposing `39<N_loc(G)`.
Proposition 31DJ closes the adjacent rows `D_23`, `D_22`, `D_21`, and
`D_20`, replayed by the same C++ certificate with
`--d-boundary-rank 23 --chain-m 18`,
`--d-boundary-rank 22 --chain-m 18`,
`--d-boundary-rank 21 --chain-m 18`, and
`--d-boundary-rank 20 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `16` rows `D_4..D_19`, before
imposing `39<N_loc(G)`.
Proposition 31DK closes the adjacent rows `D_19`, `D_18`, `D_17`, and
`D_16`, replayed by the same C++ certificate with
`--d-boundary-rank 19 --chain-m 18`,
`--d-boundary-rank 18 --chain-m 18`,
`--d-boundary-rank 17 --chain-m 18`, and
`--d-boundary-rank 16 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `12` rows `D_4..D_15`, before
imposing `39<N_loc(G)`.
Proposition 31DL closes the adjacent rows `D_15`, `D_14`, `D_13`, and
`D_12`, replayed by the same C++ certificate with
`--d-boundary-rank 15 --chain-m 18`,
`--d-boundary-rank 14 --chain-m 18`,
`--d-boundary-rank 13 --chain-m 18`, and
`--d-boundary-rank 12 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `8` rows `D_4..D_11`, before
imposing `39<N_loc(G)`.
Proposition 31DM closes the adjacent rows `D_11`, `D_10`, `D_9`, and
`D_8`, replayed by the same C++ certificate with
`--d-boundary-rank 11 --chain-m 18`,
`--d-boundary-rank 10 --chain-m 18`,
`--d-boundary-rank 9 --chain-m 18`, and
`--d-boundary-rank 8 --chain-m 18`.  Thus this fourth pre-Wong slice is
reduced from `74` rows to the remaining `4` rows `D_4..D_7`, before
imposing `39<N_loc(G)`.
Proposition 31DN closes the adjacent rows `D_7`, `D_6`, `D_5`, and
`D_4`, replayed by the same C++ certificate with
`--d-boundary-rank 7 --chain-m 18`,
`--d-boundary-rank 6 --chain-m 18`,
`--d-boundary-rank 5 --chain-m 18`, and
`--d-boundary-rank 4 --chain-m 18`.  Thus this fourth pre-Wong slice has
no remaining row before imposing `39<N_loc(G)`.
Proposition 31DP closes the fifth pre-Wong slice at `m=19`, replayed by
the same C++ certificate with
`--b-boundary-range 2 2 --chain-m 19`,
`--b-boundary-range 3 20 --chain-m 19`,
`--c-boundary-range 2 2 --chain-m 19`,
`--c-boundary-range 3 20 --chain-m 19`,
`--d-boundary-range 4 8 --chain-m 19`,
`--d-boundary-range 9 12 --chain-m 19`, and
`--d-boundary-range 13 43 --chain-m 19`.  Thus this fifth pre-Wong slice has
no remaining B/C/D row before imposing `41<N_loc(G)`.
Proposition 31DQ closes the sixth pre-Wong slice at `m=20`, replayed by the
same C++ certificate with `--chain-m 20` and `--no-predecessor-cache` over
the row set

```text
B_2,...,B_21,
C_2,...,C_21,
D_4,...,D_45.
```

The `/tmp/classical_m20_logs/*.log` parser returned `completed=82 expected=82`
with an empty missing-row list, and the error-marker grep returned no line.
Thus this sixth pre-Wong slice has no remaining B/C/D row before imposing
`43<N_loc(G)`.

The remaining analytic source chain is:

```text
NONABELIAN.md
  "SO(N) universal closure -- structural tail-dominance theorem"
  "ALL SO(N) Chain Inequality closure -- moment-stabilization theorem"
  "Wong-theorem application (Gap X1 precise statement)"
  "Closing Gap X1: Laplace tail dominance via Wong-Erdelyi"
```

Final-proof obligation: prove the remaining finite-row post-m=29 all-later
tail isolated as `CMO-PostM29AllLaterTail` in
`CLASSICAL_M3_ONE_PROOF.md`:

```text
Q3^G(n) >= 0 for every odd n > N_hit^{(29)}(G) over the
non-trace-covered finite rows after the m=29 exact bridge, after removing
the B_2/C_2 rows closed by Proposition 31FRCI and the B_3/C_3 rows closed
by Proposition 31FRCJ, the D_4 row closed by Proposition 31FRCK, and the
B_4/C_4 rows closed by Proposition 31FRCL, and the D_5 row closed by
Proposition 31FRCM, the C_5 row closed by Proposition 31FRCN, and the B_5
row closed by Proposition 31FRCO, and the C_6 row closed by Proposition
31FRCP, the D_6/B_6 rows closed by Proposition 31FRCR, and the
D_7/B_7/C_7 rows closed by Proposition 31FRCS, and the C_8 row closed by
Proposition 31FRCT, the D_8 row closed by Proposition 31FRCV, the
B_8/C_9/D_9 rows closed by Proposition 31FRCX, the C_10 row closed by
Proposition 31FRCY, the B_9 row closed by Proposition 31FRDB, the C_11
row closed by Proposition 31FRDE, the C_12 row closed by Proposition
31FRDI, the B_10 row closed by Proposition 31FRDN, the B_11 row closed by
Proposition 31FRDN3, the B_12 row closed by Proposition 31FRDN6, the C_13
row closed by Proposition 31FRDN9, the C_14 row closed by Proposition
31FRDN12, the C_15 row closed by Proposition 31FRDW, the C_16 row closed by
Proposition 31FRDW2, the C_17 row closed by Proposition 31FRDW3, the C_18
row closed by Proposition 31FRDW4, the D_10 row closed by Proposition 31FRDR, and the D_11
row closed by Proposition 31FRDV.
```

The previous local-onset candidate was

```text
N_loc(G) <= max(63, stable_odd_reach(G)+2)
```

for every finite residue row of Lemma 24.17B except the 57 rows listed in
Proposition 24.17D.
Proposition 31FRCG rules out that local-half-leading route already in the
`B_2` and `C_2` rows at `n=63`.  The first-hit replacement is certified by
Proposition 24.17F for the original 796-row first-hit replay.  After
Propositions 31FRDN3, 31FRDN6, 31FRDW6, 31FRDN9, 31FRDN12, 31FRDW,
31FRDW2, 31FRDW3, 31FRDW4, 31FRDW5, 31FRDW7, 31FRDW8, 31FRDW9,
31FRDW10, 31FRDW11, 31FRDW12, 31FRDW13, 31FRDW14, 31FRDW15, 31FRDW16,
31FRDW17, 31FRDW18, 31FRDW19, 31FRDW20, 31FRDW21, 31FRDW22,
31FRDW23, 31FRDW24, 31FRDW25, 31FRDW26, 31FRDW27, 31FRDW28, 31FRDW29, 31FRDW30, 31FRDW31, 31FRDW32, 31FRDW33, 31FRDW34, 31FRDW35, 31FRDW36, 31FRDW37, 31FRDW38, 31FRDW39, 31FRDW40, 31FRDW41, 31FRDW42, 31FRDW43, 31FRDW44, 31FRDW45, 31FRDW46, 31FRDW47, 31FRDW48, 31FRDW49, 31FRDW50, 31FRDW51, 31FRDW52, 31FRDW53, 31FRDW54, 31FRDW55, 31FRDW56, 31FRDW57, 31FRDW58, 31FRDW59, 31FRDW60, 31FRDW61, 31FRDW62, 31FRDW63, 31FRDW64, 31FRDW65, 31FRDW66, 31FRDW67, 31FRDW68, 31FRDW69, 31FRDW70, 31FRDW71, 31FRDW72, 31FRDW73, 31FRDW74, 31FRDW75, 31FRDW76, 31FRDW77, 31FRDW78, 31FRDW79, 31FRDW80, 31FRDW81, 31FRDW82, 31FRDW83, 31FRDW84, 31FRDW85, 31FRDW86, 31FRDW87, 31FRDW88, 31FRDW89, 31FRDW90, 31FRDW91, and 31FRDW92 remove
`B_11`, `B_12`, `B_13`, `C_13`, `C_14`, `C_15`, `C_16`, `C_17`, `C_18`,
`C_19`, `D_12`, `D_13`, `D_14`, `D_15`, `D_16`, `D_17`, `D_18`, `D_19`, `D_20`, `D_21`, `D_22`, `D_23`, `D_24`, `D_25`, `D_26`, `D_27`, `D_28`, `D_29`, `D_30`, `D_31`, `D_32`, `D_33`, `D_34`, `D_35`, `D_36`, `D_37`, `D_38`, `D_39`, `D_40`, `D_41`, `D_42`, `D_43`, `D_44`, `D_45`, `D_46`, `D_47`, `D_48`, `D_49`, `D_50`, `D_51`, `D_52`, `D_53`, `D_54`, `D_55`, `D_56`, `D_57`, `D_58`, `D_59`, `D_60`, `D_61`, `D_62`, `D_63`, `D_64`, `D_65`, `D_66`, `D_67`, `D_68`, `D_69`, `D_70`, `D_71`, `D_72`, `D_73`, `D_74`, `D_75`, `D_76`, `D_77`, `D_78`, `D_79`, `D_80`, `D_81`, `D_82`, `D_83`, `D_84`, `D_85`, `D_86`, `D_87`, `D_88`, `D_89`, `D_90`, `D_91`, `D_92`, `D_93`, `D_94`, `D_95`, `D_96`, and `D_97`, and Proposition 24.17E removes
`B_b,C_b` for
`218<=b<=275` or `b=277`.  Paper Propositions `prop:post29-d-interval` and
`prop:post29-b-unitary-square` remove the D interval and `B_124..B_217`,
respectively, while `prop:post29-b-tilted-mgf` and the proved correction
prefix remove `B_62..B_123`, and `prop:post29-c-tilted-mgf` with the same
prefix removes `C_62..C_217`; `prop:post29-bc-layered-mgf` removes
`B_14..B_61,C_20..C_61`, with exact Toeplitz--Hankel MGFs on `B_14,B_15`.
Thus no all-later row remains.  Proposition 24.17F alone does not prove the
strict later odd exponents; the layered-MGF theorem supplies them.  The earlier
extracted low rows are closed by Propositions 31FRCI,
31FRCJ, 31FRCK, 31FRCL, 31FRCM, 31FRCN, 31FRCO, 31FRCP, 31FRCR, and
31FRCS, 31FRCT, 31FRCV, 31FRCX, 31FRCY, 31FRDB, 31FRDE, 31FRDI, 31FRDN,
31FRDN3, 31FRDN6, 31FRDN9, 31FRDN12, 31FRDR, 31FRDV, 31FRDW, 31FRDW2,
31FRDW3, 31FRDW4, 31FRDW5, 31FRDW6, 31FRDW7, 31FRDW8, 31FRDW9,
31FRDW10, 31FRDW11, 31FRDW12, 31FRDW13, 31FRDW14, 31FRDW15, 31FRDW16,
31FRDW17, 31FRDW18, 31FRDW19, 31FRDW20, 31FRDW21, 31FRDW22,
31FRDW23, 31FRDW24, 31FRDW25, 31FRDW26, 31FRDW27, 31FRDW28, 31FRDW29, 31FRDW30, 31FRDW31, 31FRDW32, 31FRDW33, 31FRDW34, 31FRDW35, 31FRDW36, 31FRDW37, 31FRDW38, 31FRDW39, 31FRDW40, 31FRDW41, 31FRDW42, 31FRDW43, 31FRDW44, 31FRDW45, 31FRDW46, 31FRDW47, 31FRDW48, 31FRDW49, 31FRDW50, 31FRDW51, 31FRDW52, 31FRDW53, 31FRDW54, 31FRDW55, 31FRDW56, 31FRDW57, 31FRDW58, 31FRDW59, 31FRDW60, 31FRDW61, 31FRDW62, 31FRDW63, 31FRDW64, 31FRDW65, 31FRDW66, 31FRDW67, 31FRDW68, 31FRDW69, 31FRDW70, 31FRDW71, 31FRDW72, 31FRDW73, 31FRDW74, 31FRDW75, 31FRDW76, 31FRDW77, 31FRDW78, 31FRDW79, 31FRDW80, 31FRDW81, 31FRDW82, 31FRDW83, 31FRDW84, 31FRDW85, 31FRDW86, 31FRDW87, 31FRDW88, 31FRDW89, 31FRDW90, 31FRDW91, and 31FRDW92:

```text
B_b,C_b:
  N_loc <= 63       for 2 <= b <= 31,
  N_loc <= 2b+1     for 32 <= b <= 295;

D_b:
  N_loc <= 63       for 4 <= b <= 73,
  N_loc <= b-1      for even b, 80 <= b <= 294,
  N_loc <= b-2      for odd b, 81 <= b <= 295 and b=297.
```

The finite exact bridge part has been extracted into numbered propositions:

1. no unresolved odd `n=31` or `n=33` B/C/D row remains;
2. the `m=15` pre-Wong slice has no remaining B/C/D row after
   Proposition 31AY;
3. the `m=16` pre-Wong slice has no remaining row before the
   `35<N_loc(G)` gate after Proposition 31BT;
4. the `m=17` pre-Wong slice has no remaining row before the
   `37<N_loc(G)` gate after Proposition 31CQ;
5. the `m=18` pre-Wong slice has no remaining row before the
   `39<N_loc(G)` gate after Proposition 31DN.
6. the `m=19` pre-Wong slice has no remaining B/C/D row before the
   `41<N_loc(G)` gate after Proposition 31DP.
7. the `m=20` pre-Wong slice has no remaining B/C/D row before the
   `43<N_loc(G)` gate after Proposition 31DQ.
8. the `m=21` pre-Wong slice has no remaining B/C/D row before the
   `45<N_loc(G)` gate after Propositions 31DR--31DW.
9. the `m=22` pre-Wong slice has no remaining B/C/D row before the
   `47<N_loc(G)` gate after Propositions 31DX--31EI.
10. the `m=23` pre-Wong slice has no remaining B/C/D row before the
    `49<N_loc(G)` gate after Propositions 31EJ--31EV.
11. the `m=24` pre-Wong slice has no remaining B/C/D row before the
    `51<N_loc(G)` gate after Propositions 31EW--31FJ.
12. the `m=25` pre-Wong slice has no remaining B/C/D row before the
    `53<N_loc(G)` gate after Propositions 31FL--31FR.
13. the `m=26` diagnostic slice has no remaining B/C/D row before the
    `55<N_loc(G)` gate after Propositions 31FRD--31FRS.
14. the `m=27` diagnostic slice has no remaining B/C/D row before the
    `57<N_loc(G)` gate after Propositions 31FRT--31FRAI.
15. the `m=28` diagnostic slice has no remaining B/C/D row before the
    `59<N_loc(G)` gate after Propositions 31FRAJ--31FRBF.
16. the `m=29` diagnostic slice has no remaining B/C/D row before the
    `61<N_loc(G)` gate after Propositions 31FRBG--31FRCD.

Proposition 31FRCF records the now-obstructed conditional reduction through
`CMO-LocalOnsetBound29`.  Proposition 31FRCH is now conditional:
`CMO-PostM29AllLaterTail`, Proposition 24.17F at first hit, Proposition
24.17D, and these extracted row bridges discharge `CMO-ClassicalTail`.
Propositions 31FRCI,
31FRCJ, 31FRCK,
31FRCL, 31FRCM, 31FRCN, 31FRCO, 31FRCP, 31FRCR, 31FRCS, 31FRCT, and
31FRCV, 31FRCX, 31FRCY, 31FRDB, 31FRDE, 31FRDI, 31FRDN, 31FRDN3,
31FRDN6, 31FRDN9, 31FRDN12, 31FRDR, 31FRDV, 31FRDW, 31FRDW2, 31FRDW3,
31FRDW4, 31FRDW5, 31FRDW6, 31FRDW7, 31FRDW8, 31FRDW9, 31FRDW10,
31FRDW11, 31FRDW12, 31FRDW13, 31FRDW14, 31FRDW15, 31FRDW16, 31FRDW17,
31FRDW18, 31FRDW19, 31FRDW20, 31FRDW21, 31FRDW22, 31FRDW23,
31FRDW24, 31FRDW25, 31FRDW26, 31FRDW27, 31FRDW28, 31FRDW29, 31FRDW30, 31FRDW31, 31FRDW32, 31FRDW33, 31FRDW34, 31FRDW35, 31FRDW36, 31FRDW37, 31FRDW38, 31FRDW39, 31FRDW40, 31FRDW41, 31FRDW42, 31FRDW43, 31FRDW44, 31FRDW45, 31FRDW46, 31FRDW47, 31FRDW48, 31FRDW49, 31FRDW50, 31FRDW51, 31FRDW52, 31FRDW53, 31FRDW54, 31FRDW55, and 31FRDW56 supply the extracted direct
post-m29 finite-tail rows, `B_2`, `C_2`, `B_3`, `C_3`, `D_4`, `B_4`,
`C_4`, `D_5`, `C_5`, `B_5`, `C_6`, `D_6`, `B_6`, `D_7`, `B_7`, `C_7`,
`C_8`, `D_8`, `B_8`, `C_9`, `D_9`, `C_10`, `B_9`, `C_11`, `C_12`,
`B_10`, `D_10`, `D_11`, `B_11`, `B_12`, `B_13`, `C_13`, `C_14`,
`C_15`, `C_16`, `C_17`, `C_18`, `C_19`, `D_12`, `D_13`, `D_14`, `D_15`,
`D_16`, `D_17`, `D_18`, `D_19`, `D_20`, `D_21`, `D_22`, `D_23`, `D_24`, `D_25`, `D_26`, `D_27`, `D_28`, `D_29`, `D_30`, `D_31`, `D_32`, `D_33`, `D_34`, `D_35`, `D_36`, `D_37`, `D_38`, `D_39`, `D_40`, `D_41`, `D_42`, `D_43`, `D_44`, `D_45`, `D_46`, `D_47`, `D_48`, `D_49`, `D_50`, `D_51`, `D_52`, `D_53`, `D_54`, `D_55`, `D_56`, `D_57`, `D_58`, `D_59`, `D_60`, and `D_61`.

## Proposition 4 (small-rank classical bridge slot)

The historical diary records exact small-rank bridge computations for
`SO(3)`, `SO(5)`, and related low-rank overlaps.  The source locations include:

```text
NONABELIAN.md
  "SO(N) universal closure -- structural tail-dominance theorem"
  "ALL SO(N) Chain Inequality closure -- moment-stabilization theorem"
  "Sp(n) universal closure via MGF coincidence with SO(infinity)"
```

Final-proof obligation: replay or cite the exact integer computations in a
checker-style file, as was done for `verify_sun_repair_certificates.py` and
the exceptional direct-chain certificates.

## Corollary 5 (current extraction target)

The next integration target is:

```text
Classical-m3-one standalone theorem note =
  Proposition 1 (classical degree table closed in CLASSICAL_M3_ONE_PROOF.md)
  + Proposition 2 (closed in CLASSICAL_M3_ONE_PROOF.md)
  + Proposition 3 (SO/Sp stable equality, stable Chain, stable-rank
    reduction, Weyl torus formula, and Wong theorem source closed in
    CLASSICAL_M3_ONE_PROOF.md; identity-saddle leading term and negative
    region bound, positive-saddle tail threshold, classical root-data and
    negative-base constants, high-rank middle-window diagnostic,
    Macdonald-Mehta leading coefficient, leading-edge dominance reduction,
    fixed local-cap obstruction, nonlocal two-sided pushforward tail
    criterion, defining-trace deviation reduction, parameterized high-rank
    trace-deviation target, stabilized tau-tail estimate, exact Chain bridge
    through odd n=23, exact B/C/D bridge through odd n=31 closed there too)
  + Proposition 4.
```

This will put the SO/Sp side of the source map in the same format as the
SU(N) repair and the exceptional direct-chain certificates.
