# Verifiers for note item (39), 2026-09-30

Scripts behind item (39) of `../../SU2_FUNDAMENTAL_MINUS_REDUCTION_2026_09_26.md`
(the note cites them as `mech/...`, `inj/...`, `xmn/...`).  Run each one from
this directory, because several `exec` a sibling file.

| Result | Scripts |
|---|---|
| Common `phi_r` from Catalan moments | `n0check.py` |
| h_1 absorption, N0 kill | `n0check.py`, `n0scan.py` |
| R3-S1 | `kraw_check.py`, `kraw_sym.py`, `kraw_bdry.py` |
| R3-TWO | `r3defs.py` (lines 1-37 of `../verify_r3_twopower_cert.py`), `r3tail_check.py`, `r3tail_check2.py` |
| One-label factorization and Gram rows | `gram_repro.py`, `gram_link.py`, `gram_ext.py`, `gramrows.py` (`gramrows7.log`), `zj_repro.py`, `zj_print.py` |
| Theorem OL | `onelabel_repro.py`, `onelabel_indep.py` |
| Two-label reduction (E), strip | `twolabel_repro.py`, `twolabel_indep.py` |
| Theorem RF | `rootfree_repro.py`, `rootfree_indep.py` |
| (PE), two-step splitting | `pe_repro.py`, `pe_scan.py`, `split_repro.py` |
| Norm-chord kill | `normchord.py` |
| (E) screen | `escreen.py` (`escreen2.log`, `escreen_60_150.results`) |
| FM3 kernel screen | `fm3kern.py`, `fm3screen.py` (`fm3screen*.log`, `fm3screen_*.json`) |
| Theorem THREE | `parity_check.py`, `three_repro.py` |
| EM2 at r >= 3, descents | `em2r.py`, `descent_r.py` |
| Theorem LS, LR (large suffix / large r+a) | `radial_repro.py`, `ls_check.py` |
| Theorem LR4 (additive quartic cutoff) | `lr4_repro.py`, `lr4_check.py` |
| Theorem BW (balanced-wedge slice) | `wedge_repro.py` |
| Coverage map and cancellation ratios | `band.py`, `band2.py` |
| Theorem LL (three factors, large labels) | `ll_repro.py`, `ll_check.py` |
| Three-factor closed form, Theorem T3R | `closed3_repro.py`, `closed3_check.py` |
| Theorem DS (diagonal strip) | `ds_repro.py`, `ds_check.py` |
| Theorem AS (adjacent strips) | `as_repro.py`, `as_check.py` |
| Four factors, Theorem LL4 | `four_repro.py`, `four_explore.py`, `ll4_check.py` |
| Theorem LLm, LLm-S (spread labels, any m, with hat S) | `llm_check.py`, `llms_check.py` |
| Offset induction, inequality (P) | `pinduct_repro.py` |
| Theorem G3 (positive next-row minors) | `g3_repro.py`, `g3_check.py` |
| QGSI | `qgsi.py`, `qscan.py` (these need `slotHWg.py` and python-flint), `qscan_*.results` |

Files named `*_repro.py` are the agents' printed vetting code, rerun by the main
agent.  The other scripts are the main agent's independent checks.
