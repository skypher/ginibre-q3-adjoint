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
| Theorems G3X, G3O, decomposition (1) | `g3x_repro.py`, `g3x_check.py`, `g3x_region_check.py` (args: comma list of r, e.g. `26,27,...,90`) |
| Direct three-factor split evaluator | `split3.py`, `split3_selftest.py` |
| G0 branch scan, Conjecture W (real-rooted windows) | `g0_scan.py`, `window_explore.py`, `window_realrooted.py`, `window_f.py` |
| GFM3 (kernel factors `Q_t`) and controls | `wformula_check.py`, `kernel_identity.py`, `recip_test.py`, `recip_bounds.py`, `recip_quartic.py`, `general_row_words.py` (args: seed, count), `general_row_words_ctrl.py` |
| GFM3 stress screen (clustered labels, minimal e) | `gfm3_hard.py` (args: seed, count) |
| What W uses: ULC vs log-concave, signed Newton, log-concave weights | `window_ulc.py`, `logconcave_test.py` |
| Phase-polygon form of W: Theorem WS, pass criterion, G0 coverage | `polygon_lemma.py`, `window_sweep.py`, `window_passes.py`, `window_cover.py` |
| Newton certificates for W (C <= 2) and their failure for C >= 3 and for (E) | `window_lp.py` (args: C, maxdeg, t-list, extension), `e_lp.py` (args: maxdeg, t-list), `window_gf.py` |
| gamma <= N branch: scan, mirror form | `gammaN_branch_scan.py`, `mirror_check.py` |
| FM-MECH23 reproducer; long-sweep energy bound LE | `mech23_repro.py`, `uncovered_slack.py`, `long_energy.py`, `long_energy2.py`, `long_energy3.py` |
| Constant-coefficient (Sonin) model of W and (E) | `const_coeff_check.py` |
| W = (i) long discriminant + (ii) D-window mean | `sonin_route.py`, `sonin_route2.py` |
| gamma <= N branch vs the cone of OL, W, (E) forms | `basis_lp.py`, `basis_lp2.py` |
| Theorem G0E (G0 branch from (E)); window-(E) on gamma <= N | `w_from_e.py`, `window_E.py` |
| Merge move, groupings, four-wedge correction | `merge_check.py`, `merge_best.py`, `cross_telescope.py` |
| Four-factor double merge census | `merge4.py` |
| Theorem G0B (G0 branch in a growing band) | `g0b_repro.py`, `g0b_check.py` |
| Theorem LD (long discriminant (i), all real-rooted rows); residual of W | `ld_repro.py`, `ld_check.py`, `w_residual.py` |
| (E) at q = 1, gaps 2 and 3 (FM-SEC17) | `e_q1_repro.py`, `e_q1_route_kill.py` |
| QGSI | `qgsi.py`, `qscan.py` (these need `slotHWg.py` and python-flint), `qscan_*.results` |

Files named `*_repro.py` are the agents' printed vetting code, rerun by the main
agent.  The other scripts are the main agent's independent checks.
