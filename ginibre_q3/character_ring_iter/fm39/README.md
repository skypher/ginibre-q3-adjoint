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
| Theorem G0E4 (four factors from (E) on the support region) | `g0e4_check.py` |
| Row class of (E) (log-concave weights); Fourier representation of D and T | `e_class.py`, `fourier_rep.py` |
| (E) LD energy-drop region; (E) residual; psi-curve form of (E); (ii) at e=1 and plateau | `e_ld_region_repro.py`, `e_residual.py`, `psi_curve.py`, `psi_split.py`, `ii_e1_plateau_repro.py` |
| W residual with all criteria | `w_residual2.py`, `w_residual3.py` |
| FM-MECH25 (centered windows, metric chord bound); equal-label sector; W complete on grid; sampled exhaustiveness | `mech25_repro.py`, `eqlabel_check.py`, `w_residual4.py`, `w_residual_sample.py` |
| Three labels with hat S (M3): closed forms, support band; (E) with metric chords | `sp3_repro.py`, `sp3_band_check.py`, `e_metric.py` |
| Theorem EQ1 ((E) at q=1 from OL); multiplier identity for (E); OL for cyclotomic multiples; residual | `e_q1_ol.py`, `e_ol_lp.py`, `e_multiplier_identity.py`, `e_multiplier_symbolic.py`, `ol_cyclotomic.py`, `e_via_olcyc.py`, `e_residual2.py` |
| Recursive LP for (E) (kill) | `e_recursive_lp.py` |
| Binary-form method for (E); (E) complete on a,e <= 40 | `e_q2_quadratic.py`, `e_binary_form.py`, `e_residual3.py` |
| (E) coverage scans beyond the grid | `e_residual_sample.py`, `e_residual_sample2.py`, `e_residual_rows.py` |
| (E): exact certification a,e<=80; uniform-regions-only residual (only q = 2 left) | `e_cert_rows80.py`, `e_uniform_residual.py`, `e_uniform_sample.py` |
| (E) uniform criteria incl. the q = 2, W >= 0 case: residual 0 | `e_uniform_residual2.py`, `e_uniform_sample2.py` |
| (E): FM-SEC41 family is short-psi-arc; extended menu on extreme and proportional rows (exact) | `psi_arc_family.py`, `e_extreme_residual.py` |
| (E): fixed-e continuum limit, limiting menu coverage, the j = N/2 line via t = 2V | `e_continuum_limit.py`, `e_continuum_cover.py`, `e_center_mt.py` |
| (E) two-case split: sweep >= pi pairs all closed by the single-metric bound (Conjecture HT); HT margin map; where t = 2V fails; killed HT routes | `e_sweep_split.py`, `e_sweep_split_t.py`, `e_ht_margin.py`, `e_ht_2v_fail.py`, `e_spiral_probe.py`, `e_ht_monotone.py` |
| gamma <= N: FM-SEC45 certificate census; direct single-metric test (HT3) census | `gammaN_cert_lib.py`, `gammaN_cert_census.py`, `gammaN_metric_census.py` |
| (E)/(HT): exact (m, n)-coordinate identities; the permanent form (HT*); drop lemma; rigorous per-pair chain (b); two-region plan (I)+(O); true HT margin via the permanent | `e_ht_mn.py`, `e_ht_drop.py`, `e_ht_sinsum.py`, `e_ht_regions.py`, `e_ht_perm_margin.py` |
| M5 sectors: MP_2 Kostka weights = phi_2 on H-only words (any number of factors); FM-MECH28's exact evaluator | `mp2_kostka_check.py`, `mech28_eval.py` |
| M5 Hypothesis B (FM-MECH30 reproducer); randomized exact FM3 screen across levels | `mech30_hypB_repro.py`, `fm3_random_screen.py` |
| q-Gaussian deformation (FM-SEC77); two-parameter braided family screen, Fock evaluator and Bernstein tests; graph-product screen; component-model check; tightness atlas (FM-SEC83); H_AC audit (FM-MECH33); H_AC at n = sum - 4 (FM-MECH34); H_AC at d = 3 partial (FM-MECH35); H_AC dictionary census (FM-SEC89); H_AC on the {1,2} sector (FM-SEC88); canonical square-root kill; B counterexample (FM-SEC68); B for repeated odd labels (FM-MECH31); level-1 two-hat certificates a = 2 (FM-SEC66) and a = 3..6 (FM-SEC79); B restricted families | `sec77_qgauss_repro.py`, `qs_family_screen.py`, `qs_fock.py`, `qs_bernstein.py`, `qs_bernstein_fock.py`, `graph_product_screen.py`, `component_model_check.py`, `sec87_even_gf_repro.py`, `sec94_charging_repro.py`, `sec96_linear_injection_repro.py`, `sec90_braided_repro.py`, `crossing_set_restriction_kill.py`, `sec83_atlas_repro.py`, `mech33_hAC_repro.py`, `mech34_hAC_d2_repro.py`, `mech35_hAC_d3_repro.py`, `mech35_hAC_d3_budget.py`, `mech35_hAC_d3_witnesses.py`, `sec89_hAC_census_repro.py`, `sec88_hAC_12sector_repro.py`, `hAC_sqrt_test.py`, `sec86_two_hat_a6_repro.py`, `sec68_hypB_counterexample.py`, `mech31_hypB_oddlabels_repro.py`, `sec66_two_hat_a2_repro.py`, `hypB_lite.py`, `hypB_partition.py`, `hypB_necessary.py` |
| QGSI | `qgsi.py`, `qscan.py` (these need `slotHWg.py` and python-flint), `qscan_*.results` |

Files named `*_repro.py` are the agents' printed vetting code, rerun by the main
agent.  The other scripts are the main agent's independent checks.
