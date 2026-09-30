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
| QGSI | `qgsi.py`, `qscan.py` (these need `slotHWg.py` and python-flint), `qscan_*.results` |

Files named `*_repro.py` are the agents' printed vetting code, rerun by the main
agent.  The other scripts are the main agent's independent checks.
