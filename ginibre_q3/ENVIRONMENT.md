# Validated replay environment

This file records the software and hardware environment in which the July 2026
publication checks were rebuilt. It is a validated configuration, not a claim
that the sources require exactly these versions.

## Software

- Ubuntu 24.04.4 LTS (x86-64), Linux 6.8.0-90-generic
- GNU g++ 13.3.0
- GNU Make 4.3
- Python 3.12.3
- GMP 6.3.0 (`libgmp-dev` Ubuntu package
  `2:6.3.0+dfsg-2ubuntu6.1`)
- MPFR 4.2.1 (`libmpfr6` Ubuntu package `4.2.1-1build1.1`)
- GNU OpenMP runtime (`libgomp1` Ubuntu package
  `14.2.0-4ubuntu2~24.04.1`)
- pdfTeX 3.141592653-2.6-1.40.25, TeX Live 2023/Debian
- Git 2.43.0

The Makefiles contain the authoritative compiler and linker flags. In
particular, the active arithmetic programs use C++17 or C++20 as stated by
their targets, link GMP and MPFR explicitly, and enable OpenMP only for the
parallel stages.

## Validation host

- AMD EPYC 7551, 64 logical CPUs
- 125 GiB RAM

The full Parts I--II replay is documented for the separate machine-C resource
envelope in `REPLAY.md`; this `optimus` host record documents the environment
used for the publication preflight, manuscript rebuild, and the Part III
aggregate replay. Thread limits and per-stage timeout requirements remain
those in `REPLAY.md` and `FULL_Q3_DISTRIBUTED_REPLAY.md`.

## Recheck commands

```text
make -C ginibre_q3 publication-preflight
make -C ginibre_q3 clean-room-replay
make -C ginibre_q3 full-q3-extension
```

All load-bearing source identities are controlled by
`replay_sources.sha256` and
`certificates/full_q3/full_q3_source_manifest.sha256`; the prose version list
in this file does not replace those byte-level checks.
