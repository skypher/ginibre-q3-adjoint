REPLAY_THREADS ?= $(shell nproc)
# Each reverse-Pieri task carries a large GMP frontier.  The residual replay
# exposes at most four independent tasks, so using every logical CPU only
# multiplies peak memory without increasing useful parallelism.  Keep the
# reproducible default memory-safe; an audited host may override it explicitly.
FULL_Q3_REMAINING_THREADS ?= 2
FULL_Q3_BOUNDED_THREADS ?= $(shell nproc)

.PHONY: clean-room-replay full-q3-document-audit full-q3-artifact-audit full-q3-stable-audit full-q3-su2-audit full-q3-su3-audit full-q3-su45-audit full-q3-sun-ge6-audit full-q3-exceptional-audit full-q3-bcd-high-audit full-q3-bcd-mid-audit full-q3-bd-residual-audit full-q3-bcd-low-tail-audit full-q3-bcd-remaining-audit full-q3-bcd-bounded-littlewood-audit full-q3-extension

clean-room-replay:
	python3 clean_room_replay.py --threads $(REPLAY_THREADS)

full-q3-document-audit:
	python3 verify_full_q3_document.py

full-q3-artifact-audit:
	python3 -c 'from pathlib import Path; from clean_room_replay import verify_manifests; n, e = verify_manifests(Path(".")); print(f"FULL_Q3_ARTIFACT manifests={n} entries={e}"); print("FULL_Q3_ARTIFACT VERIFICATION: ALL PASS")'

full-q3-stable-audit:
	python3 verify_full_q3_stable.py

full-q3-su2-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_su2_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_su2_gmp
	/tmp/verify_full_q3_su2_gmp

full-q3-su3-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_su3_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_su3_gmp
	/tmp/verify_full_q3_su3_gmp

full-q3-su45-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_su45_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_su45_gmp
	/tmp/verify_full_q3_su45_gmp

full-q3-sun-ge6-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_sun_ge6_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_sun_ge6_gmp
	/tmp/verify_full_q3_sun_ge6_gmp

full-q3-exceptional-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_exceptional_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_exceptional_gmp
	/tmp/verify_full_q3_exceptional_gmp

full-q3-bcd-high-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/trace_square_cutoff_mpfr.cpp -l:libmpfr.so.6 -lgmp \
		-o /tmp/trace_square_cutoff_mpfr
	/tmp/trace_square_cutoff_mpfr
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_bcd_high_mpfr.cpp -l:libmpfr.so.6 -lgmp \
		-o /tmp/verify_full_q3_bcd_high_mpfr
	/tmp/verify_full_q3_bcd_high_mpfr

full-q3-bcd-mid-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_bcd_mid_mpfr.cpp -l:libmpfr.so.6 -lgmp \
		-o /tmp/verify_full_q3_bcd_mid_mpfr
	/tmp/verify_full_q3_bcd_mid_mpfr

full-q3-bd-residual-audit:
	g++ -O3 -DNDEBUG -std=c++20 -fopenmp \
		character_ring_iter/verify_full_q3_bd_residual_gmp.cpp \
		-lgmpxx -lgmp -o /tmp/verify_full_q3_bd_residual_gmp
	OMP_NUM_THREADS=$(REPLAY_THREADS) \
		/tmp/verify_full_q3_bd_residual_gmp

full-q3-bcd-low-tail-audit:
	g++ -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Werror \
		character_ring_iter/verify_full_q3_bcd_low_tail_mpfr.cpp \
		-l:libmpfr.so.6 -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_bcd_low_tail_mpfr
	/tmp/verify_full_q3_bcd_low_tail_mpfr

full-q3-bcd-remaining-audit:
	g++ -O3 -DNDEBUG -std=c++20 -fopenmp -Wall -Wextra -Wpedantic \
		-Wconversion -Wsign-conversion -Werror \
		character_ring_iter/verify_full_q3_bcd_remaining_gmp.cpp \
		-lgmpxx -lgmp -o /tmp/verify_full_q3_bcd_remaining_gmp
	OMP_NUM_THREADS=$(FULL_Q3_REMAINING_THREADS) \
		/tmp/verify_full_q3_bcd_remaining_gmp

full-q3-bcd-bounded-littlewood-audit:
	g++ -O3 -DNDEBUG -std=c++20 -fopenmp -Wall -Wextra -Wpedantic \
		-Wconversion -Wsign-conversion -Werror \
		character_ring_iter/verify_full_q3_bcd_bounded_littlewood_gmp.cpp \
		-lgmpxx -lgmp -o /tmp/verify_full_q3_bcd_bounded_littlewood_gmp
	OMP_NUM_THREADS=$(FULL_Q3_BOUNDED_THREADS) \
		/tmp/verify_full_q3_bcd_bounded_littlewood_gmp

full-q3-extension: full-q3-document-audit full-q3-artifact-audit full-q3-stable-audit full-q3-su2-audit full-q3-su3-audit full-q3-su45-audit full-q3-sun-ge6-audit full-q3-exceptional-audit full-q3-bcd-high-audit full-q3-bcd-mid-audit full-q3-bd-residual-audit full-q3-bcd-low-tail-audit full-q3-bcd-bounded-littlewood-audit
	pdflatex -interaction=nonstopmode -halt-on-error full_q3_extension.tex
	pdflatex -interaction=nonstopmode -halt-on-error full_q3_extension.tex
