REPLAY_THREADS ?= $(shell nproc)
REPLAY_MAX_STAGE_SECONDS ?= 240
REPLAY_MAX_TOTAL_SECONDS ?= 240
PUBLICATION_SOURCE_DATE_EPOCH ?= 1784419200
export SOURCE_DATE_EPOCH := $(PUBLICATION_SOURCE_DATE_EPOCH)
# Each reverse-Pieri task carries a large GMP frontier.  The residual replay
# exposes at most four independent tasks, so using every logical CPU only
# multiplies peak memory without increasing useful parallelism.  Keep the
# reproducible default memory-safe; an audited host may override it explicitly.
FULL_Q3_REMAINING_THREADS ?= 2
FULL_Q3_BOUNDED_THREADS ?= $(shell nproc)
FULL_Q3_MODULAR_THREADS ?= $(shell nproc)

.PHONY: publication-preflight clean-room-replay paper-pdf full-q3-main-pdf full-q3-pdf unified-submission publication-artifact-audit final-publication-replay final-release-archive full-q3-document-audit full-q3-artifact-audit full-q3-stable-audit full-q3-su2-audit full-q3-su3-audit full-q3-su45-audit full-q3-sun-ge6-audit full-q3-exceptional-audit full-q3-bcd-high-audit full-q3-bcd-mid-audit full-q3-bd-residual-audit full-q3-bcd-low-tail-audit full-q3-bcd-remaining-audit full-q3-bcd-independent-audit full-q3-bcd-bounded-littlewood-audit full-q3-extension-independent-controls full-q3-extension

publication-preflight:
	python3 -u run_publication_short_audits.py
	python3 -u -c 'from pathlib import Path; from clean_room_replay import verify_manifests; n, e = verify_manifests(Path(".")); print(f"PUBLICATION_PREFLIGHT manifests={n} entries={e}"); print("PUBLICATION_PREFLIGHT VERIFICATION: ALL PASS")'

clean-room-replay:
	python3 -u clean_room_replay.py --threads $(REPLAY_THREADS) \
		--max-stage-seconds $(REPLAY_MAX_STAGE_SECONDS) \
		--max-total-seconds $(REPLAY_MAX_TOTAL_SECONDS)

paper-pdf:
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error paper.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error paper.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error paper.tex
	! grep -Eq 'undefined references|undefined citations|There were undefined references' paper.log

full-q3-pdf:
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_extension.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_extension.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_extension.tex
	! grep -Eq 'undefined references|undefined citations|There were undefined references' full_q3_extension.log

full-q3-main-pdf:
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_main.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_main.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_main.tex
	! grep -Eq 'undefined references|undefined citations|There were undefined references' full_q3_main.log

unified-submission: paper-pdf full-q3-main-pdf
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error submission.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error submission.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error submission.tex
	! grep -Eq 'undefined references|undefined citations|There were undefined references' submission.log

publication-artifact-audit: unified-submission full-q3-pdf
	python3 -u verify_publication_artifacts.py

final-publication-replay:
	./run_final_publication_replay.sh

final-release-archive:
	python3 -u build_final_release_archive.py \
		--replay-log "$(FINAL_REPLAY_LOG)" \
		--output "$(FINAL_RELEASE_ARCHIVE)"

full-q3-document-audit:
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s python3 verify_full_q3_document.py

full-q3-artifact-audit:
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s python3 -c 'from pathlib import Path; from clean_room_replay import verify_manifests; n, e = verify_manifests(Path(".")); print(f"FULL_Q3_ARTIFACT manifests={n} entries={e}"); print("FULL_Q3_ARTIFACT VERIFICATION: ALL PASS")'

full-q3-stable-audit:
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s python3 verify_full_q3_stable.py

full-q3-su2-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_su2_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_su2_gmp
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_su2_gmp

full-q3-su3-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_su3_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_su3_gmp
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_su3_gmp

full-q3-su45-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_su45_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_su45_gmp
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_su45_gmp

full-q3-sun-ge6-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_sun_ge6_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_sun_ge6_gmp
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_sun_ge6_gmp

full-q3-exceptional-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_exceptional_gmp.cpp -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_exceptional_gmp
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_exceptional_gmp

full-q3-bcd-high-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/trace_square_cutoff_mpfr.cpp -l:libmpfr.so.6 -lgmp \
		-o /tmp/trace_square_cutoff_mpfr
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/trace_square_cutoff_mpfr
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_bcd_high_mpfr.cpp -l:libmpfr.so.6 -lgmp \
		-o /tmp/verify_full_q3_bcd_high_mpfr
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_bcd_high_mpfr

full-q3-bcd-mid-audit:
	g++ -O2 -std=c++17 -Wall -Wextra -Werror \
		character_ring_iter/verify_full_q3_bcd_mid_mpfr.cpp -l:libmpfr.so.6 -lgmp \
		-o /tmp/verify_full_q3_bcd_mid_mpfr
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_bcd_mid_mpfr

full-q3-bd-residual-audit:
	g++ -O3 -DNDEBUG -std=c++20 -fopenmp \
		character_ring_iter/verify_full_q3_bd_residual_gmp.cpp \
		-lgmpxx -lgmp -o /tmp/verify_full_q3_bd_residual_gmp
	OMP_NUM_THREADS=$(REPLAY_THREADS) \
		timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s \
		/tmp/verify_full_q3_bd_residual_gmp --progress

full-q3-bcd-low-tail-audit:
	g++ -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Werror \
		character_ring_iter/verify_full_q3_bcd_low_tail_mpfr.cpp \
		-l:libmpfr.so.6 -lgmpxx -lgmp \
		-o /tmp/verify_full_q3_bcd_low_tail_mpfr
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s /tmp/verify_full_q3_bcd_low_tail_mpfr

full-q3-bcd-remaining-audit:
	g++ -O3 -DNDEBUG -std=c++20 -fopenmp -Wall -Wextra -Wpedantic \
		-Wconversion -Wsign-conversion -Werror \
		character_ring_iter/verify_full_q3_bcd_remaining_gmp.cpp \
		-lgmpxx -lgmp -o /tmp/verify_full_q3_bcd_remaining_gmp
	OMP_NUM_THREADS=$(FULL_Q3_REMAINING_THREADS) \
		timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s \
		/tmp/verify_full_q3_bcd_remaining_gmp --progress

full-q3-bcd-independent-audit:
	g++ -O3 -DNDEBUG -std=c++20 -fopenmp -Wall -Wextra -Wpedantic \
		-Wconversion -Wsign-conversion -Werror \
		character_ring_iter/verify_full_q3_bcd_modular_moment_checker.cpp \
		-lgmpxx -lgmp -o /tmp/verify_full_q3_bcd_modular_moment_checker
	OMP_NUM_THREADS=$(FULL_Q3_MODULAR_THREADS) \
		timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s \
		/tmp/verify_full_q3_bcd_modular_moment_checker \
		--moments-log certificates/full_q3/fullq3bcd0010_machine_c.log \
		--reverse-log certificates/full_q3/fullq3bcd0006_machine_c_prefix40.log \
		--progress

full-q3-bcd-bounded-littlewood-audit:
	g++ -O3 -DNDEBUG -std=c++20 -fopenmp -Wall -Wextra -Wpedantic \
		-Wconversion -Wsign-conversion -Werror \
		character_ring_iter/verify_full_q3_bcd_bounded_littlewood_gmp.cpp \
		-lgmpxx -lgmp -o /tmp/verify_full_q3_bcd_bounded_littlewood_gmp
	OMP_NUM_THREADS=$(FULL_Q3_BOUNDED_THREADS) \
		timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s \
		/tmp/verify_full_q3_bcd_bounded_littlewood_gmp --progress

full-q3-extension-independent-controls: full-q3-bd-residual-audit full-q3-bcd-independent-audit

full-q3-extension: full-q3-document-audit full-q3-artifact-audit full-q3-stable-audit full-q3-su2-audit full-q3-su3-audit full-q3-su45-audit full-q3-sun-ge6-audit full-q3-exceptional-audit full-q3-bcd-high-audit full-q3-bcd-mid-audit full-q3-bcd-low-tail-audit full-q3-bcd-bounded-littlewood-audit full-q3-bcd-independent-audit
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_extension.tex
	timeout --foreground $(REPLAY_MAX_STAGE_SECONDS)s pdflatex -interaction=nonstopmode -halt-on-error full_q3_extension.tex
