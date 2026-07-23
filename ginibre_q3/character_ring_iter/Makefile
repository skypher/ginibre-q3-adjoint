CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++20 -Wall
LIBS = -lgmpxx -lgmp

# These closed-form frontier checkers finish in hundredths of a second even
# without optimization.  Avoid spending the clean replay's critical path on
# optimizing 19 tiny, one-shot translation units.
REPLAY_LIGHT_TARGETS = post_m29_bc_ninth_frontier_gmp \
	post_m29_bc_tenth_frontier_gmp \
	post_m29_bc_eleventh_frontier_gmp \
	post_m29_bc_twelfth_frontier_gmp \
	post_m29_bc_thirteenth_frontier_gmp \
	post_m29_bc_fourteenth_frontier_gmp \
	post_m29_bc_fifteenth_frontier_gmp \
	post_m29_bc_sixteenth_frontier_gmp \
	post_m29_bc_seventeenth_frontier_gmp \
	post_m29_bc_eighteenth_frontier_gmp \
	post_m29_bc_nineteenth_frontier_gmp \
	post_m29_bc_twentieth_frontier_gmp \
	post_m29_bc_twentyfirst_frontier_gmp \
	post_m29_bc_twentysecond_frontier_gmp \
	post_m29_bc_twentythird_frontier_gmp \
	post_m29_bc_twentyfourth_frontier_gmp \
	post_m29_bc_twentyfifth_frontier_gmp \
	post_m29_bc_twentysixth_frontier_gmp \
	post_m29_bc_twentyseventh_frontier_gmp
$(REPLAY_LIGHT_TARGETS): CXXFLAGS = -O0 -march=native -std=c++20 -Wall

REPLAY_TARGETS = verify_sun_repair_certificates_gmp \
	verify_character_ring_moment_sources_gmp \
	verify_bc_row_gated_bridge_bounded_littlewood_gmp \
	dunkl_exceptional_coefficients_gmp \
	direct_chain_rect_g2_f4_certificate \
	direct_chain_rect_e6_e7_rank_certificate \
	trace_square_cutoff_mpfr \
	verify_e8_negative_bounds_gmp \
	direct_chain_rect_e8_parallel_certificate \
	post_m29_b_unitary_square_trace_mpfr \
	post_m29_bc_tilted_mgf_mpfr \
	post_m29_bc_layered_mgf_mpfr \
	post_m29_d_stable_window_interval_gmp \
	post_m29_d_interval_bridge_frontier_gmp \
	classical_boundary_source_replay_gmp \
	verify_d10_exact_tail_gmp \
	verify_active_bc_b_frontiers_bounded_littlewood_gmp \
	post_m29_bc_cheb80_negative_tail_gmp \
	post_m29_bc_interval_formula_mpfr \
	post_m29_bc_interval_bridge_frontier_gmp \
	post_m29_bc_pieri_powerloss_gmp \
	post_m29_bc_b14_j37_badshape_gmp \
	post_m29_bc_frontiers_unified


all: moments moments_omp verify_sun_repair_certificates_gmp post_m29_direct_tail_log_scan post_m29_d_interval_formula_scan post_m29_d_interval_formula_mpfr post_m29_bc_interval_formula_mpfr post_m29_b_unitary_square_trace_mpfr post_m29_bc_tilted_mgf_mpfr post_m29_bc_layered_mgf_mpfr post_m29_bc_cheb80_negative_tail_gmp post_m29_d_interval_bridge_frontier_scan post_m29_d_interval_bridge_frontier_gmp post_m29_bc_interval_bridge_frontier_gmp post_m29_bc_bridge_lambda_scan post_m29_bc_pieri_compression_scan post_m29_bc_pieri_compression_gmp post_m29_bc_pieri_fullwall_gmp post_m29_bc_pieri_compression_pow2_gmp post_m29_bc_pieri_powerloss_gmp post_m29_bc_badshape_residual_probe post_m29_bc_fixed_witness_probe post_m29_bc_fixed_witness_2q_obstruction_gmp post_m29_bc_complement_packing_gmp post_m29_bc_complement_frontier_gmp post_m29_bc_fixed_even_envelope_gmp post_m29_bc_fixed_even_envelope_m29_obstruction_gmp post_m29_bc_fixed_even_envelope_m30_frontier_gmp post_m29_bc_fixed_even_envelope_m31_frontier_gmp post_m29_bc_low_low_profile_gmp post_m29_bc_b14_j37_badshape_gmp post_m29_bc_b_eighth_frontier_gmp post_m29_bc_ninth_frontier_gmp post_m29_bc_tenth_frontier_gmp post_m29_bc_eleventh_frontier_gmp post_m29_bc_twelfth_frontier_gmp post_m29_bc_thirteenth_frontier_gmp post_m29_bc_fourteenth_frontier_gmp post_m29_bc_fifteenth_frontier_gmp post_m29_bc_sixteenth_frontier_gmp post_m29_bc_seventeenth_frontier_gmp post_m29_bc_eighteenth_frontier_gmp post_m29_bc_nineteenth_frontier_gmp post_m29_bc_twentieth_frontier_gmp post_m29_bc_twentyfirst_frontier_gmp post_m29_bc_twentysecond_frontier_gmp post_m29_bc_twentythird_frontier_gmp post_m29_bc_twentyfourth_frontier_gmp post_m29_bc_twentyfifth_frontier_gmp post_m29_bc_twentysixth_frontier_gmp post_m29_bc_twentyseventh_frontier_gmp post_m29_bc_twentyeighth_frontier_gmp post_m29_bc_twentyeighth_b_absorption_gmp post_m29_bc_twentyninth_c_frontier_gmp post_m29_bc_twentyninth_b_absorption_gmp post_m29_bc_twentyninth_b_frontier_gmp

all: verify_d_high_rank_finite_tails_gmp

all: dunkl_exceptional_coefficients_gmp

replay-build: $(REPLAY_TARGETS)

verify_su2_plucker_second_syzygies: verify_su2_plucker_second_syzygies.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -Werror -fopenmp $< -o $@

verify_su2_v2_presburger: verify_su2_v2_presburger.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -lz3 -o $@

verify_su2_opd_four_suffix_z3: verify_su2_opd_four_suffix_z3.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -Werror $< -lz3 -pthread -o $@

analyze_su2_opd_four_suffix: analyze_su2_opd_four_suffix.cpp
	$(CXX) -O3 -march=native -std=c++20 -Wall -Wextra -Wpedantic \
		-Wconversion -Wshadow -Werror $< -o $@

analyze_su2_bc2_hypergroup: analyze_su2_bc2_hypergroup.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -Werror $< -o $@

analyze_su2_plucker_sr_parity: analyze_su2_plucker_sr_parity.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -Werror $< -o $@

analyze_su2_plucker_kernel: analyze_su2_plucker_kernel.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Werror $< -o $@

analyze_su2_additive_compound_sign: analyze_su2_additive_compound_sign.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Werror $< -o $@

search_su2_fundamental_transfer: search_su2_fundamental_transfer.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -Werror $< -o $@

search_su2_level_partial_character_cone: search_su2_level_partial_character_cone.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -Werror $< -o $@

search_su2_odd_orbit_gks: search_su2_odd_orbit_gks.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

analyze_su2_odd_orbit_rank3_cone: analyze_su2_odd_orbit_rank3_cone.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

verify_su2_odd_orbit_rank4_transport: verify_su2_odd_orbit_rank4_transport.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

search_su2_odd_orbit_tropical: search_su2_odd_orbit_tropical.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

analyze_su2_odd_orbit_transport_fan: analyze_su2_odd_orbit_transport_fan.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

analyze_su2_tadpole_ray: analyze_su2_tadpole_ray.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

probe_separated_pole_turan: probe_separated_pole_turan.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

analyze_su2_gt_exchange: analyze_su2_gt_exchange.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -Werror $< -o $@

verify_su2_equal_terminal_exception: verify_su2_equal_terminal_exception.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

search_su2_equal_terminal_tropical_minor: search_su2_equal_terminal_tropical_minor.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -lz3 -o $@

analyze_su2_outer_eight: analyze_su2_outer_eight.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

analyze_su2_outer_ten: analyze_su2_outer_ten.cpp
	$(CXX) -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Wshadow -Werror $< -o $@

post_m29_bc_frontiers_unified: post_m29_bc_frontiers_unified.cpp $(addsuffix .cpp,$(REPLAY_LIGHT_TARGETS)) frontier_runtime_scope.hpp
	$(CXX) -O0 -march=native -std=c++20 -Wall -fopenmp $< $(LIBS) -o $@
	@for alias in $(REPLAY_LIGHT_TARGETS); do ln -sf post_m29_bc_frontiers_unified $$alias; done

moments: moments.cpp lie_data.h
	$(CXX) $(CXXFLAGS) moments.cpp $(LIBS) -o moments

moments_omp: moments.cpp lie_data.h
	$(CXX) $(CXXFLAGS) -fopenmp moments.cpp $(LIBS) -o moments_omp

verify_character_ring_moment_sources_gmp: verify_character_ring_moment_sources_gmp.cpp lie_data.h ankerl/unordered_dense.h
	$(CXX) $(CXXFLAGS) -fopenmp verify_character_ring_moment_sources_gmp.cpp $(LIBS) -o verify_character_ring_moment_sources_gmp

verify_sun_repair_certificates_gmp: verify_sun_repair_certificates_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp verify_sun_repair_certificates_gmp.cpp $(LIBS) -o verify_sun_repair_certificates_gmp

verify_bc_row_gated_bridge_gmp: verify_bc_row_gated_bridge_gmp.cpp ankerl/unordered_dense.h
	$(CXX) $(CXXFLAGS) -fopenmp verify_bc_row_gated_bridge_gmp.cpp $(LIBS) -o verify_bc_row_gated_bridge_gmp

verify_bc_row_gated_bridge_bounded_littlewood_gmp: verify_bc_row_gated_bridge_bounded_littlewood_gmp.cpp verify_full_q3_bcd_bounded_littlewood_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Werror verify_bc_row_gated_bridge_bounded_littlewood_gmp.cpp \
		$(LIBS) -o verify_bc_row_gated_bridge_bounded_littlewood_gmp

classical_boundary_source_replay_gmp: classical_boundary_certificate.cpp ankerl/unordered_dense.h
	$(CXX) $(CXXFLAGS) classical_boundary_certificate.cpp $(LIBS) -o classical_boundary_source_replay_gmp

verify_d10_exact_tail_gmp: verify_d10_exact_tail_gmp.cpp
	$(CXX) $(CXXFLAGS) verify_d10_exact_tail_gmp.cpp $(LIBS) -o verify_d10_exact_tail_gmp

verify_full_q3_bcd_bounded_littlewood_gmp: verify_full_q3_bcd_bounded_littlewood_gmp.cpp full_q3_bcd_remaining_data.hpp
	$(CXX) $(CXXFLAGS) -fopenmp -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Werror verify_full_q3_bcd_bounded_littlewood_gmp.cpp \
		$(LIBS) -o verify_full_q3_bcd_bounded_littlewood_gmp

verify_active_bc_b_frontiers_bounded_littlewood_gmp: verify_active_bc_b_frontiers_bounded_littlewood_gmp.cpp verify_full_q3_bcd_bounded_littlewood_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp -Wall -Wextra -Wpedantic -Wconversion \
		-Wsign-conversion -Werror verify_active_bc_b_frontiers_bounded_littlewood_gmp.cpp \
		$(LIBS) -o verify_active_bc_b_frontiers_bounded_littlewood_gmp

verify_d_high_rank_finite_tails_gmp: verify_d_high_rank_finite_tails_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp verify_d_high_rank_finite_tails_gmp.cpp $(LIBS) -o verify_d_high_rank_finite_tails_gmp

dunkl_exceptional_coefficients_gmp: dunkl_exceptional_coefficients_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp dunkl_exceptional_coefficients_gmp.cpp $(LIBS) -o dunkl_exceptional_coefficients_gmp

direct_chain_rect_g2_f4_certificate: direct_chain_rect_g2_f4_certificate.cpp
	$(CXX) $(CXXFLAGS) direct_chain_rect_g2_f4_certificate.cpp $(LIBS) -o direct_chain_rect_g2_f4_certificate

direct_chain_rect_e6_e7_rank_certificate: direct_chain_rect_e6_e7_rank_certificate.cpp
	$(CXX) $(CXXFLAGS) -fopenmp direct_chain_rect_e6_e7_rank_certificate.cpp $(LIBS) -o direct_chain_rect_e6_e7_rank_certificate

trace_square_cutoff_mpfr: trace_square_cutoff_mpfr.cpp
	$(CXX) $(CXXFLAGS) trace_square_cutoff_mpfr.cpp -l:libmpfr.so.6 -lgmp -o trace_square_cutoff_mpfr

verify_e8_negative_bounds_gmp: verify_e8_negative_bounds_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp verify_e8_negative_bounds_gmp.cpp $(LIBS) -o verify_e8_negative_bounds_gmp

direct_chain_rect_e8_parallel_certificate: direct_chain_rect_e8_parallel_certificate.cpp
	$(CXX) $(CXXFLAGS) -fopenmp direct_chain_rect_e8_parallel_certificate.cpp -l:libmpfr.so.6 $(LIBS) -o direct_chain_rect_e8_parallel_certificate

post_m29_direct_tail_log_scan: post_m29_direct_tail_log_scan.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_direct_tail_log_scan.cpp -o post_m29_direct_tail_log_scan

post_m29_d_interval_formula_scan: post_m29_d_interval_formula_scan.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_d_interval_formula_scan.cpp -o post_m29_d_interval_formula_scan

post_m29_d_interval_formula_mpfr: post_m29_d_interval_formula_mpfr.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_d_interval_formula_mpfr.cpp -l:libmpfr.so.6 -lgmp -o post_m29_d_interval_formula_mpfr

post_m29_d_two_sided_obstruction_gmp: post_m29_d_two_sided_obstruction_gmp.cpp
	$(CXX) $(CXXFLAGS) post_m29_d_two_sided_obstruction_gmp.cpp $(LIBS) -o post_m29_d_two_sided_obstruction_gmp

post_m29_bc_interval_formula_mpfr: post_m29_bc_interval_formula_mpfr.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_interval_formula_mpfr.cpp -l:libmpfr.so.6 -lgmp -o post_m29_bc_interval_formula_mpfr

post_m29_b_unitary_square_trace_mpfr: post_m29_b_unitary_square_trace_mpfr.cpp
	$(CXX) $(CXXFLAGS) post_m29_b_unitary_square_trace_mpfr.cpp -l:libmpfr.so.6 -lgmp -o post_m29_b_unitary_square_trace_mpfr

post_m29_bc_tilted_mgf_mpfr: post_m29_bc_tilted_mgf_mpfr.cpp
	$(CXX) $(CXXFLAGS) post_m29_bc_tilted_mgf_mpfr.cpp -l:libmpfr.so.6 -lgmp -o post_m29_bc_tilted_mgf_mpfr

post_m29_bc_layered_mgf_mpfr: post_m29_bc_layered_mgf_mpfr.cpp
	$(CXX) $(CXXFLAGS) post_m29_bc_layered_mgf_mpfr.cpp -l:libmpfr.so.6 $(LIBS) -o post_m29_bc_layered_mgf_mpfr

post_m29_d_stable_window_interval_gmp: post_m29_d_stable_window_interval_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_d_stable_window_interval_gmp.cpp $(LIBS) -o post_m29_d_stable_window_interval_gmp

post_m29_bc_cheb80_negative_tail_gmp: post_m29_bc_cheb80_negative_tail_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_cheb80_negative_tail_gmp.cpp $(LIBS) -o post_m29_bc_cheb80_negative_tail_gmp

post_m29_d_interval_bridge_frontier_scan: post_m29_d_interval_bridge_frontier_scan.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_d_interval_bridge_frontier_scan.cpp -o post_m29_d_interval_bridge_frontier_scan

post_m29_d_interval_bridge_frontier_gmp: post_m29_d_interval_bridge_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_d_interval_bridge_frontier_gmp.cpp $(LIBS) -o post_m29_d_interval_bridge_frontier_gmp

post_m29_bc_interval_bridge_frontier_gmp: post_m29_bc_interval_bridge_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_interval_bridge_frontier_gmp.cpp $(LIBS) -o post_m29_bc_interval_bridge_frontier_gmp

post_m29_bc_bridge_lambda_scan: post_m29_bc_bridge_lambda_scan.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_bridge_lambda_scan.cpp -o post_m29_bc_bridge_lambda_scan

post_m29_bc_pieri_compression_scan: post_m29_bc_pieri_compression_scan.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_pieri_compression_scan.cpp -o post_m29_bc_pieri_compression_scan

post_m29_bc_pieri_compression_gmp: post_m29_bc_pieri_compression_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_pieri_compression_gmp.cpp $(LIBS) -o post_m29_bc_pieri_compression_gmp

post_m29_bc_pieri_fullwall_gmp: post_m29_bc_pieri_fullwall_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_pieri_fullwall_gmp.cpp $(LIBS) -o post_m29_bc_pieri_fullwall_gmp

post_m29_bc_pieri_compression_pow2_gmp: post_m29_bc_pieri_compression_pow2_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_pieri_compression_pow2_gmp.cpp $(LIBS) -o post_m29_bc_pieri_compression_pow2_gmp

post_m29_bc_pieri_powerloss_gmp: post_m29_bc_pieri_powerloss_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_pieri_powerloss_gmp.cpp $(LIBS) -o post_m29_bc_pieri_powerloss_gmp

post_m29_bc_badshape_residual_probe: post_m29_bc_badshape_residual_probe.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_badshape_residual_probe.cpp -o post_m29_bc_badshape_residual_probe

post_m29_bc_fixed_witness_probe: post_m29_bc_fixed_witness_probe.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fixed_witness_probe.cpp -o post_m29_bc_fixed_witness_probe

post_m29_bc_fixed_witness_2q_obstruction_gmp: post_m29_bc_fixed_witness_2q_obstruction_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fixed_witness_2q_obstruction_gmp.cpp $(LIBS) -o post_m29_bc_fixed_witness_2q_obstruction_gmp

post_m29_bc_complement_packing_gmp: post_m29_bc_complement_packing_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_complement_packing_gmp.cpp $(LIBS) -o post_m29_bc_complement_packing_gmp

post_m29_bc_complement_frontier_gmp: post_m29_bc_complement_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_complement_frontier_gmp.cpp $(LIBS) -o post_m29_bc_complement_frontier_gmp

post_m29_bc_fixed_even_envelope_gmp: post_m29_bc_fixed_even_envelope_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fixed_even_envelope_gmp.cpp $(LIBS) -o post_m29_bc_fixed_even_envelope_gmp

post_m29_bc_fixed_even_envelope_m29_obstruction_gmp: post_m29_bc_fixed_even_envelope_m29_obstruction_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fixed_even_envelope_m29_obstruction_gmp.cpp $(LIBS) -o post_m29_bc_fixed_even_envelope_m29_obstruction_gmp

post_m29_bc_fixed_even_envelope_m30_frontier_gmp: post_m29_bc_fixed_even_envelope_m30_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fixed_even_envelope_m30_frontier_gmp.cpp $(LIBS) -o post_m29_bc_fixed_even_envelope_m30_frontier_gmp

post_m29_bc_fixed_even_envelope_m31_frontier_gmp: post_m29_bc_fixed_even_envelope_m31_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fixed_even_envelope_m31_frontier_gmp.cpp $(LIBS) -o post_m29_bc_fixed_even_envelope_m31_frontier_gmp

post_m29_bc_low_low_profile_gmp: post_m29_bc_low_low_profile_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_low_low_profile_gmp.cpp $(LIBS) -o post_m29_bc_low_low_profile_gmp

post_m29_bc_b14_j37_badshape_gmp: post_m29_bc_b14_j37_badshape_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_b14_j37_badshape_gmp.cpp $(LIBS) -o post_m29_bc_b14_j37_badshape_gmp

post_m29_bc_b_eighth_frontier_gmp: post_m29_bc_b_eighth_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_b_eighth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_b_eighth_frontier_gmp

post_m29_bc_ninth_frontier_gmp: post_m29_bc_ninth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_ninth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_ninth_frontier_gmp

post_m29_bc_tenth_frontier_gmp: post_m29_bc_tenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_tenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_tenth_frontier_gmp

post_m29_bc_eleventh_frontier_gmp: post_m29_bc_eleventh_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_eleventh_frontier_gmp.cpp $(LIBS) -o post_m29_bc_eleventh_frontier_gmp

post_m29_bc_twelfth_frontier_gmp: post_m29_bc_twelfth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twelfth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twelfth_frontier_gmp

post_m29_bc_thirteenth_frontier_gmp: post_m29_bc_thirteenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_thirteenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_thirteenth_frontier_gmp

post_m29_bc_fourteenth_frontier_gmp: post_m29_bc_fourteenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fourteenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_fourteenth_frontier_gmp

post_m29_bc_fifteenth_frontier_gmp: post_m29_bc_fifteenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_fifteenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_fifteenth_frontier_gmp

post_m29_bc_sixteenth_frontier_gmp: post_m29_bc_sixteenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_sixteenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_sixteenth_frontier_gmp

post_m29_bc_seventeenth_frontier_gmp: post_m29_bc_seventeenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_seventeenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_seventeenth_frontier_gmp

post_m29_bc_eighteenth_frontier_gmp: post_m29_bc_eighteenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_eighteenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_eighteenth_frontier_gmp

post_m29_bc_nineteenth_frontier_gmp: post_m29_bc_nineteenth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_nineteenth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_nineteenth_frontier_gmp

post_m29_bc_twentieth_frontier_gmp: post_m29_bc_twentieth_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentieth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentieth_frontier_gmp

post_m29_bc_twentyfirst_frontier_gmp: post_m29_bc_twentyfirst_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentyfirst_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentyfirst_frontier_gmp

post_m29_bc_twentysecond_frontier_gmp: post_m29_bc_twentysecond_frontier_gmp.cpp frontier_runtime_scope.hpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentysecond_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentysecond_frontier_gmp

post_m29_bc_twentythird_frontier_gmp: post_m29_bc_twentythird_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentythird_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentythird_frontier_gmp

post_m29_bc_twentyfourth_frontier_gmp: post_m29_bc_twentyfourth_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentyfourth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentyfourth_frontier_gmp

post_m29_bc_twentyfifth_frontier_gmp: post_m29_bc_twentyfifth_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentyfifth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentyfifth_frontier_gmp

post_m29_bc_twentysixth_frontier_gmp: post_m29_bc_twentysixth_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentysixth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentysixth_frontier_gmp

post_m29_bc_twentyseventh_frontier_gmp: post_m29_bc_twentyseventh_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentyseventh_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentyseventh_frontier_gmp

post_m29_bc_twentyeighth_frontier_gmp: post_m29_bc_twentyeighth_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentyeighth_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentyeighth_frontier_gmp

post_m29_bc_twentyeighth_b_absorption_gmp: post_m29_bc_twentyeighth_b_absorption_gmp.cpp
	$(CXX) $(CXXFLAGS) post_m29_bc_twentyeighth_b_absorption_gmp.cpp $(LIBS) -o post_m29_bc_twentyeighth_b_absorption_gmp

post_m29_bc_twentyninth_c_frontier_gmp: post_m29_bc_twentyninth_c_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) post_m29_bc_twentyninth_c_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentyninth_c_frontier_gmp

post_m29_bc_twentyninth_b_absorption_gmp: post_m29_bc_twentyninth_b_absorption_gmp.cpp
	$(CXX) $(CXXFLAGS) post_m29_bc_twentyninth_b_absorption_gmp.cpp $(LIBS) -o post_m29_bc_twentyninth_b_absorption_gmp

post_m29_bc_twentyninth_b_frontier_gmp: post_m29_bc_twentyninth_b_frontier_gmp.cpp
	$(CXX) $(CXXFLAGS) -fopenmp post_m29_bc_twentyninth_b_frontier_gmp.cpp $(LIBS) -o post_m29_bc_twentyninth_b_frontier_gmp

lie_data.h: lie_data.py gen_header.py
	python3 gen_header.py

clean:
	rm -f verify_d_high_rank_finite_tails_gmp verify_d10_exact_tail_gmp verify_d10_d11_exact_tail_gmp verify_full_q3_bcd_bounded_littlewood_gmp dunkl_exceptional_coefficients_gmp direct_chain_rect_g2_f4_certificate direct_chain_rect_e6_e7_rank_certificate trace_square_cutoff_mpfr verify_e8_negative_bounds_gmp direct_chain_rect_e8_parallel_certificate post_m29_d_two_sided_obstruction_gmp classical_boundary_source_replay_gmp post_m29_bc_frontiers_unified
	rm -f search_su2_odd_orbit_gks analyze_su2_odd_orbit_rank3_cone verify_su2_odd_orbit_rank4_transport search_su2_odd_orbit_tropical analyze_su2_odd_orbit_transport_fan analyze_su2_tadpole_ray probe_separated_pole_turan
	rm -f moments moments_omp verify_sun_repair_certificates_gmp post_m29_direct_tail_log_scan post_m29_d_interval_formula_scan post_m29_d_interval_formula_mpfr post_m29_bc_interval_formula_mpfr post_m29_b_unitary_square_trace_mpfr post_m29_bc_tilted_mgf_mpfr post_m29_bc_layered_mgf_mpfr post_m29_d_stable_window_interval_gmp post_m29_bc_cheb80_negative_tail_gmp post_m29_d_interval_bridge_frontier_scan post_m29_d_interval_bridge_frontier_gmp post_m29_bc_interval_bridge_frontier_gmp post_m29_bc_bridge_lambda_scan post_m29_bc_pieri_compression_scan post_m29_bc_pieri_compression_gmp post_m29_bc_pieri_fullwall_gmp post_m29_bc_pieri_compression_pow2_gmp post_m29_bc_pieri_powerloss_gmp post_m29_bc_badshape_residual_probe post_m29_bc_fixed_witness_probe post_m29_bc_fixed_witness_2q_obstruction_gmp post_m29_bc_complement_packing_gmp post_m29_bc_complement_frontier_gmp post_m29_bc_fixed_even_envelope_gmp post_m29_bc_fixed_even_envelope_m29_obstruction_gmp post_m29_bc_fixed_even_envelope_m30_frontier_gmp post_m29_bc_fixed_even_envelope_m31_frontier_gmp post_m29_bc_low_low_profile_gmp post_m29_bc_b14_j37_badshape_gmp post_m29_bc_b_eighth_frontier_gmp post_m29_bc_ninth_frontier_gmp post_m29_bc_tenth_frontier_gmp post_m29_bc_eleventh_frontier_gmp post_m29_bc_twelfth_frontier_gmp post_m29_bc_thirteenth_frontier_gmp post_m29_bc_fourteenth_frontier_gmp post_m29_bc_fifteenth_frontier_gmp post_m29_bc_sixteenth_frontier_gmp post_m29_bc_seventeenth_frontier_gmp post_m29_bc_eighteenth_frontier_gmp post_m29_bc_nineteenth_frontier_gmp post_m29_bc_twentieth_frontier_gmp post_m29_bc_twentyfirst_frontier_gmp post_m29_bc_twentysecond_frontier_gmp post_m29_bc_twentythird_frontier_gmp post_m29_bc_twentyfourth_frontier_gmp post_m29_bc_twentyfifth_frontier_gmp post_m29_bc_twentysixth_frontier_gmp post_m29_bc_twentyseventh_frontier_gmp post_m29_bc_twentyeighth_frontier_gmp post_m29_bc_twentyeighth_b_absorption_gmp post_m29_bc_twentyninth_c_frontier_gmp post_m29_bc_twentyninth_b_absorption_gmp post_m29_bc_twentyninth_b_frontier_gmp *.log

.PHONY: all clean replay-build
