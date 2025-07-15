
#include "igusa.h"



int thomae_correct_signs(slong* perm, slong* signs, acb_srcptr roots,
		acb_srcptr I, slong prec) {
  time_pair start; timestamp_mark(&start);
	// p in {0..719} s in {0..15}
	slong nb_candidates = 720 * 16;
	slong weights[4] = IGUSA_HALFWEIGHTS;
	slong correct_perm = -1;
	slong correct_signs = -1;
	int res = 1;
	int v = get_thomae_verbose();

	int removed_ps[720][16] = { 0 }; // all elements 0
	int removed_p[720] = { 0 }; // all elements 0



	#pragma omp parallel
	{ // this initializes the private	variables
		acb_ptr new_roots;
		acb_ptr ros;
		acb_ptr th4, th2;
		acb_ptr I_test;
		acb_mat_t tau;
		new_roots = _acb_vec_init(6);
		ros = _acb_vec_init(3);
		th4 = _acb_vec_init(16);
		th2 = _acb_vec_init(16);
		I_test = _acb_vec_init(4);
		acb_mat_init(tau, 2, 2);
		slong current_prec = thomae_startprec(prec);


		// put a barrier around the while guard
		#pragma omp barrier
		while((nb_candidates > 1) && (current_prec < prec)) {
			#pragma omp barrier
			#pragma omp single
			{
				if(v) flint_printf("(thomae_correct_signs) Trying precision %wd\n", current_prec);
				nb_candidates = 0;
			}
			#pragma omp barrier
			#pragma omp for schedule(static)
			for (slong p = 0; p < 720; ++p) {
				if (removed_p[p]) continue;
				removed_p[p] = 1;
				thomae_reorder(new_roots, roots, p);
				thomae_rosenhain(ros, new_roots, current_prec);
				thomae_theta4(th4, ros, current_prec);
				for (slong s = 0; s < 16; ++s) {
					if ( removed_ps[p][s] ) continue;
					removed_ps[p][s] = 1;
					thomae_theta2(th2, th4, ros, s, current_prec);
					if (thomae_discard(th2, current_prec)) continue;
					int tau_success = theta2_inverse(tau, th2, current_prec);
					/* We have to keep them in case of failure to compute. */
					if (!tau_success || thomae_keep_candidate(tau, I, current_prec)) {
						#pragma omp atomic write
						correct_perm = p;
						#pragma omp atomic write
						correct_signs = s;
						removed_p[p] = 0;
						removed_ps[p][s] = 0;
						#pragma omp atomic update
						++nb_candidates;
					}
					/* Else, th2 was discarded: do nothing */
				}
				/* End for(s) loop */
			}
			/* End for(p) loop */
			// there is an implicit barrier at the end of the for loop
			#pragma omp single
			{
				if(v) flint_printf("(thomae_correct_signs) Remaining candidates: %wd\n", nb_candidates);
			}
			current_prec *= THOMAE_MULPREC;
			// put a barrier around the while guard
			#pragma omp barrier
		}

		// put a barrier around the if guard
		#pragma omp barrier
		/* Then we have a last run at prec: this time we want to succeed in computing tau */
		if (nb_candidates > 1) {
		  #pragma omp barrier

		  #pragma omp single
			{
				if(v) flint_printf("(thomae_correct_signs) Last run at precision %wd\n", prec);
				nb_candidates = 0;
			}
		  #pragma omp barrier

		  #pragma omp for schedule(static)
			for (slong p = 0; p < 720; ++p) {
				if(nb_candidates) continue;
				if(removed_p[p]) continue;
				thomae_reorder(new_roots, roots, p);
				thomae_rosenhain(ros, new_roots, prec);
				thomae_theta4(th4, ros, prec);
				for (slong s = 0; s < 16; ++s) {
					if(nb_candidates) continue;
					if(removed_ps[p][s]) continue;
					// initialize
					thomae_theta2(th2, th4, ros, s, current_prec);
					if(thomae_discard(th2, current_prec)) continue;
					int tau_success = theta2_inverse(tau, th2, current_prec);
					if (!tau_success) continue;
					tau_success = igusa_from_tau(I_test, tau, prec);
					if (!tau_success) continue;
					tau_success = !cov_distinct(I_test, I, 4, weights, prec);
					if (!tau_success) continue;
					#pragma omp critical
					{
						if(!nb_candidates) {
							nb_candidates = 1;
							correct_perm = p;
							correct_signs = s;

							*perm = correct_perm;
							*signs = correct_signs;
							/* acb_mat_printd(tau, 10); */
							if(v) flint_printf("(thomae_correct_signs) Found one! p=%wd s=%wd\n", p, s);
						}
					}
				}
			}
		}
		// this clears private variables
		_acb_vec_clear(new_roots, 6);
		_acb_vec_clear(ros, 3);
		_acb_vec_clear(th4, 16);
		_acb_vec_clear(th2, 16);
		_acb_vec_clear(I_test, 4);
		acb_mat_clear(tau);
	}
	if (nb_candidates == 0) {
		res = 0;
	} else {
		*perm = correct_perm;
		*signs = correct_signs;
	}
	if (v) report_end(start);
	return res;
}
