
#include "hecke.h"

int hecke_basis_matrix(acb_mat_t basis, slong nb, const acb_mat_struct* pts,
		       slong wt, slong prec)
{
  acb_ptr I;
  acb_ptr M;
  
  slong weights[4] = IGUSA_WEIGHTS;
  slong m = cov_nb_monomials(wt, 4, weights);
  slong j, k;
  int res = 1;

  I = _acb_vec_init(4);
  M = _acb_vec_init(m);

  for (k = 0; k < nb; k++)
    {
      res = igusa_from_tau(I, &pts[k], prec);
      if (!res) break;
      cov_eval_all_monomials(M, I, wt, 4, weights, prec);
      for (j = 0; j < m; j++)
	{
	  acb_set(acb_mat_entry(basis, j, k), &M[j]);
	}
    }

  _acb_vec_clear(I, 4);
  _acb_vec_clear(M, m);
  return res;
}
