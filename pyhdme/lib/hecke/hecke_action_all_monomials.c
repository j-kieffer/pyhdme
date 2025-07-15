
#include "hecke.h"

void hecke_action_all_monomials(acb_ptr r, const hecke_t H, slong wt,
				slong ell, slong prec)
{
  slong weights[4] = IGUSA_WEIGHTS;
  slong m = cov_nb_monomials(wt, 4, weights);
  slong nb = hecke_nb(H);
  acb_mat_t monomials;
  acb_ptr row;
  acb_ptr col;
  slong k, j;

  acb_mat_init(monomials, nb, m);
  row = _acb_vec_init(m);
  col = _acb_vec_init(nb);
  
  for (j = 0; j < nb; j++)
    {
      cov_eval_all_monomials(row, hecke_I(H, j), wt, 4, weights, prec);
      for (k = 0; k < m; k++)
	{
	  acb_set(acb_mat_entry(monomials, j, k), &row[k]);
	}
    }

  for (k = 0; k < m; k++)
    {
      for (j = 0; j < nb; j++)
	{
	  acb_set(&col[j], acb_mat_entry(monomials, j, k));
	}
      hecke_operator(&r[k], H, col, ell, wt, 0, prec);
    }
  
  acb_mat_clear(monomials);
  _acb_vec_clear(row, m);
  _acb_vec_clear(col, nb);
}
