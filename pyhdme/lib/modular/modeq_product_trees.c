
#include "modular.h"

void modeq_product_trees(modeq_acb_t E, const hecke_t H,
			 const modeq_ctx_t ctx, slong prec)
{
  acb_ptr dens;
  acb_ptr nums;
  acb_ptr aux;
  acb_t x;
  slong d = hecke_nb(H);
  slong j, k;

  dens = _acb_vec_init(d);
  nums = _acb_vec_init(d);
  aux = _acb_vec_init(d);
  acb_init(x);

  modeq_nb(E) = modeq_ctx_nb(ctx);
  modeq_degree(E) = d;
  
  /* Collect values of parameter num/den at all isog matrices */
  for (k = 0; k < d; k++)
    {
      cov_mpoly_eval(&dens[k], modeq_ctx_den(ctx), hecke_I(H, k), modeq_ctx_ctx(ctx), prec);
      cov_mpoly_eval(&nums[k], modeq_ctx_num(ctx), hecke_I(H, k), modeq_ctx_ctx(ctx), prec);
      acb_neg(&nums[k], &nums[k]);
    }

  /* Compute denominator and equation */
  acb_one(x);
  for (k = 0; k < d; k++)
  {
      acb_mul(x, x, &dens[k], prec);
  }
  acb_set(modeq_den(E), x);
  acb_poly_product_tree_1(modeq_equation(E), dens, nums, d, prec);

  /* Construct equations for other monomials */
  for (j = 0; j < modeq_ctx_nb(ctx); j++)
    {
      for (k = 0; k < d; k++)
	{
	  cov_mpoly_eval(&aux[k], modeq_ctx_monomial(ctx, j), hecke_I(H, k),
			 modeq_ctx_ctx(ctx), prec);
	}
      acb_poly_product_tree_2(modeq_interpolate(E, j), dens, nums, aux, d, prec);
    }  

  _acb_vec_clear(dens, d);
  _acb_vec_clear(nums, d);
  _acb_vec_clear(aux, d);
  acb_clear(x);
}
