
#include "modular.h"

int hilbert_modeq_eval_split(modeq_t R1, modeq_t R2, modeq_ctx_t ctx1,
			     modeq_ctx_t ctx2, fmpz* I, slong ell, slong delta)
{
  hecke_t H;
  modeq_acb_t E;
  acb_t c;
  fmpz_poly_t beta;
  slong nb = hilbert_nb_cosets(ell, delta);
  slong prec = hilbert_modeq_startprec(I, ell, delta);
  int res;
  int v = get_modeq_verbose();
  int stop = 0;
  
  hecke_init(H, nb);
  modeq_acb_init(E);
  fmpz_poly_init(beta);
  acb_init(c);

  res = hilbert_splits(beta, ell, delta);
  if (!res)
    {
      flint_printf("(hilbert_modeq_eval_split) Error: prime does not split\n");
      fflush(stdout);
      flint_abort();
    }     
  
  while (!stop)
    {
      if (v) modeq_verbose_start(prec);
      
      res = hecke_set_I_fmpz_hilbert(H, I, delta, prec);
      if (res) res = hecke_collect_hilbert(H, beta, ell, delta, prec);
      if (res) res = modeq_ctx_choose(ctx1, hecke_all_I(H), nb, prec);
      if (res) modeq_product_trees(E, H, ctx1, prec);
      if (res) res = modeq_rationalize(R1, E, prec);

      hilbert_conjugate(beta, beta, delta);
            
      if (res) res = hecke_collect_hilbert(H, beta, ell, delta, prec);
      if (res) res = modeq_ctx_choose(ctx2, hecke_all_I(H), nb, prec);
      if (res) modeq_product_trees(E, H, ctx2, prec);
      if (res) res = modeq_rationalize(R2, E, prec);

      prec = modeq_nextprec_generic(prec);
      stop = modeq_stop(res, prec);
    }
  
  hecke_clear(H);
  modeq_acb_clear(E);
  fmpz_poly_clear(beta);
  acb_clear(c);
  return res;
}
