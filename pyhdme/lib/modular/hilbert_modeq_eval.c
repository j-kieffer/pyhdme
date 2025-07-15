
#include "modular.h"

int hilbert_modeq_eval(modeq_t R, modeq_ctx_t ctx, fmpz* I,
		       slong ell, slong delta)
{  
  hecke_t H;
  modeq_acb_t E;
  acb_t c;
  slong nb = 2*hilbert_nb_cosets(ell, delta);
  slong prec = hilbert_modeq_startprec(I, ell, delta);
  slong gap;
  int res;
  int v = get_modeq_verbose();
  int stop = 0;
  
  hecke_init(H, nb);
  modeq_acb_init(E);
  acb_init(c);
  
  while (!stop)
    {
      if (v) modeq_verbose_start(prec);
      
      res = hecke_set_I_fmpz_hilbert(H, I, delta, prec);
      if (res) res = hecke_collect_hilbert_sym(H, ell, delta, prec);
      if (res) res = modeq_ctx_choose(ctx, hecke_all_I(H), nb, prec);
      if (res) modeq_product_trees(E, H, ctx, prec);
      if (res && (delta == 5))
	{
	  modeq_scalar(c, H, I, ctx, prec);
	  modeq_rescale(E, E, c, prec);
	  res = modeq_round(R, &gap, E);
	  prec = modeq_nextprec_precise(prec, gap);
	}
      else
	{
	  prec = modeq_nextprec_generic(prec);
	  if (res) res = modeq_rationalize(R, E, prec);
	}

      stop = modeq_stop(res, prec);
    }
  
  hecke_clear(H);
  modeq_acb_clear(E);
  acb_clear(c);
  return res;
}
