
#include "modular.h"

int siegel_modeq_eval(modeq_t R, modeq_ctx_t ctx, fmpz* I, slong ell)
{
  hecke_t H;
  modeq_acb_t E;
  acb_t c;
  slong nb = siegel_nb_cosets(ell);
  slong prec = siegel_modeq_startprec(I, ell);
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
      
      res = hecke_set_I_fmpz(H, I, prec);
      if (res) res = hecke_collect_siegel(H, ell, prec);
      if (res) res = modeq_ctx_choose(ctx, hecke_all_I(H), nb, prec);
      if (res)
	{
	  modeq_scalar(c, H, I, ctx, prec);
	  modeq_product_trees(E, H, ctx, prec);
	  modeq_rescale(E, E, c, prec);
	  res = modeq_round(R, &gap, E);
	  prec = modeq_nextprec_precise(prec, gap);
	}
      else
	{	  
	  prec = modeq_nextprec_generic(prec);
	}

      stop = modeq_stop(res, prec);
    }
  
  hecke_clear(H);
  modeq_acb_clear(E);
  acb_clear(c);
  return res;
}
