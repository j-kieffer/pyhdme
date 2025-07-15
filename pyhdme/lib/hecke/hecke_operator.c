
#include "hecke.h"

void hecke_operator(acb_ptr im, const hecke_t H, acb_srcptr val,
		   slong ell, slong k, slong j, slong prec)
{
  slong nb = hecke_nb(H);
  slong len = j+1;
  slong i;
  acb_ptr term;
  acb_ptr res;
  acb_t scal;

  term = _acb_vec_init(len);
  res = _acb_vec_init(len);
  acb_init(scal);

  if (k < 0 || j < 0)
    {
      flint_printf("(hecke_operator) Error: weights %wd, %wd must be nonnegative\n", k, j);
      fflush(stdout);
      flint_abort();
    }

  for (i = 0; i < nb; i++)
    {
      hecke_slash(term, hecke_star(H, i), &val[i*len], k, j, prec);
      _acb_vec_add(res, res, term, len, prec);
    }
  acb_set_si(scal, ell);
  acb_pow_si(scal, scal, 2*k + j - 3, prec);
  _acb_vec_scalar_mul(res, res, len, scal, prec);

  _acb_vec_set(im, res, len);

  _acb_vec_clear(term, len);
  _acb_vec_clear(res, len);
}
