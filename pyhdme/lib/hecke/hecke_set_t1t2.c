
#include "hecke.h"

int hecke_set_t1t2(hecke_t H, acb_srcptr t, slong delta, slong prec)
{
  int res;
  
  _acb_vec_set(hecke_t1t2(H), t, 2);
  hilbert_map(hecke_tau(H), hecke_t1t2(H), delta, prec);
  res = hecke_set_tau(H, hecke_tau(H), prec);

  return res;
}
