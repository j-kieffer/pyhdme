
#include "igusa.h"

static void perm_nb_to_images(slong* images, slong nb, slong n)
{
  slong fact;
  slong k, q, r;
  
  if (n == 0) images[0] = 0;
  else if (n < 0)
    {
      flint_printf("(perm_nb_to_images) Negative length\n");
      flint_abort();
    }
  else
    {
      fact = 1;
      for (k = 2; k < n; k++) fact *= k; /* fact = (n-1)! */
      q = nb / fact;
      r = nb % fact;
      images[0] = q;
      perm_nb_to_images(&images[1], r, n-1);
      for (k = 1; k < n; k++)
	{
	  if (images[k] >= q) images[k] += 1;
	}
    }
}

void thomae_reorder(acb_ptr new_roots, acb_srcptr roots, slong perm)
{
  acb_ptr aux;
  slong images[6];
  slong k;
  
  perm_nb_to_images(images, perm, 6);
  aux = _acb_vec_init(6); /* Do this AFTER computing images, otherwise
			     segfault: what is this mess? Compiler optimization? */
  
  for (k = 0; k < 6; k++)
    {
      acb_set(&aux[images[k]], &roots[k]);
    }

  _acb_vec_set(new_roots, aux, 6);
  _acb_vec_clear(aux, 6);
}
