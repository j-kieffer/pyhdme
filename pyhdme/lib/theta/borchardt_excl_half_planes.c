
#include "theta.h"

/* Half planes are defined by an angle: an angle t defines the open
   half plane containing all complex numbers with arguments in the
   open interval ]t-pi, t[ of S^1. 
   
   Subintervals of S1 are of the form: ]-pi, b0[ \cup ]b1, b2[ \cup
   ]b3, pi]. */

void borchardt_excl_half_planes(arf_struct* b, const acb_t z, slong prec)
{
  arb_t pi;
  acb_t temp;
  arb_t arg;
  arf_t gap;
  slong k;

  arb_init(pi);
  acb_init(temp);
  arb_init(arg);
  arf_init(gap);
  
  arb_const_pi(pi, prec);
  arf_set_si(gap, 1);
  /* arf_mul_2exp_si(gap, gap, BORCHARDT_ARG_GAP_EXP); */
  /* Lengthen the excluded interval by gap? NO. */
  arf_zero(gap);
  
  /* Init to zero */
  for (k = 0; k < 4; k++) arf_zero(&b[k]);
  
  if (acb_contains_zero(z))
    {
      /* Empty interval */
      arf_set_si(&b[0], -1);
      arf_set_si(&b[3], 1);
    }
  else if (arb_contains_zero(acb_imagref(z))
	   && arb_is_negative(acb_realref(z)))
    {
      /* z intersects the negative real axis: multiply by -1 to get
	 argument, and -1 is not excluded */
      acb_neg(temp, z);
      acb_arg(arg, temp, prec);
      arb_div(arg, arg, pi, prec);
      arb_get_ubound_arf(&b[1], arg, prec); /* \geq 0 */
      arf_sub(&b[1], &b[1], gap, prec, ARF_RND_CEIL);
      
      arb_get_lbound_arf(&b[2], arg, prec); /* \leq 0 */
      arf_add_si(&b[2], &b[2], 1, prec, ARF_RND_FLOOR);
      arf_add(&b[2], &b[2], gap, prec, ARF_RND_FLOOR);
      
      arf_set_si(&b[0], -1);
      arf_set_si(&b[3], 1);
    }
  else if (arb_is_negative(acb_imagref(z)))
    {
      /* -1 is excluded */
      acb_arg(arg, z, prec);
      arb_div(arg, arg, pi, prec);
      arb_get_lbound_arf(&b[0], arg, prec);
      arf_add(&b[0], &b[0], gap, prec, ARF_RND_FLOOR);
      
      arb_get_ubound_arf(&b[3], arg, prec);
      arf_add_si(&b[3], &b[3], 1, prec, ARF_RND_CEIL);
      arf_sub(&b[3], &b[3], gap, prec, ARF_RND_CEIL);
    }
  else /* -1 is not excluded */
    {
      acb_arg(arg, z, prec);
      arb_div(arg, arg, pi, prec);
      arb_get_ubound_arf(&b[1], arg, prec);
      arf_sub_si(&b[1], &b[1], 1, prec, ARF_RND_CEIL);
      arf_sub(&b[1], &b[1], gap, prec, ARF_RND_CEIL);
      
      arb_get_lbound_arf(&b[2], arg, prec);
      arf_add(&b[2], &b[2], gap, prec, ARF_RND_FLOOR);
      
      arf_set_si(&b[0], -1);
      arf_set_si(&b[3], 1);
    }

  arb_clear(pi);
  acb_clear(temp);
  arb_clear(arg);
  arf_clear(gap);
}
