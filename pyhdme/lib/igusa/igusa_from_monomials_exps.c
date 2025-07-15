
#include "igusa.h"

/* See igusa_base_exps */

/* Weight 20 */
static slong i_0020 = 0;
static slong i_1110 = 1;
static slong i_2001 = 2;
static slong i_5000 = 3;
static slong i_2200 = 4;

/* Weight 30 */
static slong i_0500 = 0;
static slong i_3300 = 1;
static slong i_0030 = 2;
static slong i_5010 = 3;
static slong i_0301 = 4;
static slong i_1120 = 5;
static slong i_2011 = 6;

/* Weight 60 */
static slong i_0060 = 0;
static slong i_0005 = 1;
static slong i_0530 = 2;
static slong i_0204 = 3;
static slong i_3004 = 4;
static slong i_0X00 = 5;
static slong i_3800 = 6;
static slong i_5040 = 7;
static slong i_1150 = 8;
static slong i_3203 = 9;

#define SET_VEC_30(e, x0, x1, x2, x3, x4, x5, x6) {	\
    e[i_0500] = (x0);					\
    e[i_3300] = (x1);					\
    e[i_0030] = (x2);					\
    e[i_5010] = (x3);					\
    e[i_0301] = (x4);					\
    e[i_1120] = (x5);					\
    e[i_2011] = (x6);					\
  }

#define SET_VEC_60(e, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9) { \
    e[i_0060] = (x0);						\
    e[i_0005] = (x1);						\
    e[i_0530] = (x2);						\
    e[i_0204] = (x3);						\
    e[i_3004] = (x4);						\
    e[i_0X00] = (x5);						\
    e[i_3800] = (x6);						\
    e[i_5040] = (x7);						\
    e[i_1150] = (x8);						\
    e[i_3203] = (x9);						\
  }


void igusa_from_monomials_exps(slong* e4, slong* e6, slong* e10, slong* e12,
			       int z4, int z6, int z10, int z12, slong wt)
{
  slong nb = igusa_nb_base_monomials(wt);
  slong k;

  for (k = 0; k < nb; k++)
    {
      e4[k] = 0;
      e6[k] = 0;
      e10[k] = 0;
      e12[k] = 0;
    }
  
  if (z10 && z12)
    {	      
      flint_printf("(igusa_from_monomials_exps) Error: chi10, chi12 cannot be simultaneously zero\n");
      fflush(stdout);
      flint_abort();
    }
 
  if (wt == 20)
    {
      if (z4)
	{
	  flint_printf("(igusa_from_monomials_exps) Error: psi4 cannot be zero in case of weight 20\n");
	  fflush(stdout);
	  flint_abort();	  
	}
      if (!z10)
	{
	  e4[i_0020] = 1;
	  e4[i_5000] = 1;

	  e6[i_0020] = 1;
	  e6[i_1110] = 1;
	  e6[i_5000] = 1;

	  e10[i_5000] = 2;
	  e10[i_0020] = 3;

	  e12[i_5000] = 2;
	  e12[i_2001] = 1;;
	  e12[i_0020] = 3;
	}
      else if (z10 && !z6)
	{
	  e4[i_5000] = 1;
	  e4[i_2200] = 1;
	  
	  e6[i_2200] = 2;
	  e6[i_5000] = 1;
	  
	  e12[i_2001] = 1;
	  e12[i_2200] = 3;
	  e12[i_5000] = 2;
	}
      else /* z10 and z6 */	    
	{
	  e4[i_5000] = 1;
	  
	  e12[i_5000] = 2;
	  e12[i_2001] = 1;
	}
    }

  else if (wt == 30)
    {
      if (!z4 && !z6)
	{
	  SET_VEC_30(e4, 1, 1, 0, 0, 0, 0, 0);
	  SET_VEC_30(e6, 2, 1, 0, 0, 0, 0, 0);
	  SET_VEC_30(e10, 4, 0, 0, 1, 0, 0, 0);
	  SET_VEC_30(e12, 3, 2, 0, 0, 1, 0, 0);
	}
      else if (!z4 && !z10)
	{
	  SET_VEC_30(e4, 0, 0, 1, 1, 0, 0, 0);
	  SET_VEC_30(e6, 0, 0, 1, 1, 0, 1, 0);
	  SET_VEC_30(e10, 0, 0, 3, 2, 0, 0, 0);
	  SET_VEC_30(e12, 0, 0, 3, 2, 0, 0, 1);
	}
      else if (!z6 && !z10)
	{
	  SET_VEC_30(e4, 1, 0, 0, 0, 0, 1, 0);
	  SET_VEC_30(e6, 2, 0, 1, 0, 0, 0, 0);
	  SET_VEC_30(e10, 3, 0, 2, 0, 0, 0, 0);
	  SET_VEC_30(e12, 3, 0, 2, 0, 1, 0, 0);
	}
      else if (z4 && z10)
	{
	  e6[i_0500] = 1;
	  
	  e12[i_0500] = 1;
	  e12[i_0301] = 1;
	}
      else
	{
	  flint_printf("(igusa_from_monomials) Error: unsupported z4, z6, z10, z12 in weight 30: %i %i %i %i\n", z4, z6, z10, z12);
	  fflush(stdout);
	  flint_abort();
	}
    }

  else if (wt == 60)
    {
      if (z4 && z6 && z10)
	{
	  e12[i_0005] = 1;
	}
      else if (z4 && z6 && z12)
	{
	  e10[i_0060] = 1;
	}      
      else if (!z4 && !z10 && !z12)
	{
	  SET_VEC_60(e4, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0);
	  SET_VEC_60(e6, 1, 0, 0, 0, 3, 0, 0, 1, 1, 0);
	  SET_VEC_60(e10, 1, 4, 0, 0, 0, 0, 0, 5, 0, 0);
	  SET_VEC_60(e12, 3, 1, 0, 0, 5, 0, 0, 3, 0, 0);

	}
      else if (z4 && !z6 && !z10 && !z12)
	{
	  SET_VEC_60(e6, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0);
	  SET_VEC_60(e10, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0);
	  SET_VEC_60(e12, 0, 3, 2, 1, 0, 0, 0, 0, 0, 0);
	}
      else if (z10 && !z4 && !z6 && !z12)
	{
	  SET_VEC_60(e4, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0);
	  SET_VEC_60(e6, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0);
	  SET_VEC_60(e12, 0, 2, 0, 2, 1, 0, 0, 0, 0, 1);
	}

      else if (z12 && !z4 && !z10)
	{
	  SET_VEC_60(e4, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0);
	  SET_VEC_60(e6, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0);
	  SET_VEC_60(e10, 3, 0, 0, 0, 0, 0, 0, 2, 0, 0);
	}

      else if (z4 && z6 && !z10 && !z12)
	{
	  SET_VEC_60(e10, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0);
	  SET_VEC_60(e12, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0);
	}

      else if (z4 && z10 && !z6 && !z12)
	{
	  SET_VEC_60(e6, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
	  SET_VEC_60(e12, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0);
	}

      else if (z4 && z12 && !z6 && !z10)
	{
	  SET_VEC_60(e6, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0);
	  SET_VEC_60(e10, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0);
	}

      else if (z6 && z10 && !z4 && !z12)
	{
	  SET_VEC_60(e4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
	  SET_VEC_60(e12, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0);
	}	
      else
	{
	  flint_printf("(igusa_from_monomials) Error: unsupported z4, z6, z10, z12 in weight 60: %i %i %i %i\n", z4, z6, z10, z12);
	  fflush(stdout);
	  flint_abort();
	}
    }
  
  else
    {
      flint_printf("(igusa_from_monomials) Error: weight %wd not implemented\n", wt);
      fflush(stdout);
      flint_abort();
    }
}
