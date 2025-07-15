
#include "igusa.h"

/* See igusa_base_exps and igusa_from_monomials_exps */

/* Weight 20 */
static slong i_0020 = 0;
//static slong i_1110 = 1;
static slong i_2001 = 2;
static slong i_5000 = 3;
static slong i_2200 = 4;

/* Weight 30 */
static slong i_0500 = 0;
static slong i_3300 = 1;
static slong i_0030 = 2;
static slong i_5010 = 3;
static slong i_0301 = 4;
//static slong i_1120 = 5;
static slong i_2011 = 6;

/* Weight 60 */
static slong i_0060 = 0;
static slong i_0005 = 1;
//static slong i_0530 = 2;
//static slong i_0204 = 3;
static slong i_3004 = 4;
static slong i_0X00 = 5;
static slong i_3800 = 6;
static slong i_5040 = 7;
//static slong i_1150 = 8;
//static slong i_3203 = 9;

void igusa_from_monomials_zeroes(int* z4, int* z6, int* z10, int* z12,
				 fmpz* M, slong wt)
{
  if (wt == 20)
    {
      if (fmpz_is_zero(&M[i_5000]))
	{
	  flint_printf("(igusa_from_monomials_zeroes) Error: psi4 must be nonzero if weight 20 is used\n");
	  fflush(stdout);
	  flint_abort();
	}
      *z10 = fmpz_is_zero(&M[i_0020]);
      *z6 = fmpz_is_zero(&M[i_2200]);
      *z12 = fmpz_is_zero(&M[i_2001]);
      *z4 = 0;
    }

  else if (wt == 30)
    {
      *z6 = fmpz_is_zero(&M[i_0500]);
      *z10 = fmpz_is_zero(&M[i_0030]);
      if (!*z6)
	{
	  *z4 = fmpz_is_zero(&M[i_3300]);
	  *z12 = fmpz_is_zero(&M[i_0301]);
	}
      else /* z6; then z4 and z10 have to be 0 */
	{
	  if (fmpz_is_zero(&M[i_5010]) || *z10)
	    {	      
	      flint_printf("(igusa_from_monomials_zeroes) Error: psi4, chi10 must be nonzero if psi6 = 0 when weight 30 is used\n");
	      fflush(stdout);
	      flint_abort();
	    }
	  *z12 = fmpz_is_zero(&M[i_2011]);
	  *z4 = 0;
	}
    }

  else if (wt == 60)
    {
      *z6 = fmpz_is_zero(&M[i_0X00]);
      *z10 = fmpz_is_zero(&M[i_0060]);
      *z12 = fmpz_is_zero(&M[i_0005]);
      if (!*z6) *z4 = fmpz_is_zero(&M[i_3800]);
      else if (!*z10) *z4 = fmpz_is_zero(&M[i_5040]);
      else if (!*z12) *z4 = fmpz_is_zero(&M[i_3004]);
      else *z4 = 0;      
    }    
}
