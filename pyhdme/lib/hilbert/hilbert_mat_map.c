
#include "hilbert.h"

void hilbert_mat_map(fmpz_mat_t eta, const fmpz_poly_mat_t m, slong delta)
{
  slong prec = 100;
  int success = 0;
  int b;
  arf_t radius;
  acb_t sqrtd;
  acb_t msqrtd;

  acb_mat_t Ri, Rt;
  acb_mat_t M, N, P;
  slong j, k;

  arf_init(radius);
  acb_mat_init(Ri, 2, 2);
  acb_mat_init(Rt, 2, 2);
  acb_mat_init(M, 4, 4);
  acb_mat_init(N, 4, 4);
  acb_mat_init(P, 4, 4);
  acb_init(sqrtd);
  acb_init(msqrtd);


  while (!success)
    {
      success = 1;
      
      acb_zero(sqrtd);
      arb_sqrt_ui(acb_realref(sqrtd), delta, prec);
      acb_neg(msqrtd, sqrtd);
      
      hilbert_R(Ri, delta, prec);
      acb_mat_transpose(Rt, Ri);
      acb_mat_inv(Ri, Ri, prec);
      
      acb_mat_zero(M);
      acb_mat_set_window(M, 0, 0, Rt);
      acb_mat_set_window(M, 2, 2, Ri);
      acb_mat_inv(P, M, prec);

      hilbert_sigma1(acb_mat_entry(N, 0, 0),
		     fmpz_poly_mat_entry(m, 0, 0), delta, prec);
      hilbert_sigma2(acb_mat_entry(N, 1, 1),
		     fmpz_poly_mat_entry(m, 0, 0), delta, prec);
      
      hilbert_sigma1(acb_mat_entry(N, 0, 2),
		     fmpz_poly_mat_entry(m, 0, 1), delta, prec);
      acb_div(acb_mat_entry(N, 0, 2), acb_mat_entry(N, 0, 2), sqrtd, prec);
      hilbert_sigma2(acb_mat_entry(N, 1, 3),
		     fmpz_poly_mat_entry(m, 0, 1), delta, prec);
      acb_div(acb_mat_entry(N, 1, 3), acb_mat_entry(N, 1, 3), msqrtd, prec);
      
      hilbert_sigma1(acb_mat_entry(N, 2, 0),
		     fmpz_poly_mat_entry(m, 1, 0), delta, prec);
      acb_mul(acb_mat_entry(N, 2, 0), acb_mat_entry(N, 2, 0), sqrtd, prec);
      hilbert_sigma2(acb_mat_entry(N, 3, 1),
		     fmpz_poly_mat_entry(m, 1, 0), delta, prec);
      acb_mul(acb_mat_entry(N, 3, 1), acb_mat_entry(N, 3, 1), msqrtd, prec);
      
      hilbert_sigma1(acb_mat_entry(N, 2, 2),
		     fmpz_poly_mat_entry(m, 1, 1), delta, prec);
      hilbert_sigma2(acb_mat_entry(N, 3, 3),
		     fmpz_poly_mat_entry(m, 1, 1), delta, prec);
      
      acb_mat_mul(M, M, N, prec);
      acb_mat_mul(M, M, P, prec);
      for (j = 0; j < 4; j++)
	{
	  for (k = 0; k < 4; k++)
	    {
	      b = acb_round(fmpz_mat_entry(eta, j, k), radius,
			    acb_mat_entry(M, j, k));
	      success = success && b;
	    }
	}
      prec *= 2;
      if (!success)
	{
	  flint_printf("(hilbert_mat_map) Info: increasing precision\n");
	  acb_mat_printd(M, 10); flint_printf("\n");
	}
    }

  arf_clear(radius);
  acb_mat_clear(Ri);
  acb_mat_clear(Rt);
  acb_mat_clear(M);
  acb_mat_clear(N);
  acb_mat_clear(P);
  acb_clear(sqrtd);
  acb_clear(msqrtd);
}
