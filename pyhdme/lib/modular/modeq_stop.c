
#include "modular.h"

int modeq_stop(int res, slong prec)
{
  int v = get_modeq_verbose();
  int stop;
  
  if (res)
    {
      if (v) flint_printf("(modeq_stop) Success at current working precision\n", prec);
      stop = 1;
    }
  else if (prec > MODEQ_MAX_PREC)    
    {
      flint_printf("(modeq_stop) Reached maximal allowed precision %wd, abandon.\n", MODEQ_MAX_PREC);
      stop = 1;      
    }
  else stop = 0;

  fflush(stdout);
  return stop;
}
