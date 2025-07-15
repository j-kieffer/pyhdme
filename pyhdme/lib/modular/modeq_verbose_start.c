
#include "modular.h"

void modeq_verbose_start(slong prec)
{
  flint_printf("(modeq_start) Start new run at precision %wd\n", prec);
  fflush(stdout);
}
