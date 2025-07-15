
#include "hecke.h"

void hecke_collect_verbose_start(slong nb)
{
  flint_printf("(hecke_collect) Computing theta constants (%wd)\n", nb);
  fflush(stdout);
}
