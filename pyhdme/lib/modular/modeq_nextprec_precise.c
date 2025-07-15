
#include "modular.h"

slong modeq_nextprec_precise(slong current_prec, slong gap)
{
  slong res = current_prec + gap + 25;
  return 100 * (res/100 + 1);
}
