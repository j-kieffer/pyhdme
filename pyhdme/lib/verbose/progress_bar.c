#include <stdio.h> // printf
#include <string.h> // memset
#include <unistd.h> // stdout
#include <flint/flint.h> // slong
//#include <omp.h>
//What I really wanted ot use '▰' '▱'

#define bar_length 30
void progress_bar(slong current, slong total, const char prefix[]) {
  static char donebar[bar_length];
  static char notdonebar[bar_length];
  float proportion = ((float) current)/total;
  int done = proportion * bar_length;
  int notdone = bar_length - done;
  memset(donebar, '#', done);
  donebar[done] = '\0';
  memset(notdonebar, '-', notdone);
  notdonebar[notdone] = '\0';
  printf("\r%s [%s%s] %.2f%% (%ld/%ld)", prefix, donebar, notdonebar, 100*proportion, current, total);
  if (current == total) printf("\n");
  fflush(stdout);
}
