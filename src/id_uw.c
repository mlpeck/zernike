#include "puw.h"

#define unwrap(dp) {          \
puw[kn] = puw[kh] + dp;   \
flags[kn] |= UNWRAPPED;   \
uw[kn] = MASK;            \
todo[ntodo++] = kn;       \
}

void id_dxy_uw(int *nr, int *nc, double *phase, double *mask, 
               double *dx, double *dy, double *puw, int *uw) 
{
  
  int size=(*nr)*(*nc);
  char *flags = malloc(size*sizeof(char));
  int *todo =malloc(size*sizeof(int));
  int k, kh, kn, x, y, ntodo=0;
  
  for (k=0; k<size; k++) flags[k] = isnan(mask[k]);
  for (k=0; k<size; k++) {
    if (!flags[k]) {
      kn = kh = k;
      break;
    }
  }
  
  puw[kh] = phase[kh];
  unwrap(0.0);
  
  while (ntodo) {
    kh = todo[--ntodo];
    x = kh % (*nr);
    y = kh / (*nr);
    kn = kh-1;
    if (x > 0 && !flags[kn])
      unwrap(-dx[kn]);
    kn = kh+1;
    if (x < (*nr-1) && !flags[kn])
      unwrap(dx[kh]);
    kn = kh - (*nr);
    if (y > 0 && !flags[kn])
      unwrap(-dy[kn]);
    kn = kh + (*nr);
    if (y < (*nc-1) && !flags[kn])
      unwrap(dy[kh]);
  }
}

// same routine as above, but without dx and dy passed.

#define unwrap2(dp) {          \
puw[kn] = puw[kh] + dp;   \
flags[kn] |= UNWRAPPED;   \
todo[ntodo++] = kn;       \
}

void id_uw(int *nr, int *nc, double *phase, double *puw)
{
  int size=(*nr)*(*nc);
  char *flags = malloc(size*sizeof(char));
  int *todo =malloc(size*sizeof(int));
  int k, kh, kn, x, y, ntodo=0;
  
  for (k=0; k<size; k++) flags[k] = isnan(phase[k]);
  for (k=0; k<size; k++) {
    if (!flags[k]) {
      kn = kh = k;
      break;
    }
  }
  puw[kh] = phase[kh];
  unwrap2(0.0);
  
  while (ntodo) {
    kh = todo[--ntodo];
    x = kh % (*nr);
    y = kh / (*nr);
    kn = kh-1;
    if (x > 0 && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
    kn = kh+1;
    if (x < (*nr-1) && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
    kn = kh - (*nr);
    if (y > 0 && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
    kn = kh + (*nr);
    if (y < (*nc-1) && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
  }
}
