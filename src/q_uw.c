#include "puw.h"

/* Heap management code courtesy Steve K. */

#define SWAP(a, i, j)                           \
{                                           \
  int t = a[i];                           \
  a[i] = a[j];                            \
  a[j] = t;                               \
}

void todo_push (int ndx, double *qual, int *todo, int *end)
{
  int child;
  todo[*end] = ndx;
  child = (*end)++;
  while (child > 0) {
    int parent = (child-1) / 2;
    if (qual[todo[parent]] < qual[todo[child]]) {
      SWAP (todo, parent, child);
      child = parent;
    }
    else
      break;
  }
}

int todo_pop(double *qual, int *todo, int *end)
{
  int result = todo[0], root;
  --(*end);
  SWAP (todo, 0, *end);
  root = 0;
  while (root*2+1 < *end) {
    int child = root*2+1; // left child
    if (child+1 < *end && qual[todo[child]] < qual[todo[child+1]])
      ++child;
    if (qual[todo[root]] < qual[todo[child]]) {
      SWAP (todo, root, child);
      root = child;
    }
    else
      break;
  }
  return (result);
}


void q_uw(int *nr, int *nc, double *phase, double *qual, double *puw) {
  
  int size=(*nr)*(*nc);
  char *flags = malloc(size*sizeof(char));
  int *todo =malloc(size*sizeof(int));
  double qm= -HUGE;
  int k, kh, kn, x, y, end=0;
  int nk[4] = {-1, 1, -(*nr), (*nr)};
  int bxy[4] = {0, -(*nr-1), 0, -(*nc-1)}, cxy[4];
  
  for (k=0; k<size; k++) {
    flags[k] = isnan(phase[k]);
    if (qual[k] > qm && !flags[k])
      qm = qual[kh = k];
  };
  
  // start at the highest quality point.
  
  puw[kh] = phase[kh];
  flags[kh] |= UNWRAPPED;
  todo_push(kh, qual, todo, &end);
  
  while (end) {
    kh = todo_pop(qual, todo, &end);
    x = kh % (*nr);
    y = kh / (*nr);
    cxy[1] = -(cxy[0]=x);
    cxy[3] = -(cxy[2]=y);
    
    //these are the neighbors of the currently selected point.
    
    for (k=0; k<4; k++) {
      kn = kh+nk[k];
      if ((cxy[k] > bxy[k]) && !flags[kn]) {
        puw[kn] = puw[kh] + WRAP(phase[kn]-phase[kh]);
        flags[kn] |= UNWRAPPED;
        todo_push(kn, qual, todo, &end);
      }
    }
  }
}
