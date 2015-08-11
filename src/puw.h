#include <math.h>
#include <stdlib.h>

#ifndef HUGE
#define HUGE 1.0E+20  //doesn't matter what this number is, just want -HUGE < 0.
#endif

#define WRAP(x) (((x) > 0.5) ? ((x)-1.0) : (((x) <= -0.5) ? ((x)+1.0) : (x)))

#define MASK 		0x1
#define UNWRAPPED 	0x2
