#ifndef MACROS_H__
#define MACROS_H__

#include "constants.h"

#ifndef NDIM
#error NDIM needs to be defined
#endif

#if ((NDIM!=1)&&(NDIM!=2)&&(NDIM!=3))
#error NDIM needs to be in range 1-3
#endif


#if NDIM==1
// Make a comma-separated list of NDIM elements from a list of 3
#define LIST_NDIM(x,y,z) x

// Make a space-separated list of NDIM elements from a list of 3
#define CONCAT_NDIM(x,y,z) x

// Pick the NDIM element of a list of 3
#define SELECT_NDIM(x,y,z) x

#elif NDIM==2
#define LIST_NDIM(x,y,z) x,y
#define CONCAT_NDIM(x,y,z) x y
#define SELECT_NDIM(x,y,z) y
#elif NDIM==3
#define LIST_NDIM(x,y,z) x,y,z
#define CONCAT_NDIM(x,y,z) x y z
#define SELECT_NDIM(x,y,z) z
#endif

#endif
