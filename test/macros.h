#ifndef MACROS_H__
#define MACROS_H__

#include "Bittree_constants.h"

#ifndef BTDIM
#error BTDIM needs to be defined
#endif

#if ((BTDIM!=1)&&(BTDIM!=2)&&(BTDIM!=3))
#error BTDIM needs to be in range 1-3
#endif


#if BTDIM==1
// Make a comma-separated list of BTDIM elements from a list of 3
#define LIST_NDIM(x,y,z) x

// Make a space-separated list of BTDIM elements from a list of 3
#define CONCAT_NDIM(x,y,z) x

// Pick the BTDIM element of a list of 3
#define SELECT_NDIM(x,y,z) x

#elif BTDIM==2
#define LIST_NDIM(x,y,z) x,y
#define CONCAT_NDIM(x,y,z) x y
#define SELECT_NDIM(x,y,z) y
#elif BTDIM==3
#define LIST_NDIM(x,y,z) x,y,z
#define CONCAT_NDIM(x,y,z) x y z
#define SELECT_NDIM(x,y,z) z
#endif

#endif
