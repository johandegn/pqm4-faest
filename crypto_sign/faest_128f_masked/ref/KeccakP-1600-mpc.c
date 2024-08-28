#include "config.h"

#ifndef KECCAK_MASK_NONE

#ifdef STM32F4
#include "KeccakP-1600-inplace-32bi-armv7m-le-gcc-mpc.c.i"
#else
#include "KeccakP-1600-opt64.c.i"
#include "KeccakP-1600-ref-mpc.c.i"
#endif

#else

#ifdef STM32F4
#else
#include "KeccakP-1600-opt64.c.i"
#endif

#endif
