#include "config.h"

#ifndef KECCAK_MASK_NONE

#if defined(STM32F4) || defined(STM32L4)
#include "KeccakP-1600-inplace-32bi-armv7m-le-gcc-mpc.c.i"
#else
#include "KeccakP-1600-opt64.c.i"
#include "KeccakP-1600-ref-mpc.c.i"
#endif

#else

#if defined(STM32F4) || defined(STM32L4)
#else
#include "KeccakP-1600-opt64.c.i"
#endif

#endif
