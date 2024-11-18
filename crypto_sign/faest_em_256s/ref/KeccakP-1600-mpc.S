#include "config.h"

#ifndef KECCAK_MASK_NONE

#if defined(STM32F4) || defined(STM32L4)
#include "KeccakP-1600-inplace-32bi-armv7m-le-gcc.S.i"
#include "KeccakP-1600-inplace-32bi-armv7m-le-gcc-mpc.S.i"
#endif

#else

#if defined(STM32F4) || defined(STM32L4)
#include "KeccakP-1600-inplace-32bi-armv7m-le-gcc.S.i"
#endif

#endif
