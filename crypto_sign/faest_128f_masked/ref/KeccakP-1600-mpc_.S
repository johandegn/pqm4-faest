#include "config.h"

#ifndef KECCAK_MASK_NONE
//#warning "Can i see this warning Here"
#ifdef STM32F4
//#warning "This one also works for STM32F4"
#include "KeccakP-1600-inplace-32bi-armv7m-le-gcc.S.i"
#include "KeccakP-1600-inplace-32bi-armv7m-le-gcc-mpc.S.i"
#endif

#endif

