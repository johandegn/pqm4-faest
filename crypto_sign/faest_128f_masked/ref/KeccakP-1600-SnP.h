/* SHA3 configuration. */
#include "config.h"
#include "KeccakP-1600-SnP-mpc.h"

#ifdef STM32F4
#include "KeccakP-1600-SnP-armv7m.h"
#else
#include "KeccakP-1600-SnP-opt64.h"
#endif
