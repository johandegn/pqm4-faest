/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef RANDOMNESS_H
#define RANDOMNESS_H

#include "macros.h"

#include <stddef.h>
#include <stdint.h>

FAEST_BEGIN_C_DECL

int rand_bytes(uint8_t* dst, size_t num_bytes);

FAEST_END_C_DECL

#if defined(STM32F4) || defined(MUPQ_NAMESPACE)

/* If we are building under pqm4 for either host or board, use external PRNG. */
#define HAVE_RANDOMBYTES

/* Also define another RNG for mask generation. */
int rand_mask(uint8_t* dst, size_t len);

#else

#define rand_mask rand_bytes

#endif

#endif /* RANDOMNESS_H */