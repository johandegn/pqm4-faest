#ifndef FAEST_AES_MASKED_H
#define FAEST_AES_MASKED_H

#include <stdint.h>
#include "instances.h"

uint8_t* aes_extend_witness_masked(const uint8_t* key, const uint8_t* in, const faest_paramset_t* params, uint8_t* w);

#endif