#ifndef FAEST_FAEST_AES_MASKED_H
#define FAEST_FAEST_AES_MASKED_H

#include <stdint.h>
#include <assert.h>

#include "instances.h"
#include "aes.h"
#include "vbb.h"

void aes_prove_masked(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
               const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params);

#endif