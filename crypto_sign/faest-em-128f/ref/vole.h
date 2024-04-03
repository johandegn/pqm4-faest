/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_VOLE_H
#define FAEST_VOLE_H

#include "vc.h"
#include <stdbool.h>

FAEST_BEGIN_C_DECL

// k_b is at most 12, so chalout needs to point to an array of at most 12 bytes
int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout);

void stream_vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, uint8_t* hcom, stream_vec_com_t* sVecCom, uint8_t* c,
                 uint8_t* u, uint8_t** v);

void stream_vole_reconstruct(const uint8_t* iv, const uint8_t* chall, const uint8_t* const* pdec,
                             const uint8_t* const* com_j, uint8_t* hcom, uint8_t** q, unsigned int ellhat,
                             const faest_paramset_t* params);

FAEST_END_C_DECL

#endif
