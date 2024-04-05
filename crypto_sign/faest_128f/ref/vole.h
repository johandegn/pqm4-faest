/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_VOLE_H
#define FAEST_VOLE_H

#include <stdbool.h>
#include <assert.h>

#include "vc.h"
#include "macros.h"

FAEST_BEGIN_C_DECL

typedef enum { EXCLUDE_U_HCOM_C, EXCLUDE_HCOM, EXCLUDE_V, EXCLUDE_Q, INCLUDE_ALL } vole_mode_t;

typedef struct sign_vole_mode_ctx_t {
  vole_mode_t mode;
  uint8_t* v;
  uint8_t* u;
  uint8_t* hcom;
  uint8_t* c;
} sign_vole_mode_ctx_t;

ATTR_CONST ATTR_ALWAYS_INLINE inline sign_vole_mode_ctx_t vole_mode_all_sign(uint8_t* v, uint8_t* u,
                                                               uint8_t* hcom, uint8_t* c) {
  assert(v != NULL);
  assert(u != NULL);
  assert(hcom != NULL);
  assert(c != NULL);
  sign_vole_mode_ctx_t ctx = {INCLUDE_ALL, v, u, hcom, c};
  return ctx;
}

ATTR_CONST ATTR_ALWAYS_INLINE inline sign_vole_mode_ctx_t vole_mode_u_hcom_c(uint8_t* u, uint8_t* hcom,
                                                                    uint8_t* c) {
  assert(u != NULL);
  assert(hcom != NULL);
  assert(c != NULL);
  sign_vole_mode_ctx_t ctx = {EXCLUDE_V, NULL, u, hcom, c};
  return ctx;
}

ATTR_CONST ATTR_ALWAYS_INLINE inline sign_vole_mode_ctx_t vole_mode_v(uint8_t* v) {
  assert(v != NULL);
  sign_vole_mode_ctx_t ctx = {EXCLUDE_U_HCOM_C, v, NULL, NULL, NULL};
  return ctx;
}

typedef struct verify_vole_mode_ctx_t {
  vole_mode_t mode;
  uint8_t* q;
  uint8_t* hcom;
} verify_vole_mode_ctx_t;

ATTR_CONST ATTR_ALWAYS_INLINE inline verify_vole_mode_ctx_t vole_mode_all_verify(uint8_t* q, uint8_t* hcom) {
  assert(q != NULL);
  assert(hcom != NULL);
  verify_vole_mode_ctx_t ctx = {INCLUDE_ALL, q, hcom};
  return ctx;
}

ATTR_CONST ATTR_ALWAYS_INLINE inline verify_vole_mode_ctx_t vole_mode_q(uint8_t* q) {
  assert(q != NULL);
  verify_vole_mode_ctx_t ctx = {EXCLUDE_HCOM, q, NULL};
  return ctx;
}

ATTR_CONST ATTR_ALWAYS_INLINE inline verify_vole_mode_ctx_t vole_mode_hcom(uint8_t* hcom) {
  assert(hcom != NULL);
  verify_vole_mode_ctx_t ctx = {EXCLUDE_Q, NULL, hcom};
  return ctx;
}

// k_b is at most 12, so chalout needs to point to an array of at most 12 bytes
int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout);

// Signer
void partial_vole_commit_cmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             unsigned int chunk_start, unsigned int chunk_end,
                             sign_vole_mode_ctx_t vole_mode, const faest_paramset_t* params);

void partial_vole_commit_rmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int start,
                             unsigned int len, const faest_paramset_t* params, uint8_t* v);

// Verifier
void partial_vole_reconstruct_cmo(const uint8_t* iv, const uint8_t* chall, const uint8_t* const* pdec, 
                                  const uint8_t* const* com_j, unsigned int ellhat,
                                  unsigned int start, unsigned int len,
                                  verify_vole_mode_ctx_t vole_mode, 
                                  const faest_paramset_t* params);

void partial_vole_reconstruct_rmo(const uint8_t* iv, const uint8_t* chall,
                                  const uint8_t* const* pdec, const uint8_t* const* com_j,
                                  uint8_t* q, unsigned int ellhat, const faest_paramset_t* params,
                                  unsigned int start, unsigned int len);

FAEST_END_C_DECL

#endif
