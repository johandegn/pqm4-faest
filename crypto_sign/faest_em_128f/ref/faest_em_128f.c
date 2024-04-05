/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest_em_128f.h"
#include "compat.h"
#include "randomness.h"
#include "owf.h"
#include "instances.h"
#include "faest.h"
#include "parameters.h"

#include <stdlib.h>
#include <string.h>

// memory layout of the public key: OWF input || OWF output
#define PK_INPUT(pk) (pk)
#define PK_OUTPUT(pk) (&pk[32 / 2])

// memory layout of the secret key: OWF input || OWF key
#define SK_INPUT(sk) (sk)
#define SK_KEY(sk) (&sk[32 / 2])

int FAEST_CALLING_CONVENTION faest_em_128f_keygen(uint8_t* pk, uint8_t* sk) {
  if (!pk || !sk) {
    return -1;
  }

  bool done = false;
  while (!done) {
    rand_bytes(sk, 32);
    // declassify OWF input
    faest_declassify(SK_INPUT(sk), 32 / 2);
    done = faest_em_128f_owf(SK_KEY(sk), SK_INPUT(sk), PK_OUTPUT(pk));
    faest_declassify(&done, sizeof(done));
  }
  memcpy(PK_INPUT(pk), SK_INPUT(sk), 32 / 2);

  // declassify public key
  faest_declassify(pk, 32);
  return 0;
}

int FAEST_CALLING_CONVENTION faest_em_128f_validate_keypair(const uint8_t* pk, const uint8_t* sk) {
  if (!sk || !pk) {
    return -1;
  }

  uint8_t pk_check[32];
  if (!faest_em_128f_owf(SK_KEY(sk), SK_INPUT(sk), PK_OUTPUT(pk_check))) {
    // zero bytes in SubBytes input
    return 1;
  }
  memcpy(PK_INPUT(pk_check), SK_INPUT(sk), 32 / 2);

  return faest_timingsafe_bcmp(pk_check, pk, sizeof(pk_check)) == 0 ? 0 : 2;
}

int FAEST_CALLING_CONVENTION faest_em_128f_sign_with_randomness(const uint8_t* sk, const uint8_t* message, size_t message_len, const uint8_t* rho, size_t rho_len, uint8_t* signature, size_t* signature_len) {
  if (!sk || !signature || !signature_len || *signature_len < FAEST_EM_128F_SIGNATURE_SIZE || (!rho && rho_len)) {
    return -1;
  }

  uint8_t owf_output[32 / 2];
  if (!faest_em_128f_owf(SK_KEY(sk), SK_INPUT(sk), owf_output)) {
    // invalid key
    return -1;
  }
  // declassify OWF output
  faest_declassify(owf_output, sizeof(owf_output));

  const faest_paramset_t params = faest_get_paramset(FAEST_EM_128F);
  faest_sign(signature, message, message_len, SK_KEY(sk), SK_INPUT(sk), owf_output, rho, rho_len, &params);
  *signature_len = FAEST_EM_128F_SIGNATURE_SIZE;

  return 0;
}

int FAEST_CALLING_CONVENTION faest_em_128f_sign(const uint8_t* sk, const uint8_t* message, size_t message_len, uint8_t* signature, size_t* signature_len) {
  if (!sk || !signature || !signature_len || *signature_len < FAEST_EM_128F_SIGNATURE_SIZE) {
    return -1;
  }

  uint8_t rho[FAEST_EM_128F_LAMBDA / 8];
  rand_bytes(rho, sizeof(rho));

  return faest_em_128f_sign_with_randomness(sk, message, message_len, rho, sizeof(rho), signature, signature_len);
}

int FAEST_CALLING_CONVENTION faest_em_128f_verify(const uint8_t* pk, const uint8_t* message, size_t message_len, const uint8_t* signature, size_t signature_len) {
  if (!pk || !signature || signature_len != FAEST_EM_128F_SIGNATURE_SIZE) {
    return -1;
  }

  const faest_paramset_t params = faest_get_paramset(FAEST_EM_128F);
  return faest_verify(message, message_len, signature, PK_INPUT(pk), PK_OUTPUT(pk), &params);
}

void FAEST_CALLING_CONVENTION faest_em_128f_clear_private_key(uint8_t* key) {
  faest_explicit_bzero(key, FAEST_EM_128F_PRIVATE_KEY_SIZE);
}

// vim: ft=c
