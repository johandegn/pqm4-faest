/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef SUPERCOP
#include "crypto_sign.h"
#else
#include "api.h"
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest_192s.h"

#include <string.h>

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk) {
  return faest_192s_keygen(pk, sk);
}

int crypto_sign(unsigned char* sm, size_t* smlen, const unsigned char* m,
                size_t mlen, const unsigned char* sk) {
  *smlen = mlen + FAEST_192S_SIGNATURE_SIZE;
  memmove(sm, m, mlen);

  size_t signature_len = FAEST_192S_SIGNATURE_SIZE;
  return faest_192s_sign(sk, sm, mlen, sm + mlen, &signature_len);
}

int crypto_sign_open(unsigned char* m, size_t* mlen, const unsigned char* sm,
                     size_t smlen, const unsigned char* pk) {
  if (smlen < FAEST_192S_SIGNATURE_SIZE) {
    // signature too short
    return -1;
  }
  size_t m_length = smlen - FAEST_192S_SIGNATURE_SIZE;
  if (faest_192s_verify(pk, sm, m_length, sm + m_length, FAEST_192S_SIGNATURE_SIZE)) {
    return -1;
  }

  *mlen = m_length;
  memmove(m, sm, m_length);
  return 0;
}

// vim: ft=c
