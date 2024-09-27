/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef CRYPTO_SIGN_EM_192F_H
#define CRYPTO_SIGN_EM_192F_H

#define CRYPTO_SECRETKEYBYTES 48
#define CRYPTO_PUBLICKEYBYTES 48
#define CRYPTO_BYTES 13912
#define CRYPTO_ALGNAME "faest_em_192f"

#include <stdlib.h>

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk);
int crypto_sign(unsigned char* sm, size_t* smlen, const unsigned char* m,
                size_t mlen, const unsigned char* sk);
int crypto_sign_open(unsigned char* m, size_t* mlen, const unsigned char* sm,
                     size_t smlen, const unsigned char* pk);

#endif

// vim: ft=c
