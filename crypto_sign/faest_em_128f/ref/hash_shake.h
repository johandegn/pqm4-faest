#include "config.h"
/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef HASH_SHAKE_H
#define HASH_SHAKE_H

#include <stdint.h>
#include <stdio.h>
#include <alloca.h>

#include "macros.h"
#include "endian_compat.h"
#include "randomness.h"

#if defined(WITH_SHAKE_S390_CPACF)
/* use the KIMD/KLMD instructions from CPACF for SHAKE support on S390 */
#include "sha3/s390_cpacf.h"
#elif 0 //defined(OQS) || defined(PQCLEAN)
#if defined(OQS)
/* use OQS's SHAKE implementation */
#include <oqs/sha3.h>
#include <oqs/sha3x4.h>
#elif defined(PQCLEAN)
/* PQClean's SHAKE implementation
 *
 * PQClean current does not expose the AVX2-optimized version of Keccakx4.
 */
#include <fips202.h>
#define OQS_SHA3_shake128_inc_ctx shake128incctx
#define OQS_SHA3_shake128_inc_init shake128_inc_init
#define OQS_SHA3_shake128_inc_absorb shake128_inc_absorb
#define OQS_SHA3_shake128_inc_finalize shake128_inc_finalize
#define OQS_SHA3_shake128_inc_squeeze shake128_inc_squeeze
#define OQS_SHA3_shake128_inc_ctx_release shake128_inc_ctx_release
#define OQS_SHA3_shake256_inc_ctx shake256incctx
#define OQS_SHA3_shake256_inc_init shake256_inc_init
#define OQS_SHA3_shake256_inc_absorb shake256_inc_absorb
#define OQS_SHA3_shake256_inc_finalize shake256_inc_finalize
#define OQS_SHA3_shake256_inc_squeeze shake256_inc_squeeze
#define OQS_SHA3_shake256_inc_ctx_release shake256_inc_ctx_release
#endif

typedef struct hash_context_oqs_s {
  union {
    OQS_SHA3_shake128_inc_ctx shake128_ctx;
    OQS_SHA3_shake256_inc_ctx shake256_ctx;
  };
  unsigned char shake256;
} hash_context;

/**
 * Initialize hash context based on the security parameter. If the security parameter is 128,
 * SHAKE128 is used, otherwise SHAKE256 is used.
 */
static inline void hash_init(hash_context* ctx, unsigned int security_param) {
  if (security_param == 128) {
    OQS_SHA3_shake128_inc_init(&ctx->shake128_ctx);
    ctx->shake256 = 0;
  } else {
    OQS_SHA3_shake256_inc_init(&ctx->shake256_ctx);
    ctx->shake256 = 1;
  }
}

static inline void hash_update(hash_context* ctx, const uint8_t* data, size_t size) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_inc_absorb(&ctx->shake256_ctx, data, size);
  } else {
    OQS_SHA3_shake128_inc_absorb(&ctx->shake128_ctx, data, size);
  }
}

static inline void hash_final(hash_context* ctx) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_inc_finalize(&ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_inc_finalize(&ctx->shake128_ctx);
  }
}

static inline void hash_squeeze(hash_context* ctx, uint8_t* buffer, size_t buflen) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_inc_squeeze(buffer, buflen, &ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_inc_squeeze(buffer, buflen, &ctx->shake128_ctx);
  }
}

static inline void hash_clear(hash_context* ctx) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_inc_ctx_release(&ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_inc_ctx_release(&ctx->shake128_ctx);
  }
}
#if defined(OQS)
/* Instances that work with 4 states in parallel using the base implementation. */
typedef struct hash_context_x4_oqs_s {
  union {
    OQS_SHA3_shake128_x4_inc_ctx shake128_ctx;
    OQS_SHA3_shake256_x4_inc_ctx shake256_ctx;
  };
  unsigned int shake256;
} hash_context_x4;

static inline void hash_init_x4(hash_context_x4* ctx, unsigned int security_param) {
  if (security_param == 128) {
    OQS_SHA3_shake128_x4_inc_init(&ctx->shake128_ctx);
    ctx->shake256 = 0;
  } else {
    OQS_SHA3_shake256_x4_inc_init(&ctx->shake256_ctx);
    ctx->shake256 = 1;
  }
}

static inline void hash_update_x4(hash_context_x4* ctx, const uint8_t** data, size_t size) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_x4_inc_absorb(&ctx->shake256_ctx, data[0], data[1], data[2], data[3], size);
  } else {
    OQS_SHA3_shake128_x4_inc_absorb(&ctx->shake128_ctx, data[0], data[1], data[2], data[3], size);
  }
}

static inline void hash_update_x4_4(hash_context_x4* ctx, const uint8_t* data0,
                                    const uint8_t* data1, const uint8_t* data2,
                                    const uint8_t* data3, size_t size) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_x4_inc_absorb(&ctx->shake256_ctx, data0, data1, data2, data3, size);
  } else {
    OQS_SHA3_shake128_x4_inc_absorb(&ctx->shake128_ctx, data0, data1, data2, data3, size);
  }
}

static inline void hash_update_x4_1(hash_context_x4* ctx, const uint8_t* data, size_t size) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_x4_inc_absorb(&ctx->shake256_ctx, data, data, data, data, size);
  } else {
    OQS_SHA3_shake128_x4_inc_absorb(&ctx->shake128_ctx, data, data, data, data, size);
  }
}

static inline void hash_init_prefix_x4(hash_context_x4* ctx, unsigned int security_parameter,
                                       const uint8_t prefix) {
  if (security_param == 128) {
    OQS_SHA3_shake128_x4_inc_init(&ctx->shake128_ctx);
    OQS_SHA3_shake128_x4_inc_absorb(&ctx->shake128_ctx, &prefix, &prefix, &prefix, &prefix,
                                    sizeof(prefix));
    ctx->shake256 = 0;
  } else {
    OQS_SHA3_shake256_x4_inc_init(&ctx->shake256_ctx);
    OQS_SHA3_shake256_x4_inc_absorb(&ctx->shake256_ctx, &prefix, &prefix, &prefix, &prefix,
                                    sizeof(prefix));
    ctx->shake256 = 1;
  }
}

static inline void hash_final_x4(hash_context_x4* ctx) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_x4_inc_finalize(&ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_x4_inc_finalize(&ctx->shake128_ctx);
  }
}

static inline void hash_squeeze_x4(hash_context_x4* ctx, uint8_t** buffer, size_t buflen) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_x4_inc_squeeze(buffer[0], buffer[1], buffer[2], buffer[3], buflen,
                                     &ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_x4_inc_squeeze(buffer[0], buffer[1], buffer[2], buffer[3], buflen,
                                     &ctx->shake128_ctx);
  }
}

static inline void hash_squeeze_x4_4(hash_context_x4* ctx, uint8_t* buffer0, uint8_t* buffer1,
                                     uint8_t* buffer2, uint8_t* buffer3, size_t buflen) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_x4_inc_squeeze(buffer0, buffer1, buffer2, buffer3, buflen,
                                     &ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_x4_inc_squeeze(buffer0, buffer1, buffer2, buffer3, buflen,
                                     &ctx->shake128_ctx);
  }
}

static inline void hash_clear_x4(hash_context_x4* ctx) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_x4_inc_ctx_release(&ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_x4_inc_ctx_release(&ctx->shake128_ctx);
  }
}
#endif
#else
#if !defined(SUPERCOP)
/* use SHAKE implementation in sha3/ */
#include "KeccakHash.h"
#include "KeccakMPCHash.h"
#else
/* use SUPERCOP implementation */
#include <libkeccak.a.headers/KeccakHash.h>
#endif

#if defined(WITH_KECCAK_X4)
/* use the Keccakx4 implementation */
#include "KeccakHashtimes4.h"
#endif

#ifdef KECCAK_MASK_ALL    /* When KECCAK_MASK_ALL is set, ALL calls to SHA-3 are masked. */
typedef KeccakMPC_HashInstance hash_context ATTR_ALIGNED(32);
typedef KeccakMPC_HashInstance masked_hash_context ATTR_ALIGNED(32);
#else                    /* Otherwise, we only mask those calls that use masked_hash_context */
typedef Keccak_HashInstance hash_context;
#ifdef KECCAK_MASK_NONE
typedef Keccak_HashInstance masked_hash_context;
#else
typedef KeccakMPC_HashInstance masked_hash_context ATTR_ALIGNED(32);
#endif
#endif

#ifdef KECCAK_MASK_NONE

static inline void array_xor(uint8_t* a, const uint8_t* b, const uint8_t* c, size_t len ) {
  size_t i;
  for( i=0; i< len; i++) {
    a[i] = b[i] ^ c[i];
  }
}

static inline void masked_hash_init(masked_hash_context *ctx, unsigned int security_param, size_t shareCount) {
  (void)shareCount;
  if (security_param == 128) {
    Keccak_HashInitialize_SHAKE128(ctx);
  } else {
    Keccak_HashInitialize_SHAKE256(ctx);
  }
}

static inline void masked_hash_update(masked_hash_context* ctx, const uint8_t** data, size_t byteLen) {
  uint8_t *buffer = (uint8_t *)alloca(byteLen);
  array_xor(buffer, data[0], data[1], byteLen);
  Keccak_HashUpdate(ctx, buffer, byteLen * 8);
}

static inline void masked_hash_final(masked_hash_context* ctx) {
  Keccak_HashFinal(ctx, NULL);
}

static inline void masked_hash_squeeze(masked_hash_context* ctx, uint8_t** digestShares, size_t byteLen) {
  uint8_t *buffer = (uint8_t *)alloca(byteLen);
  Keccak_HashSqueeze(ctx, buffer, byteLen * 8);
  memcpy(digestShares[0], buffer, byteLen);
  memset(digestShares[1], 0, byteLen);
}

static inline void masked_hash_update_uint16_le(masked_hash_context* ctx, uint16_t ** shares) {
  Keccak_HashUpdate(ctx, (const uint8_t*)shares[0], 16);
}

#else /* KECCAK_MASK_NONE */

static inline void masked_hash_init(masked_hash_context *ctx, unsigned int security_param, size_t shareCount) {
  if (shareCount != 2) {
    return;
  }
  if( security_param == 128) {
    KeccakMPC_HashInitialize_2SHARE_SHAKE128(ctx);
  } else {
    KeccakMPC_HashInitialize_2SHARE_SHAKE256(ctx);
  }
}

static inline void masked_hash_update(masked_hash_context* ctx, const uint8_t** data, size_t byteLen) {
  KeccakMPC_HashUpdate(ctx, (uint8_t**)data, byteLen * 8);
}

static inline void masked_hash_final(masked_hash_context* ctx) {
  KeccakMPC_HashFinal(ctx, NULL);
}

static inline void masked_hash_squeeze(masked_hash_context* ctx, uint8_t** digestShares, size_t byteLen) {
  KeccakMPC_HashSqueeze(ctx, digestShares, byteLen * 8);
#ifdef KECCAK_SNI_SECURE
  uint8_t *buffer = (uint8_t *)alloca(byteLen);
  rand_mask(buffer, byteLen);
  array_xor(digestShares[0], digestShares[0], buffer, byteLen);
  array_xor(digestShares[1], digestShares[1], buffer, byteLen);
#endif
}

static inline void masked_hash_update_uint16_le(masked_hash_context* ctx, uint16_t ** shares) {
  KeccakMPC_HashUpdate(ctx, (BitSequence**)shares, 16);
}

#define masked_hash_clear(ctx)
#endif /* KECCAK_MASK_NONE */

static inline void masked_hash_init_prefix(masked_hash_context* ctx, unsigned int security_param, const uint8_t prefix, size_t shareCount) {
  masked_hash_init(ctx, security_param, shareCount);
  uint8_t **shares = (uint8_t **)alloca(shareCount * sizeof(uint8_t**));
  uint8_t * slab = (uint8_t *)alloca(shareCount);
  shares[0] = slab;
  shares[0][0] = prefix;
  for(size_t i = 1; i < shareCount; ++i) {
    shares[i] = slab + i;
    shares[i][0] = 0;
  }
  masked_hash_update(ctx, (const uint8_t **)shares, sizeof(prefix));
}

#ifdef KECCAK_MASK_ALL

static inline void array_xor(uint8_t* a, const uint8_t* b, const uint8_t* c, size_t len ) {
  size_t i;
  for( i=0; i< len; i++) {
    a[i] = b[i] ^ c[i];
  }
}

/**
 * Initialize hash context based on the digest size used by FAEST.
 */
static inline void hash_init(hash_context* ctx, unsigned int security_param) {
  masked_hash_init(ctx, security_param, 2);
}

static inline void hash_update(hash_context* ctx, const uint8_t* data, size_t size) {
  const int shareCount = 2;
  uint8_t ** shares = (uint8_t **)alloca(shareCount * sizeof(uint8_t*));
  for(int i = 0; i < shareCount; i++) {
    shares[i] = (uint8_t *)alloca(size);
  }
  memcpy(shares[shareCount-1], data, size);

  //set pointers of "outer" array and fill last share with xor
  for(int i = 0; i < shareCount-1; i++) {
    rand_mask(shares[i], size);
    array_xor(shares[shareCount-1], shares[shareCount-1], shares[i], size);
  }

  masked_hash_update(ctx, (const uint8_t**)shares, size);
}

static inline void hash_final(hash_context* ctx) {
  masked_hash_final(ctx);
}

static inline void hash_squeeze(hash_context* ctx, uint8_t* buffer, size_t buflen) {
  uint8_t ** shares = (uint8_t **)alloca(ctx->shareCount * sizeof(uint8_t*));

  for(size_t i = 0; i < ctx->shareCount; i++) {
	shares[i] = (uint8_t *)alloca(buflen);
  }

  masked_hash_squeeze(ctx, shares, buflen);

  memset(buffer, 0, buflen);
  for(size_t i = 0; i < ctx->shareCount; i++) {
    array_xor(buffer, buffer, shares[i], buflen);
  }
}

static inline void hash_update_uint16_le(hash_context* ctx, uint16_t data) {
  uint16_t zero = 0;
  uint16_t *shares[2] = {&data, &zero};
  masked_hash_update_uint16_le(ctx, (uint16_t**) shares);
}

static inline void hash_init_prefix(hash_context* ctx, size_t digest_size, const uint8_t prefix) {
  masked_hash_init_prefix(ctx, digest_size, prefix, 2);
}

#define hash_clear(ctx)
#else /* !KECCAK_MASK_ALL */

/**
 * Initialize hash context based on the digest size used by Picnic. If the size is 32 bytes,
 * SHAKE128 is used, otherwise SHAKE256 is used.
 */
static inline void hash_init(hash_context* ctx, unsigned int security_param) {
  if (security_param == 128) {
    Keccak_HashInitialize_SHAKE128(ctx);
  } else {
    Keccak_HashInitialize_SHAKE256(ctx);
  }
}

static inline void hash_update(hash_context* ctx, const uint8_t* data, size_t size) {
  Keccak_HashUpdate(ctx, data, size << 3);
}

static inline void hash_final(hash_context* ctx) {
  Keccak_HashFinal(ctx, NULL);
}

static inline void hash_squeeze(hash_context* ctx, uint8_t* buffer, size_t buflen) {
  Keccak_HashSqueeze(ctx, buffer, buflen << 3);
}

static inline void hash_update_uint16_le(hash_context* ctx, uint16_t data) {
  const uint16_t data_le = htole16(data);
  hash_update(ctx, (const uint8_t*)&data_le, sizeof(data_le));
}

static inline void hash_init_prefix(hash_context* ctx, unsigned int security_param,
                                    const uint8_t prefix) {
  hash_init(ctx, security_param);
  hash_update(ctx, &prefix, sizeof(prefix));
}

#define hash_clear(ctx)
#endif

#if defined(WITH_KECCAK_X4)
/* Instances that work with 4 states in parallel. */
typedef Keccak_HashInstancetimes4 hash_context_x4;

static inline void hash_init_x4(hash_context_x4* ctx, unsigned int security_param) {
  if (security_param == 128) {
    Keccak_HashInitializetimes4_SHAKE128(ctx);
  } else {
    Keccak_HashInitializetimes4_SHAKE256(ctx);
  }
}

static inline void hash_update_x4(hash_context_x4* ctx, const uint8_t** data, size_t size) {
  Keccak_HashUpdatetimes4(ctx, data, size << 3);
}

static inline void hash_update_x4_4(hash_context_x4* ctx, const uint8_t* data0,
                                    const uint8_t* data1, const uint8_t* data2,
                                    const uint8_t* data3, size_t size) {
  const uint8_t* data[4] = {data0, data1, data2, data3};
  hash_update_x4(ctx, data, size);
}

static inline void hash_update_x4_1(hash_context_x4* ctx, const uint8_t* data, size_t size) {
  const uint8_t* tmp[4] = {data, data, data, data};
  hash_update_x4(ctx, tmp, size);
}

static inline void hash_init_prefix_x4(hash_context_x4* ctx, unsigned int security_param,
                                       const uint8_t prefix) {
  hash_init_x4(ctx, security_param);
  hash_update_x4_1(ctx, &prefix, sizeof(prefix));
}

static inline void hash_final_x4(hash_context_x4* ctx) {
  Keccak_HashFinaltimes4(ctx, NULL);
}

static inline void hash_squeeze_x4(hash_context_x4* ctx, uint8_t** buffer, size_t buflen) {
  Keccak_HashSqueezetimes4(ctx, buffer, buflen << 3);
}

static inline void hash_squeeze_x4_4(hash_context_x4* ctx, uint8_t* buffer0, uint8_t* buffer1,
                                     uint8_t* buffer2, uint8_t* buffer3, size_t buflen) {
  uint8_t* buffer[4] = {buffer0, buffer1, buffer2, buffer3};
  hash_squeeze_x4(ctx, buffer, buflen);
}

#define hash_clear_x4(ctx)
#elif !defined(OQS)
/* Instances that work with 4 states in parallel using the base Keccak implementation. */
typedef struct hash_context_x4_s {
  hash_context instances[4];
} hash_context_x4;

static inline void hash_init_x4(hash_context_x4* ctx, unsigned int security_param) {
  for (unsigned int i = 0; i < 4; ++i) {
    hash_init(&ctx->instances[i], security_param);
  }
}

static inline void hash_update_x4(hash_context_x4* ctx, const uint8_t** data, size_t size) {
  for (unsigned int i = 0; i < 4; ++i) {
    hash_update(&ctx->instances[i], data[i], size);
  }
}

static inline void hash_update_x4_4(hash_context_x4* ctx, const uint8_t* data0,
                                    const uint8_t* data1, const uint8_t* data2,
                                    const uint8_t* data3, size_t size) {
  hash_update(&ctx->instances[0], data0, size);
  hash_update(&ctx->instances[1], data1, size);
  hash_update(&ctx->instances[2], data2, size);
  hash_update(&ctx->instances[3], data3, size);
}

static inline void hash_update_x4_1(hash_context_x4* ctx, const uint8_t* data, size_t size) {
  for (unsigned int i = 0; i < 4; ++i) {
    hash_update(&ctx->instances[i], data, size);
  }
}

static inline void hash_init_prefix_x4(hash_context_x4* ctx, unsigned int security_param,
                                       const uint8_t prefix) {
  for (unsigned int i = 0; i < 4; ++i) {
    hash_init_prefix(&ctx->instances[i], security_param, prefix);
  }
}

static inline void hash_final_x4(hash_context_x4* ctx) {
  for (unsigned int i = 0; i < 4; ++i) {
    hash_final(&ctx->instances[i]);
  }
}

static inline void hash_squeeze_x4(hash_context_x4* ctx, uint8_t** buffer, size_t buflen) {
  for (unsigned int i = 0; i < 4; ++i) {
    hash_squeeze(&ctx->instances[i], buffer[i], buflen);
  }
}

static inline void hash_squeeze_x4_4(hash_context_x4* ctx, uint8_t* buffer0, uint8_t* buffer1,
                                     uint8_t* buffer2, uint8_t* buffer3, size_t buflen) {
  hash_squeeze(&ctx->instances[0], buffer0, buflen);
  hash_squeeze(&ctx->instances[1], buffer1, buflen);
  hash_squeeze(&ctx->instances[2], buffer2, buflen);
  hash_squeeze(&ctx->instances[3], buffer3, buflen);
}

#define hash_clear_x4(ctx)
#endif

static inline void hash_update_x4_uint16_le(hash_context_x4* ctx, uint16_t data) {
  const uint16_t data_le = htole16(data);
  hash_update_x4_1(ctx, (const uint8_t*)&data_le, sizeof(data_le));
}

static inline void hash_update_x4_uint16s_le(hash_context_x4* ctx, const uint16_t data[4]) {
  const uint16_t data0_le = htole16(data[0]);
  const uint16_t data1_le = htole16(data[1]);
  const uint16_t data2_le = htole16(data[2]);
  const uint16_t data3_le = htole16(data[3]);
  hash_update_x4_4(ctx, (const uint8_t*)&data0_le, (const uint8_t*)&data1_le,
                   (const uint8_t*)&data2_le, (const uint8_t*)&data3_le, sizeof(data[0]));
}

#endif

#endif /* HASH_SHAKE_H */