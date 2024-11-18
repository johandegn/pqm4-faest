/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest.h"
#include "faest_aes.h"
#include "fields.h"
#include "vole.h"
#include "universal_hashing.h"
#include "utils.h"
#include "parameters.h"
#include <stdio.h>

#include <string.h>
#include <stdlib.h>

static_assert(FAEST_128F_L == FAEST_128S_L, "Invalid parameters");
static_assert(FAEST_128F_LAMBDA == FAEST_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_128F_Lke == FAEST_128S_Lke, "Invalid parameters");
static_assert(FAEST_128F_Nwd == FAEST_128S_Nwd, "Invalid parameters");
static_assert(FAEST_128F_R == FAEST_128S_R, "Invalid parameters");
static_assert(FAEST_128F_Senc == FAEST_128S_Senc, "Invalid parameters");
static_assert(FAEST_128F_Ske == FAEST_128S_Ske, "Invalid parameters");

static_assert(FAEST_192F_L == FAEST_192S_L, "Invalid parameters");
static_assert(FAEST_192F_LAMBDA == FAEST_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_192F_Lke == FAEST_192S_Lke, "Invalid parameters");
static_assert(FAEST_192F_Nwd == FAEST_192S_Nwd, "Invalid parameters");
static_assert(FAEST_192F_R == FAEST_192S_R, "Invalid parameters");
static_assert(FAEST_192F_Senc == FAEST_192S_Senc, "Invalid parameters");
static_assert(FAEST_192F_Ske == FAEST_192S_Ske, "Invalid parameters");

static_assert(FAEST_256F_L == FAEST_256S_L, "Invalid parameters");
static_assert(FAEST_256F_LAMBDA == FAEST_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_256F_Lke == FAEST_256S_Lke, "Invalid parameters");
static_assert(FAEST_256F_Nwd == FAEST_256S_Nwd, "Invalid parameters");
static_assert(FAEST_256F_R == FAEST_256S_R, "Invalid parameters");
static_assert(FAEST_256F_Senc == FAEST_256S_Senc, "Invalid parameters");
static_assert(FAEST_256F_Ske == FAEST_256S_Ske, "Invalid parameters");

static_assert(FAEST_EM_128F_LAMBDA == FAEST_EM_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_128F_Lenc == FAEST_EM_128S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_128F_Nwd == FAEST_EM_128S_Nwd, "Invalid parameters");
static_assert(FAEST_EM_128F_R == FAEST_EM_128S_R, "Invalid parameters");
static_assert(FAEST_EM_128F_Senc == FAEST_EM_128S_Senc, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_128F_Nwd * (FAEST_EM_128F_R + 1),
              "Invalid parameters");

static_assert(FAEST_EM_192F_LAMBDA == FAEST_EM_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_192F_Lenc == FAEST_EM_192S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_192F_Nwd == FAEST_EM_192S_Nwd, "Invalid parameters");
static_assert(FAEST_EM_192F_R == FAEST_EM_192S_R, "Invalid parameters");
static_assert(FAEST_EM_192F_Senc == FAEST_EM_192S_Senc, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_192F_Nwd * (FAEST_EM_192F_R + 1),
              "Invalid parameters");

static_assert(FAEST_EM_256F_LAMBDA == FAEST_EM_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_256F_Lenc == FAEST_EM_256S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_256F_Nwd == FAEST_EM_256S_Nwd, "Invalid parameters");
static_assert(FAEST_EM_256F_R == FAEST_EM_256S_R, "Invalid parameters");
static_assert(FAEST_EM_256F_Senc == FAEST_EM_256S_Senc, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_256F_Nwd * (FAEST_EM_256F_R + 1),
              "Invalid parameters");

static const bf8_t Rcon[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf128_t* new_v = faest_aligned_alloc(BF128_ALIGN, (ell + FAEST_128F_LAMBDA) * sizeof(bf128_t));
  for (unsigned int row = 0; row != ell + FAEST_128F_LAMBDA; ++row) {
    uint8_t new_row[BF128_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_128F_LAMBDA; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf128_load(new_row);
  }

  return new_v;
}

bf192_t* column_to_row_major_and_shrink_V_192(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf192_t* new_v = faest_aligned_alloc(BF192_ALIGN, (ell + FAEST_192F_LAMBDA) * sizeof(bf192_t));
  for (unsigned int row = 0; row != ell + FAEST_192F_LAMBDA; ++row) {
    uint8_t new_row[BF192_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_192F_LAMBDA; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf192_load(new_row);
  }

  return new_v;
}

bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf256_t* new_v = faest_aligned_alloc(BF256_ALIGN, (ell + FAEST_256F_LAMBDA) * sizeof(bf256_t));
  for (unsigned int row = 0; row != ell + FAEST_256F_LAMBDA; ++row) {
    uint8_t new_row[BF256_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_256F_LAMBDA; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf256_load(new_row);
  }

  return new_v;
}

// m == 1 implementations

static void aes_key_schedule_forward_1(const uint8_t* x, uint8_t* out,
                                       const faest_paramset_t* params) {
  // Step: 1 skipped (sanity check)

  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int R           = params->faest_param.R;
  const unsigned int Nwd         = params->faest_param.Nwd;
  const unsigned int lambdaBytes = lambda / 8;

  const unsigned int out_len = (R + 1) * 128 / 8;
  // Step 3
  memcpy(out, x, lambdaBytes);
  memset(out + lambdaBytes, 0, out_len - lambdaBytes);

  // Step: 4
  unsigned int i_wd = lambda;
  // Step: 5..10
  for (unsigned int j = Nwd; j < 4 * (R + 1); j++) {
    if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
      memcpy(out + 32 * j / 8, x + i_wd / 8, 4);
      i_wd += 32;
    } else {
      for (unsigned int i = 0; i < 32; i += 8) {
        // bit spliced
        out[(32 * j + i) / 8] |= out[(32 * (j - Nwd) + i) / 8] ^ out[(32 * (j - 1) + i) / 8];
      }
    }
  }
}

static void aes_key_schedule_backward_1_round(const uint8_t* x, const uint8_t* xk, uint8_t* out,
                                              unsigned int j, const faest_paramset_t* params) {
  // Step: 1 skipped (sanity check)
  unsigned int lambda = params->faest_param.lambda;

  // Step: 2
  unsigned int iwd;
  if (lambda == 192){
    iwd   = 192/8 * j/4;
  } else {
    iwd   = 128/8 * j/4;
  }
  bool rmvRcon = true;
  unsigned int ircon = j/4;
  if (lambda == 256){
    rmvRcon       = ((j/4) %2) == 0;
    ircon = j/8;
  }

  for (unsigned int c = 0; c < 4; c++) {
    // Step 7 (bit sliced)
    uint8_t x_tilde = x[c + j] ^ xk[iwd + c];

    // Step 8
    // this function is only called with Mtag == Mkey == 0
    if (/* Mtag == 0 && */ rmvRcon == true && c == 0) {
      // Steps 12 and 13, bitsliced; delta is always 0
      x_tilde ^= Rcon[ircon];
      ++ircon;
    }

    // Step: 15..19 (bit spliced)
    const uint8_t y_tilde = rotr8(x_tilde, 7) ^ rotr8(x_tilde, 5) ^ rotr8(x_tilde, 2);
    // this function is only called with Mtag == Mkey == 0
    // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
    out[c] = y_tilde ^ 0x5;

    // Step: 20
    // TODO: fix this for all varaints
    /*
    ++c;
    if (c == 4) {
      c = 0;
      if (lambda == 192) {
        iwd += 192 / 8;
      } else {
        iwd += 128 / 8;
        if (lambda == 256) {
          rmvRcon = !rmvRcon;
        }
      }
    }
    */
  }
}

// lambda == 128 implementation
static void aes_key_schedule_backward_128_vbb_vk_round(vbb_t* vbb, uint8_t Mtag, uint8_t Mkey,
                                                 const uint8_t* delta, bf128_t* bf_out, unsigned int j) {
  // Step: 1
  assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();

  unsigned int iwd   = 128 * j/4;
  unsigned int ircon = j/4;

  bf128_t bf_minus_mkey       = bf128_from_bit(1 ^ Mkey);
  uint8_t minus_mtag          = 1 ^ Mtag;
  bf128_t bf_mkey_times_delta = bf128_mul_bit(bf_delta, Mkey);
  bf_mkey_times_delta         = bf128_add(bf_mkey_times_delta, bf_minus_mkey);

  // Step 7
  for(unsigned int c = 0; c < 4; c++) {
    bf128_t bf_x_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      bf128_t VOLE = *get_vole_aes_128(vbb, (8 * (j+c) + i) + FAEST_128F_LAMBDA);
      add_vole_to_vk_cache(vbb, (8 * (j+c) + i) + FAEST_128F_LAMBDA, &VOLE);
      bf_x_tilde[i] = VOLE;
    }
    for(unsigned int i = 0; i < 8; i++){
      bf_x_tilde[i] = bf128_add(bf_x_tilde[i], *get_vk_128(vbb, iwd + 8 * c + i));
    }

    if (Mtag == 0 && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf128_t bf_r[8];
      for (unsigned int i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf128_mul_bit(bf_mkey_times_delta, get_bit(r, i));
        // Step 13
        bf_x_tilde[i] = bf128_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    for (unsigned int i = 0; i < 8; ++i) {
      bf_out[i + 8 * c] = bf128_add(bf128_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
    }
    bf_out[0 + 8 * c] =
        bf128_add(bf_out[0 + 8 * c], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
    bf_out[2 + 8 * c] =
        bf128_add(bf_out[2 + 8 * c], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
  }
}

static void aes_key_schedule_128(const uint8_t* w, vbb_t* vbb,
                                                    zk_hash_128_ctx* a0_ctx,
                                                    zk_hash_128_ctx* a1_ctx, uint8_t* k,
                                                    const faest_paramset_t* params) {
  uint8_t w_dash[4];
  bf128_t v_w_dash[4 * 8];
  
  aes_key_schedule_forward_1(w, k, params);

  const unsigned int Ske    = params->faest_param.Ske;

  // backward 1
  unsigned int iwd   = 0;
  for (unsigned int j = 0; j < Ske; j+=4) {
    // backward 1
    aes_key_schedule_backward_1_round(w + FAEST_128F_LAMBDA / 8, k, w_dash, j, params);
    aes_key_schedule_backward_128_vbb_vk_round(vbb, 1, 0, NULL, v_w_dash, j);

    bf128_t bf_k_hat[4];
    bf128_t bf_v_k_hat[4];
    bf128_t bf_w_dash_hat[4];
    bf128_t bf_v_w_dash_hat[4];

    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 10..11
      bf_k_hat[(r + 3) % 4]   = bf128_byte_combine_bits(k[(96 + iwd*8 + 8 * r) / 8]);
      bf_v_k_hat[(r + 3) % 4] = bf128_byte_combine_vk(vbb, (96 + iwd*8 + 8 * r));
      bf_w_dash_hat[r]        = bf128_byte_combine_bits(w_dash[(8 * r) / 8]);
      bf_v_w_dash_hat[r]      = bf128_byte_combine(v_w_dash + (8 * r));
    }
    for (unsigned int r = 0; r <= 3; r++) {
      const bf128_t tmp = bf128_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
      zk_hash_128_update(a0_ctx, tmp);
      zk_hash_128_update(
          a1_ctx, bf128_add(bf128_add(bf128_mul(bf128_add(bf_k_hat[r], bf_v_k_hat[r]),
                                                bf128_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                      bf128_one()),
                            tmp));
    }
    iwd += 128 / 8;
  }
}

static void aes_key_schedule_constraints_Mkey_1_128_new(vbb_t* vbb, const uint8_t* delta,
                                                    zk_hash_128_ctx* b0_ctx) {
  // Step: 19..20
  // aes_key_schedule_forward_128_vbb(vbb, qk);
  bf128_t q_w_dash[4 * 8];

  const bf128_t bf_delta      = bf128_load(delta);
  const bf128_t delta_squared = bf128_mul(bf_delta, bf_delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_128F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_128F_Ske; j += 4) {
    aes_key_schedule_backward_128_vbb_vk_round(vbb, 0, 1, delta, q_w_dash, j);
    bf128_t bf_q_hat_k[4];
    bf128_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf128_byte_combine_vk(vbb, ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf128_byte_combine(q_w_dash + ((8 * r)));
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf128_t bf_tmp = bf128_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      // instead of storing B, hash it
      zk_hash_128_update(b0_ctx, bf128_add(bf_tmp, delta_squared));
    }
    iwd = iwd + 128;
  }
}

static void aes_enc_forward_128_1_round(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf128_t* bf_y, int round) {
  if (round == 0){
    for (unsigned int i = 0; i < 16; i++) {
      const uint8_t xin = in[i];
      bf_y[i] = bf128_add(bf128_byte_combine_bits(xin), bf128_byte_combine_bits(xk[i]));
    }
  }

  if (round > 0){
    const bf128_t bf_two   = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);
    unsigned int j = round;
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf128_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf128_byte_combine_bits(xk[(ik + 8 * r) / 8]);
      }

      // Step : 14
      bf_y[iy + 0] = bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
    }
  }
}

static void aes_enc_forward_backward_128(vbb_t* vbb, unsigned int offset, const uint8_t* in, const uint8_t* out,
                                         uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                         bf128_t* vs, bf128_t* vs_old, bf128_t* vs_dash, unsigned int round) {
  const bf128_t bf_delta  = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t bf_factor = bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey));
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  if (round == FAEST_128F_R -1){
    for (unsigned int c = 0; c <= 3; c++) {
        for (unsigned int r = 0; r <= 3; r++) {
          bf128_t bf_x_tilde[8];
          unsigned int ird = (128 * round) + (32 * ((c - r + 4) % 4)) + (8 * r);

          // Step: 10
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 11
            bf128_t bf_xout =
                bf128_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_128F_R - 1)) / 8], i));
            // Step: 12
            bf_x_tilde[i] = bf128_add(bf_xout, *get_vk_128(vbb, 128 + ird + i));
          }

          // Step: 13..17
          bf128_t bf_y_tilde[8];
          for (unsigned int i = 0; i < 8; ++i) {
            bf_y_tilde[i] = bf128_add(bf128_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                      bf_x_tilde[(i + 2) % 8]);
          }
          bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
          bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

          // Step: 18
          vs_dash[4 * c + r] = bf128_byte_combine(bf_y_tilde);
        }
      }
    return;
  }

  if (round == 0){
    // Step: 2..4
    for (unsigned int i = 0; i < 16; i++) {
      bf128_t bf_xin[8];
      for (unsigned int j = 0; j < 8; j++) {
        bf_xin[j] = bf128_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
      }
      // Step: 5
      vs[i] = bf128_add(bf128_byte_combine(bf_xin), bf128_byte_combine_vk(vbb, (8 * i)));
    }
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  unsigned int i_counter = 0;
  unsigned int j = round;
  for (unsigned int c = 0; c <= 3; c++) {
    const unsigned int ix = 128 * (j) + 32 * c;
    const unsigned int ik = 128 * (j+1) + 32 * c;
    const unsigned int iy = 4 * c;

    bf128_t bf_x_hat[4];
    bf128_t bf_xk_hat[4];
    for (unsigned int r = 0; r <= 3; r++) {
      bf128_t bf_x_tilde[8];
      for (unsigned int i = 0; i < 8; i++){
        bf128_t t = *get_vole_aes_128(vbb, offset + ix + 8 * r + i);
        memcpy(bf_x_tilde + i, &t, sizeof(bf128_t));
        i_counter++;

        // forward part
      }
      bf_x_hat[r]  =  bf128_byte_combine(bf_x_tilde);
      bf_xk_hat[r] = bf128_byte_combine_vk(vbb, (ik + 8 * r));

        // backwards part
      unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
      unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % 4;
      bf128_t bf_y_tilde[8];
      for (unsigned int k = 0; k < 8; ++k) {
        bf_y_tilde[k] = bf128_add(bf128_add(bf_x_tilde[(k + 7) % 8], bf_x_tilde[(k + 5) % 8]),
                                  bf_x_tilde[(k + 2) % 8]);
      }
      bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);
      vs_dash[4 * c_bkwd + r_bkwd] = bf128_byte_combine(bf_y_tilde);
    }

    // Step : 14
    vs_old[iy + 0] = bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
    vs_old[iy + 0] = bf128_add(vs_old[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
    vs_old[iy + 0] = bf128_add(vs_old[iy + 0], bf_x_hat[2]);
    vs_old[iy + 0] = bf128_add(vs_old[iy + 0], bf_x_hat[3]);

    // Step: 15
    vs_old[iy + 1] = bf128_add(bf_xk_hat[1], bf_x_hat[0]);
    vs_old[iy + 1] = bf128_add(vs_old[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
    vs_old[iy + 1] = bf128_add(vs_old[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
    vs_old[iy + 1] = bf128_add(vs_old[iy + 1], bf_x_hat[3]);

    // Step: 16
    vs_old[iy + 2] = bf128_add(bf_xk_hat[2], bf_x_hat[0]);
    vs_old[iy + 2] = bf128_add(vs_old[iy + 2], bf_x_hat[1]);
    vs_old[iy + 2] = bf128_add(vs_old[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
    vs_old[iy + 2] = bf128_add(vs_old[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

    // Step: 17
    vs_old[iy + 3] = bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
    vs_old[iy + 3] = bf128_add(vs_old[iy + 3], bf_x_hat[1]);
    vs_old[iy + 3] = bf128_add(vs_old[iy + 3], bf_x_hat[2]);
    vs_old[iy + 3] = bf128_add(vs_old[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
  }
}

static void aes_enc_backward_128_1_round(const uint8_t* x, const uint8_t* xk, const uint8_t* out,
                                   bf128_t* y_out, unsigned int round) {
  // called only with Mtag == Mkey == 0

  uint8_t xtilde;
  // Step:2..4
  unsigned int j = round;
  for (unsigned int c = 0; c <= 3; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 5..6
      unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
      if (j < (FAEST_128F_R - 1)) {
        // Step: 7
        xtilde = x[ird / 8];
      } else {
        // Step: 9..11 (bit spliced)
        // -((1 ^ Mtag) & (1 ^ Mkey)) == 0xff
        const uint8_t xout = out[(ird - 128 * (FAEST_128F_R - 1)) / 8];
        xtilde             = xout ^ xk[(128 + ird) / 8];
      }
      

      // Step: 12..17 (bit spliced)
      // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
      const uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2) ^ 0x5;

      // Step: 18
      y_out[4 * c + r] = bf128_byte_combine_bits(ytilde);
    }
  }
}

static void aes_constraints_0_128(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           vbb_t* vbb, unsigned int offset, const uint8_t* k,
                                           zk_hash_128_ctx* a0_ctx, zk_hash_128_ctx* a1_ctx) {
  unsigned int w_offset = offset / 8;
  w += w_offset;

  bf128_t s[16];
  bf128_t s_dash[16];
  bf128_t vs[16];
  bf128_t vs_old[16];
  bf128_t vs_dash[16];

  for (unsigned int j = 0; j < FAEST_128F_R; j++) {
    // forward 1
    aes_enc_forward_128_1_round(w, k, in, s, j);
    // backward 1
    aes_enc_backward_128_1_round(w, k, out, s_dash, j);

    if (j != 0){
      memcpy(vs, vs_old, sizeof(vs));
    }
    aes_enc_forward_backward_128(vbb, offset, in, out, 1, 0, NULL, vs, vs_old, vs_dash, j);
    for (unsigned int i = 0; i < 16; i++) {
        const bf128_t tmp = bf128_mul(vs[i], vs_dash[i]);
        zk_hash_128_update(a0_ctx, tmp);
        zk_hash_128_update(a1_ctx, bf128_add(bf128_add(bf128_mul(bf128_add(s[i], vs[i]),
                                                                bf128_add(s_dash[i], vs_dash[i])),
                                                      tmp),
                                            bf128_one()));
      }
  }
}

static void aes_enc_constraints_Mkey_1_128(const uint8_t* in, const uint8_t* out, vbb_t* vbb,
                                           unsigned int offset, const uint8_t* delta,
                                           zk_hash_128_ctx* b0_ctx) {

  // Step: 11..12
  bf128_t qs[16];
  bf128_t qs_old[16];
  bf128_t qs_dash[16];
  for (unsigned int j = 0; j < FAEST_128F_R; j++){
    /*
    aes_enc_forward_128_vbb_vk(vbb, offset, in, 0, 1, delta, qs);
    aes_enc_backward_128_vbb_linear_access(vbb, offset, 0, 1, delta, out, qs_dash);
    */
    if(j != 0){
      memcpy(qs, qs_old, sizeof(qs));
    }
    aes_enc_forward_backward_128(vbb, offset, in, out, 0, 1, delta, qs, qs_old, qs_dash, j);

    // Step: 13..14
    bf128_t minus_part = bf128_mul(bf128_load(delta), bf128_load(delta));
    for (unsigned int i = 0; i < 16; i++){
      // instead of storing it, hash it
      zk_hash_128_update(b0_ctx, bf128_add(bf128_mul(qs[i], qs_dash[i]), minus_part));
    }
  }
}

static void aes_prove_128(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                          const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
                          const faest_paramset_t* params) {
  // Step: 1..2
  // compute on the fly

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k = alloca((FAEST_128F_R + 1) * 128 / 8);
  // bf128_t* vk = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((FAEST_128F_R + 1) * 128));
  zk_hash_128_ctx a0_ctx;
  zk_hash_128_ctx a1_ctx;

  zk_hash_128_init(&a0_ctx, chall);
  zk_hash_128_init(&a1_ctx, chall);
  aes_key_schedule_128(w, vbb, &a0_ctx, &a1_ctx, k, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  unsigned int offset = FAEST_128F_Lke;
  aes_constraints_0_128(in, out, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // Step: 12 (beta == 1)
  // faest_aligned_free(vk);
  //free(k);

  // Step: 16..18
  zk_hash_128_finalize(a_tilde, &a1_ctx, bf128_load(get_vole_u(vbb) + FAEST_128F_L / 8));
  zk_hash_128_finalize(b_tilde, &a0_ctx, bf128_sum_poly_vbb(vbb, FAEST_128F_L));
}

static uint8_t* aes_verify_128(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                               const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out, uint8_t* q_tilde) {
  // Step: 1
  const uint8_t* delta = chall_3;

  // Step: 13 + 21
  // bf128_t* qk = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((FAEST_128F_R + 1) * 128));
  // instead of storing B_0 in an array, we process the values with zk_hash_128
  zk_hash_128_ctx b0_ctx;
  zk_hash_128_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_128_new(vbb, delta, &b0_ctx);

  // Step: 14
  aes_enc_constraints_Mkey_1_128(in, out, vbb, FAEST_128F_Lke, delta, &b0_ctx);
  // Step: 18 (beta == 1)
  // faest_aligned_free(qk);

  // Step: 20+21
  //uint8_t* q_tilde = malloc(FAEST_128F_LAMBDA / 8);
  zk_hash_128_finalize(q_tilde, &b0_ctx, bf128_sum_poly_vbb(vbb, FAEST_128F_L));

  bf128_t bf_qtilde = bf128_load(q_tilde);
  bf128_store(q_tilde, bf128_add(bf_qtilde, bf128_mul(bf128_load(a_tilde), bf128_load(delta))));

  return q_tilde;
}

// lambda == 192 implementation

static void aes_key_schedule_backward_192_vbb_vk_round(vbb_t* vbb, uint8_t Mtag, uint8_t Mkey,
                                                 const uint8_t* delta, bf192_t* bf_out, unsigned int j) {
  // Step: 1
  assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  unsigned int iwd       = 192 * j/4;
  unsigned int ircon     = j/4;

  bf192_t bf_minus_mkey       = bf192_from_bit(1 ^ Mkey);
  uint8_t minus_mtag          = 1 ^ Mtag;
  bf192_t bf_mkey_times_delta = bf192_mul_bit(bf_delta, Mkey);
  bf_mkey_times_delta         = bf192_add(bf_mkey_times_delta, bf_minus_mkey);

  for (unsigned int c = 0; c < 4; c++) {
    // Step 7
    bf192_t bf_x_tilde[8];
    for (unsigned int i = 0; i < 8; i++){
      bf192_t VOLE = *get_vole_aes_192(vbb, (8 * (j+c) + i) + FAEST_192F_LAMBDA);
      add_vole_to_vk_cache_192(vbb, (8 * (j+c) + i) + FAEST_192F_LAMBDA, &VOLE);
      bf_x_tilde[i] = VOLE;
    }

    for (unsigned int i = 0; i < 8; i++) {
      bf_x_tilde[i] = bf192_add(bf_x_tilde[i], *get_vk_192(vbb, iwd + 8 * c + i));
    }

    if (Mtag == 0 && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf192_t bf_r[8];
      for (unsigned int i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf192_mul_bit(bf_mkey_times_delta, get_bit(r, i));
        // Step 13
        bf_x_tilde[i] = bf192_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    for (unsigned int i = 0; i < 8; ++i) {
      bf_out[i + 8 * c] = bf192_add(bf192_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
    }
    bf_out[0 + 8 * c] =
        bf192_add(bf_out[0 + 8 * c], bf192_mul_bit(bf_mkey_times_delta, minus_mtag));
    bf_out[2 + 8 * c] =
        bf192_add(bf_out[2 + 8 * c], bf192_mul_bit(bf_mkey_times_delta, minus_mtag));
  }
}

static void aes_key_schedule_192(const uint8_t* w, vbb_t* vbb,
                                                    zk_hash_192_ctx* a0_ctx,
                                                    zk_hash_192_ctx* a1_ctx, uint8_t* k,
                                                    const faest_paramset_t* params) {
  uint8_t w_dash[4];
  bf192_t v_w_dash[4 * 8];

  aes_key_schedule_forward_1(w, k, params);

  unsigned int iwd = 32 * (FAEST_192F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_192F_Ske; j+= 4) {
    aes_key_schedule_backward_1_round(w + FAEST_192F_LAMBDA / 8, k, w_dash, j, params);
    aes_key_schedule_backward_192_vbb_vk_round(vbb, 1, 0, NULL, v_w_dash, j);

    bf192_t bf_k_hat[4];
    bf192_t bf_v_k_hat[4];
    bf192_t bf_w_dash_hat[4];
    bf192_t bf_v_w_dash_hat[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 10..11
      bf_k_hat[(r + 3) % 4]   = bf192_byte_combine_bits(k[(iwd + 8 * r) / 8]);
      bf_v_k_hat[(r + 3) % 4] = bf192_byte_combine_vk(vbb, (iwd + 8 * r));
      bf_w_dash_hat[r]        = bf192_byte_combine_bits(w_dash[(8 * r) / 8]);
      bf_v_w_dash_hat[r]      = bf192_byte_combine(v_w_dash + (8 * r));
    }
    // Step: 13..17
    for (unsigned int r = 0; r <= 3; r++) {
      // instead of storing in A0, A1, hash it
      const bf192_t tmp = bf192_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
      zk_hash_192_update(a0_ctx, tmp);
      zk_hash_192_update(
          a1_ctx, bf192_add(bf192_add(bf192_mul(bf192_add(bf_k_hat[r], bf_v_k_hat[r]),
                                                bf192_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                      bf192_one()),
                            tmp));
    }
    iwd = iwd + 192;
  }
}

static void aes_key_schedule_constraints_Mkey_1_192_new(vbb_t* vbb, const uint8_t* delta,
                                                    zk_hash_192_ctx* b0_ctx) {
  // Step: 19..20
  // aes_key_schedule_forward_192_vbb(vbb, qk);
  bf192_t q_w_dash[4 * 8];

  const bf192_t bf_delta      = bf192_load(delta);
  const bf192_t delta_squared = bf192_mul(bf_delta, bf_delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_192F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_192F_Ske; j += 4) {
    aes_key_schedule_backward_192_vbb_vk_round(vbb, 0, 1, delta, q_w_dash, j);
    bf192_t bf_q_hat_k[4];
    bf192_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf192_byte_combine_vk(vbb, ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf192_byte_combine(q_w_dash + ((8 * r)));
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf192_t bf_tmp = bf192_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      // instead of storing B, hash it
      zk_hash_192_update(b0_ctx, bf192_add(bf_tmp, delta_squared));
    }
    iwd = iwd + 192;
  }
}

static void aes_enc_forward_192_1_round(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  uint8_t Mtag, uint8_t Mkey, bf192_t* bf_y, unsigned int round) {
  // Step: 2
  if (round == 0) {
    for (unsigned int i = 0; i < 16; i++) {
      // Step: 3,4 (bit spliced)
      const uint8_t xin = in[i] & -((1 ^ Mtag) & (1 ^ Mkey));
      // Step: 5
      bf_y[i] = bf192_add(bf192_byte_combine_bits(xin), bf192_byte_combine_bits(xk[i]));
    }
  } else {
    const bf192_t bf_two   = bf192_byte_combine_bits(2);
    const bf192_t bf_three = bf192_byte_combine_bits(3);
    unsigned int j = round;
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf192_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf192_byte_combine_bits(xk[(ik + 8 * r) / 8]);
      }

      // Step : 14
      bf_y[iy + 0] = bf192_add(bf_xk_hat[0], bf192_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf192_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf192_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf192_add(bf_xk_hat[3], bf192_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[3], bf_two));
    }
  }
}

static void aes_enc_backward_192_1_round(const uint8_t* x, const uint8_t* xk, uint8_t Mtag, uint8_t Mkey,
                                   const uint8_t* out, bf192_t* y_out, unsigned int round) {
  uint8_t xtilde;
  // Step:2..4
  unsigned int j = round;
  for (unsigned int c = 0; c <= 3; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 5..6
      unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
      if (j < (FAEST_192F_R - 1)) {
        // Step: 7
        xtilde = x[ird / 8];
      } else {
        // Step: 9..11 (bit spliced)
        uint8_t xout = out[(ird - 128 * (FAEST_192F_R - 1)) / 8] & -((1 ^ Mtag) & (1 ^ Mkey));
        xtilde       = xout ^ xk[(128 + ird) / 8];
      }

      // Step: 12..17 (bit spliced)
      uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2);
      ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 0);
      ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 2);

      // Step: 18
      y_out[4 * c + r] = bf192_byte_combine_bits(ytilde);
    }
  }
}

static void aes_enc_forward_backward_192(vbb_t* vbb, unsigned int offset, const uint8_t* in, const uint8_t* out,
                                         uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                         bf192_t* vs, bf192_t* vs_old, bf192_t* vs_dash, unsigned int round) {
  const bf192_t bf_delta  = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t bf_factor = bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey));
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  if (round == FAEST_192F_R -1){
    for (unsigned int c = 0; c <= 3; c++) {
        for (unsigned int r = 0; r <= 3; r++) {
          bf192_t bf_x_tilde[8];
          unsigned int ird = (128 * round) + (32 * ((c - r + 4) % 4)) + (8 * r);

          // Step: 10
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 11
            bf192_t bf_xout =
                bf192_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_192F_R - 1)) / 8], i));
            // Step: 12
            bf_x_tilde[i] = bf192_add(bf_xout, *get_vk_192(vbb, 128 + ird + i));
          }

          // Step: 13..17
          bf192_t bf_y_tilde[8];
          for (unsigned int i = 0; i < 8; ++i) {
            bf_y_tilde[i] = bf192_add(bf192_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                      bf_x_tilde[(i + 2) % 8]);
          }
          bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
          bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

          // Step: 18
          vs_dash[4 * c + r] = bf192_byte_combine(bf_y_tilde);
        }
      }
    return;
  }

  if (round == 0){
    // Step: 2..4
    for (unsigned int i = 0; i < 16; i++) {
      bf192_t bf_xin[8];
      for (unsigned int j = 0; j < 8; j++) {
        bf_xin[j] = bf192_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
      }
      // Step: 5
      vs[i] = bf192_add(bf192_byte_combine(bf_xin), bf192_byte_combine_vk(vbb, (8 * i)));
    }
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  unsigned int i_counter = 0;
  unsigned int j = round;
  for (unsigned int c = 0; c <= 3; c++) {
    const unsigned int ix = 128 * (j) + 32 * c;
    const unsigned int ik = 128 * (j+1) + 32 * c;
    const unsigned int iy = 4 * c;

    bf192_t bf_x_hat[4];
    bf192_t bf_xk_hat[4];
    for (unsigned int r = 0; r <= 3; r++) {
      bf192_t bf_x_tilde[8];
      for (unsigned int i = 0; i < 8; i++){
        bf192_t t = *get_vole_aes_192(vbb, offset + ix + 8 * r + i);
        memcpy(bf_x_tilde + i, &t, sizeof(bf192_t));
        i_counter++;

        // forward part
      }
      bf_x_hat[r]  =  bf192_byte_combine(bf_x_tilde);
      bf_xk_hat[r] = bf192_byte_combine_vk(vbb, (ik + 8 * r));

        // backwards part
      unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
      unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % 4;
      bf192_t bf_y_tilde[8];
      for (unsigned int k = 0; k < 8; ++k) {
        bf_y_tilde[k] = bf192_add(bf192_add(bf_x_tilde[(k + 7) % 8], bf_x_tilde[(k + 5) % 8]),
                                  bf_x_tilde[(k + 2) % 8]);
      }
      bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);
      vs_dash[4 * c_bkwd + r_bkwd] = bf192_byte_combine(bf_y_tilde);
    }

    // Step : 14
    vs_old[iy + 0] = bf192_add(bf_xk_hat[0], bf192_mul(bf_x_hat[0], bf_two));
    vs_old[iy + 0] = bf192_add(vs_old[iy + 0], bf192_mul(bf_x_hat[1], bf_three));
    vs_old[iy + 0] = bf192_add(vs_old[iy + 0], bf_x_hat[2]);
    vs_old[iy + 0] = bf192_add(vs_old[iy + 0], bf_x_hat[3]);

    // Step: 15
    vs_old[iy + 1] = bf192_add(bf_xk_hat[1], bf_x_hat[0]);
    vs_old[iy + 1] = bf192_add(vs_old[iy + 1], bf192_mul(bf_x_hat[1], bf_two));
    vs_old[iy + 1] = bf192_add(vs_old[iy + 1], bf192_mul(bf_x_hat[2], bf_three));
    vs_old[iy + 1] = bf192_add(vs_old[iy + 1], bf_x_hat[3]);

    // Step: 16
    vs_old[iy + 2] = bf192_add(bf_xk_hat[2], bf_x_hat[0]);
    vs_old[iy + 2] = bf192_add(vs_old[iy + 2], bf_x_hat[1]);
    vs_old[iy + 2] = bf192_add(vs_old[iy + 2], bf192_mul(bf_x_hat[2], bf_two));
    vs_old[iy + 2] = bf192_add(vs_old[iy + 2], bf192_mul(bf_x_hat[3], bf_three));

    // Step: 17
    vs_old[iy + 3] = bf192_add(bf_xk_hat[3], bf192_mul(bf_x_hat[0], bf_three));
    vs_old[iy + 3] = bf192_add(vs_old[iy + 3], bf_x_hat[1]);
    vs_old[iy + 3] = bf192_add(vs_old[iy + 3], bf_x_hat[2]);
    vs_old[iy + 3] = bf192_add(vs_old[iy + 3], bf192_mul(bf_x_hat[3], bf_two));
  }
}

static void aes_enc_constraints_Mkey_0_192(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           vbb_t* vbb, unsigned int offset, const uint8_t* k,
                                           zk_hash_192_ctx* a0_ctx, zk_hash_192_ctx* a1_ctx) {
  unsigned int w_offset = offset / 8;
  w += w_offset;
  bf192_t s[16];
  bf192_t s_dash[16];
  bf192_t vs[16];
  bf192_t vs_old[16];
  bf192_t vs_dash[16];
  for (unsigned int j = 0; j < FAEST_192F_R; j++) {
    aes_enc_forward_192_1_round(w, k, in, 0, 0, s, j);
    aes_enc_backward_192_1_round(w, k, 0, 0, out, s_dash, j);
    if (j != 0){
      memcpy(vs, vs_old, sizeof(vs));
    }
    aes_enc_forward_backward_192(vbb, offset, in, out, 1, 0, NULL, vs, vs_old, vs_dash, j);
    for(unsigned int i = 0; i < 16; i++){
      const bf192_t tmp = bf192_mul(vs[i], vs_dash[i]);
      zk_hash_192_update(a0_ctx, tmp);
      zk_hash_192_update(a1_ctx, bf192_add(bf192_add(bf192_mul(bf192_add(s[i], vs[i]),
                                                              bf192_add(s_dash[i], vs_dash[i])),
                                                    tmp),
                                          bf192_one()));
    }
  }
}

static void aes_enc_constraints_Mkey_1_192(const uint8_t* in, const uint8_t* out, vbb_t* vbb,
                                           unsigned int offset, const uint8_t* delta,
                                           zk_hash_192_ctx* b0_ctx) {
  // Step: 11..12
  bf192_t qs[16];
  bf192_t qs_old[16];
  bf192_t qs_dash[16];

  // Step: 13..14
  for (unsigned int j = 0; j < FAEST_192F_R; j++) {
    if(j != 0){
      memcpy(qs, qs_old, sizeof(qs));
    }
    aes_enc_forward_backward_192(vbb, offset, in, out, 0, 1, delta, qs, qs_old, qs_dash, j);

    bf192_t minus_part = bf192_mul(bf192_load(delta), bf192_load(delta));
    for (unsigned int i = 0; i < 16; i++){
      zk_hash_192_update(b0_ctx, bf192_add(bf192_mul(qs[i], qs_dash[i]), minus_part));
    }
  }
}

static void aes_prove_192(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                          const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
                          const faest_paramset_t* params) {
  // Step: 1..2
  // compute on the fly

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k = alloca((FAEST_192F_R + 1) * 128 / 8);
  // bf192_t* vk = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  zk_hash_192_ctx a0_ctx;
  zk_hash_192_ctx a1_ctx;

  zk_hash_192_init(&a0_ctx, chall);
  zk_hash_192_init(&a1_ctx, chall);
  aes_key_schedule_192(w, vbb, &a0_ctx, &a1_ctx, k, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  unsigned int offset = FAEST_192F_Lke;
  aes_enc_constraints_Mkey_0_192(in, out, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // Step: 12-15
  offset = FAEST_192F_Lke + FAEST_192F_Lenc;
  aes_enc_constraints_Mkey_0_192(in + 16, out + 16, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // faest_aligned_free(vk);
  //free(k);

  // Step: 16..18
  zk_hash_192_finalize(a_tilde, &a1_ctx, bf192_load(get_vole_u(vbb) + FAEST_192F_L / 8));
  zk_hash_192_finalize(b_tilde, &a0_ctx, bf192_sum_poly_vbb(vbb, FAEST_192F_L));
}

static uint8_t* aes_verify_192(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                               const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out, uint8_t* q_tilde) {
  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 13 + 21
  // bf192_t* qk = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  // instead of storing B_0 in an array, we process the values with zk_hash_128
  zk_hash_192_ctx b0_ctx;
  zk_hash_192_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_192_new(vbb, delta, &b0_ctx);

  // Step: 14
  aes_enc_constraints_Mkey_1_192(in, out, vbb, FAEST_192F_Lke, delta, &b0_ctx);

  // Step: 18
  aes_enc_constraints_Mkey_1_192(in + 16, out + 16, vbb, FAEST_192F_Lke + FAEST_192F_Lenc, delta,
                                 &b0_ctx);
  // faest_aligned_free(qk);

  // Step: 20+21
  //uint8_t* q_tilde = malloc(FAEST_192F_LAMBDA / 8);
  zk_hash_192_finalize(q_tilde, &b0_ctx, bf192_sum_poly_vbb(vbb, FAEST_192F_L));

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

// lambda == 256 implementation
static void aes_key_schedule_backward_256_vbb_vk_round(vbb_t* vbb, uint8_t Mtag, uint8_t Mkey,
                                                 const uint8_t* delta, bf256_t* bf_out, unsigned int j) {
  // Step: 1
  assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

  unsigned int iwd   = 128 * j/4;
  bool rmvRcon       = (j/4)%2 == 0;//true;
  unsigned int ircon = j/8;

  const bf256_t bf_delta      = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t bf_minus_mkey = bf256_from_bit(1 ^ Mkey);
  const uint8_t minus_mtag    = 1 ^ Mtag;
  bf256_t bf_mkey_times_delta = bf256_mul_bit(bf_delta, Mkey);
  bf_mkey_times_delta         = bf256_add(bf_mkey_times_delta, bf_minus_mkey);

  for(unsigned int c = 0; c < 4; c++) {
    bf256_t bf_x_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      bf256_t VOLE = *get_vole_aes_256(vbb, (8 * (j+c) + i) + FAEST_256F_LAMBDA);
      add_vole_to_vk_cache_256(vbb, (8 * (j+c) + i) + FAEST_256F_LAMBDA, &VOLE);
      bf_x_tilde[i] = VOLE;
    }
    for(unsigned int i = 0; i < 8; i++){
      bf_x_tilde[i] = bf256_add(bf_x_tilde[i], *get_vk_256(vbb, iwd + 8 * c + i));
    }

    if (Mtag == 0 && rmvRcon == true && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf256_t bf_r[8];
      for (unsigned int i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf256_mul_bit(bf_mkey_times_delta, get_bit(r, i));
        // Step 13
        bf_x_tilde[i] = bf256_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    for (unsigned int i = 0; i < 8; ++i) {
      bf_out[i + 8 * c] = bf256_add(bf256_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
    }
    bf_out[0 + 8 * c] =
        bf256_add(bf_out[0 + 8 * c], bf256_mul_bit(bf_mkey_times_delta, minus_mtag));
    bf_out[2 + 8 * c] =
        bf256_add(bf_out[2 + 8 * c], bf256_mul_bit(bf_mkey_times_delta, minus_mtag));
  }
}

static void aes_key_schedule_256(const uint8_t* w, vbb_t* vbb,
                                                    zk_hash_256_ctx* a0_ctx,
                                                    zk_hash_256_ctx* a1_ctx, uint8_t* k,
                                                    const faest_paramset_t* params) {
  aes_key_schedule_forward_1(w, k, params);
  uint8_t w_dash[4];
  bf256_t v_w_dash[4 * 8];

  // Step: 6..8
  bool rotate_word = true;
  unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_256F_Ske; j += 4) {
    aes_key_schedule_backward_1_round(w + FAEST_256F_LAMBDA / 8, k, w_dash, j, params);
    aes_key_schedule_backward_256_vbb_vk_round(vbb, 1, 0, NULL, v_w_dash, j);
    
    bf256_t bf_k_hat[4];
    bf256_t bf_v_k_hat[4];
    bf256_t bf_w_dash_hat[4];
    bf256_t bf_v_w_dash_hat[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 10..11
      if (rotate_word) {
        bf_k_hat[(r + 3) % 4]   = bf256_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[(r + 3) % 4] = bf256_byte_combine_vk(vbb, (iwd + 8 * r));
        bf_w_dash_hat[r]        = bf256_byte_combine_bits(w_dash[(8 * r) / 8]);
        bf_v_w_dash_hat[r]      = bf256_byte_combine(v_w_dash + (8 * r));
      } else {
        bf_k_hat[r]        = bf256_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[r]      = bf256_byte_combine_vk(vbb, (iwd + 8 * r));
        bf_w_dash_hat[r]   = bf256_byte_combine_bits(w_dash[(8 * r) / 8]);
        bf_v_w_dash_hat[r] = bf256_byte_combine(v_w_dash + (8 * r));
      }
    }
    // Step: 13..17
    for (unsigned int r = 0; r <= 3; r++) {
      const bf256_t tmp = bf256_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
      zk_hash_256_update(a0_ctx, tmp);
      zk_hash_256_update(
          a1_ctx, bf256_add(bf256_add(bf256_mul(bf256_add(bf_k_hat[r], bf_v_k_hat[r]),
                                                bf256_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                      bf256_one()),
                            tmp));
    }
    iwd         = iwd + 128;
    rotate_word = !rotate_word;
  }
}

static void aes_key_schedule_constraints_Mkey_1_256_new(vbb_t* vbb, const uint8_t* delta,
                                                    zk_hash_256_ctx* b0_ctx) {
  // Step: 19..20
  // aes_key_schedule_forward_256_vbb(vbb, qk);
  bf256_t q_w_dash[4 * 8];

  const bf256_t bf_delta         = bf256_load(delta);
  const bf256_t bf_delta_squared = bf256_mul(bf_delta, bf_delta);
  bool rotate_word               = true;

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_256F_Ske; j += 4) {
    aes_key_schedule_backward_256_vbb_vk_round(vbb, 0, 1, delta, q_w_dash, j);
    bf256_t bf_q_hat_k[4];
    bf256_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      if (rotate_word) {
        bf_q_hat_k[(r + 3) % 4] = bf256_byte_combine_vk(vbb, ((iwd + 8 * r)));
        bf_q_hat_w_dash[r]      = bf256_byte_combine(q_w_dash + ((8 * r)));
      } else {
        bf_q_hat_k[r]      = bf256_byte_combine_vk(vbb, ((iwd + 8 * r)));
        bf_q_hat_w_dash[r] = bf256_byte_combine(q_w_dash + ((8 * r)));
      }
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf256_t bf_tmp = bf256_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      zk_hash_256_update(b0_ctx, bf256_add(bf_tmp, bf_delta_squared));
    }
    iwd         = iwd + 128;
    rotate_word = !rotate_word;
  }
}

static void aes_enc_forward_256_1_round(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  uint8_t Mtag, uint8_t Mkey, bf256_t* bf_y, unsigned int round) {
  // Step: 2
  if(round == 0){
    for (unsigned int i = 0; i < 16; i++) {
      // Step: 3,4 (bit spliced)
      const uint8_t xin = in[i] & -((1 ^ Mtag) & (1 ^ Mkey));
      // Step: 5
      bf_y[i] = bf256_add(bf256_byte_combine_bits(xin), bf256_byte_combine_bits(xk[i]));
    }
  } else {
    const bf256_t bf_two   = bf256_byte_combine_bits(2);
    const bf256_t bf_three = bf256_byte_combine_bits(3);
    unsigned int j = round;
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf256_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf256_byte_combine_bits(xk[(ik + 8 * r) / 8]);
      }

      // Step : 14
      bf_y[iy + 0] = bf256_add(bf_xk_hat[0], bf256_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf256_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf256_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf256_add(bf_xk_hat[3], bf256_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[3], bf_two));
    }
  }
}

static void aes_enc_backward_256_1_round(const uint8_t* x, const uint8_t* xk, uint8_t Mtag, uint8_t Mkey,
                                   const uint8_t* out, bf256_t* y_out, unsigned int round) {
  uint8_t xtilde;
  // Step:2..4
  unsigned int j = round;
  for (unsigned int c = 0; c <= 3; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 5..6
      unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
      if (j < (FAEST_256F_R - 1)) {
        // Step: 7
        xtilde = x[ird / 8];
      } else {
        // Step: 9..11 (bit spliced)
        uint8_t xout = out[(ird - 128 * (FAEST_256F_R - 1)) / 8] & -((1 ^ Mtag) & (1 ^ Mkey));
        xtilde       = xout ^ xk[(128 + ird) / 8];
      }

      // Step: 12..17 (bit spliced)
      uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2);
      ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 0);
      ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 2);

      // Step: 18
      y_out[4 * c + r] = bf256_byte_combine_bits(ytilde);
    }
  }
}

static void aes_enc_forward_backward_256(vbb_t* vbb, unsigned int offset, const uint8_t* in, const uint8_t* out,
                                         uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                         bf256_t* vs, bf256_t* vs_old, bf256_t* vs_dash, unsigned int round) {
  const bf256_t bf_delta  = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t bf_factor = bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey));
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  if (round == FAEST_256F_R -1){
    for (unsigned int c = 0; c <= 3; c++) {
        for (unsigned int r = 0; r <= 3; r++) {
          bf256_t bf_x_tilde[8];
          unsigned int ird = (128 * round) + (32 * ((c - r + 4) % 4)) + (8 * r);

          // Step: 10
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 11
            bf256_t bf_xout =
                bf256_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_256F_R - 1)) / 8], i));
            // Step: 12
            bf_x_tilde[i] = bf256_add(bf_xout, *get_vk_256(vbb, 128 + ird + i));
          }

          // Step: 13..17
          bf256_t bf_y_tilde[8];
          for (unsigned int i = 0; i < 8; ++i) {
            bf_y_tilde[i] = bf256_add(bf256_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                      bf_x_tilde[(i + 2) % 8]);
          }
          bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
          bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

          // Step: 18
          vs_dash[4 * c + r] = bf256_byte_combine(bf_y_tilde);
        }
      }
    return;
  }

  if (round == 0){
    // Step: 2..4
    for (unsigned int i = 0; i < 16; i++) {
      bf256_t bf_xin[8];
      for (unsigned int j = 0; j < 8; j++) {
        bf_xin[j] = bf256_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
      }
      // Step: 5
      vs[i] = bf256_add(bf256_byte_combine(bf_xin), bf256_byte_combine_vk(vbb, (8 * i)));
    }
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  unsigned int i_counter = 0;
  unsigned int j = round;
  for (unsigned int c = 0; c <= 3; c++) {
    const unsigned int ix = 128 * (j) + 32 * c;
    const unsigned int ik = 128 * (j+1) + 32 * c;
    const unsigned int iy = 4 * c;

    bf256_t bf_x_hat[4];
    bf256_t bf_xk_hat[4];
    for (unsigned int r = 0; r <= 3; r++) {
      bf256_t bf_x_tilde[8];
      for (unsigned int i = 0; i < 8; i++){
        bf256_t t = *get_vole_aes_256(vbb, offset + ix + 8 * r + i);
        memcpy(bf_x_tilde + i, &t, sizeof(bf256_t));
        i_counter++;

        // forward part
      }
      bf_x_hat[r]  =  bf256_byte_combine(bf_x_tilde);
      bf_xk_hat[r] = bf256_byte_combine_vk(vbb, (ik + 8 * r));

        // backwards part
      unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
      unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % 4;
      bf256_t bf_y_tilde[8];
      for (unsigned int k = 0; k < 8; ++k) {
        bf_y_tilde[k] = bf256_add(bf256_add(bf_x_tilde[(k + 7) % 8], bf_x_tilde[(k + 5) % 8]),
                                  bf_x_tilde[(k + 2) % 8]);
      }
      bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);
      vs_dash[4 * c_bkwd + r_bkwd] = bf256_byte_combine(bf_y_tilde);
    }

    // Step : 14
    vs_old[iy + 0] = bf256_add(bf_xk_hat[0], bf256_mul(bf_x_hat[0], bf_two));
    vs_old[iy + 0] = bf256_add(vs_old[iy + 0], bf256_mul(bf_x_hat[1], bf_three));
    vs_old[iy + 0] = bf256_add(vs_old[iy + 0], bf_x_hat[2]);
    vs_old[iy + 0] = bf256_add(vs_old[iy + 0], bf_x_hat[3]);

    // Step: 15
    vs_old[iy + 1] = bf256_add(bf_xk_hat[1], bf_x_hat[0]);
    vs_old[iy + 1] = bf256_add(vs_old[iy + 1], bf256_mul(bf_x_hat[1], bf_two));
    vs_old[iy + 1] = bf256_add(vs_old[iy + 1], bf256_mul(bf_x_hat[2], bf_three));
    vs_old[iy + 1] = bf256_add(vs_old[iy + 1], bf_x_hat[3]);

    // Step: 16
    vs_old[iy + 2] = bf256_add(bf_xk_hat[2], bf_x_hat[0]);
    vs_old[iy + 2] = bf256_add(vs_old[iy + 2], bf_x_hat[1]);
    vs_old[iy + 2] = bf256_add(vs_old[iy + 2], bf256_mul(bf_x_hat[2], bf_two));
    vs_old[iy + 2] = bf256_add(vs_old[iy + 2], bf256_mul(bf_x_hat[3], bf_three));

    // Step: 17
    vs_old[iy + 3] = bf256_add(bf_xk_hat[3], bf256_mul(bf_x_hat[0], bf_three));
    vs_old[iy + 3] = bf256_add(vs_old[iy + 3], bf_x_hat[1]);
    vs_old[iy + 3] = bf256_add(vs_old[iy + 3], bf_x_hat[2]);
    vs_old[iy + 3] = bf256_add(vs_old[iy + 3], bf256_mul(bf_x_hat[3], bf_two));
  }
}

static void aes_enc_constraints_Mkey_0_256(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           vbb_t* vbb, unsigned int offset, const uint8_t* k,
                                           zk_hash_256_ctx* a0_ctx, zk_hash_256_ctx* a1_ctx) {
  unsigned int w_offset = offset / 8;
  w += w_offset;
  bf256_t s[16];
  bf256_t s_dash[16];
  bf256_t vs[16];
  bf256_t vs_old[16];
  bf256_t vs_dash[16];
  for (unsigned int i = 0; i < FAEST_256F_R; i++){
    aes_enc_forward_256_1_round(w, k, in, 0, 0, s, i);
    aes_enc_backward_256_1_round(w, k, 0, 0, out, s_dash, i);

    if (i != 0){
      memcpy(vs, vs_old, sizeof(vs));
    }
    aes_enc_forward_backward_256(vbb, offset, in, out, 1, 0, NULL, vs, vs_old, vs_dash, i);
    for (unsigned int j = 0; j < 16; j++) {
      const bf256_t tmp = bf256_mul(vs[j], vs_dash[j]);
      zk_hash_256_update(a0_ctx, tmp);
      zk_hash_256_update(a1_ctx, bf256_add(bf256_add(bf256_mul(bf256_add(s[j], vs[j]),
                                                              bf256_add(s_dash[j], vs_dash[j])),
                                                    tmp),
                                          bf256_one()));
    }
  }
}

static void aes_enc_constraints_Mkey_1_256(const uint8_t* in, const uint8_t* out, vbb_t* vbb,
                                           unsigned int offset, const uint8_t* delta,
                                           zk_hash_256_ctx* b0_ctx) {
  // Step: 11..12
  bf256_t qs[16];
  bf256_t qs_old[16];
  bf256_t qs_dash[16];

  // Step: 13..14
  for (unsigned int j = 0; j < FAEST_256F_R; j++) {
    if(j != 0){
      memcpy(qs, qs_old, sizeof(qs));
    }
    aes_enc_forward_backward_256(vbb, offset, in, out, 0, 1, delta, qs, qs_old, qs_dash, j);

    bf256_t minus_part = bf256_mul(bf256_load(delta), bf256_load(delta));
    for (unsigned int i = 0; i < 16; i++){
      zk_hash_256_update(b0_ctx, bf256_add(bf256_mul(qs[i], qs_dash[i]), minus_part));
    }
  }
}

static void aes_prove_256(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                          const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
                          const faest_paramset_t* params) {
  // Step: 1..2
  // compute on the fly

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k = alloca((FAEST_256F_R + 1) * 128 / 8);
  // bf256_t* vk = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  zk_hash_256_ctx a0_ctx;
  zk_hash_256_ctx a1_ctx;

  zk_hash_256_init(&a0_ctx, chall);
  zk_hash_256_init(&a1_ctx, chall);
  aes_key_schedule_256(w, vbb, &a0_ctx, &a1_ctx, k, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  unsigned int offset = FAEST_256F_Lke;
  aes_enc_constraints_Mkey_0_256(in, out, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // Step: 12-15
  offset = FAEST_256F_Lke + FAEST_256F_Lenc;
  aes_enc_constraints_Mkey_0_256(in + 16, out + 16, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // faest_aligned_free(vk);
  //free(k);

  // Step: 16..18
  zk_hash_256_finalize(a_tilde, &a1_ctx, bf256_load(get_vole_u(vbb) + FAEST_256F_L / 8));
  zk_hash_256_finalize(b_tilde, &a0_ctx, bf256_sum_poly_vbb(vbb, FAEST_256F_L));
}

static uint8_t* aes_verify_256(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                               const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out, uint8_t* q_tilde) {
  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 13, 21
  // bf256_t* qk = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  zk_hash_256_ctx b0_ctx;

  zk_hash_256_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_256_new(vbb, delta, &b0_ctx);

  // Step: 14
  aes_enc_constraints_Mkey_1_256(in, out, vbb, FAEST_256F_Lke, delta, &b0_ctx);
  // Step: 18
  aes_enc_constraints_Mkey_1_256(in + 16, out + 16, vbb, FAEST_256F_Lke + FAEST_256F_Lenc, delta,
                                 &b0_ctx);
  // faest_aligned_free(qk);

  // Step: 20, 21
  //uint8_t* q_tilde = malloc(FAEST_256F_LAMBDA / 8);
  zk_hash_256_finalize(q_tilde, &b0_ctx, bf256_sum_poly_vbb(vbb, FAEST_256F_L));

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}

// EM-128

static void em_enc_forward_128_1_round(const uint8_t* z, const uint8_t* x, bf128_t* bf_y, unsigned int round) { // Step: 2
  if(round == 0){
    for (unsigned int j = 0; j < 4 * FAEST_EM_128F_Nwd; j++) {
      bf_y[j] = bf128_add(bf128_byte_combine_bits(z[j]), bf128_byte_combine_bits(x[j]));
    }
  } else {
    const bf128_t bf_two   = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);

    unsigned int j = round;
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_128F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf128_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf128_byte_combine_bits(x[(i + 8 * r) / 8]);
      }

      bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_two), bf128_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf128_add(bf_z_hat[0], bf128_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf128_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf128_add(bf128_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[3]);
    }
  }
}

static void em_enc_forward_backward_128_vbb_round(vbb_t* vbb, const bf128_t* bf_x, bf128_t* bf_y_fwd, bf128_t* bf_y_bkwd, bf128_t* bf_y_bkwd_final, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, unsigned int round) {
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  
  if (round == 0){
    for (unsigned int chunk_idx = 0; chunk_idx < FAEST_EM_128F_LAMBDA; chunk_idx += 8) {
      unsigned int j   = FAEST_EM_128F_R - 1;
      unsigned int ird = chunk_idx + 32 * FAEST_EM_128F_Nwd * (j + 1);
      unsigned int r   = chunk_idx % 32 / 8;
      unsigned int c   = ((chunk_idx - 8 * r) / 32 + r) % FAEST_EM_128F_Nwd;

      // Setup VOLE
      bf128_t bf_z_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        bf_z_tilde[i] = *get_vole_aes_128(vbb, chunk_idx + i);
      }

      // Forward
      unsigned int j_fwd = chunk_idx / 8;
      bf_y_fwd[j_fwd] = bf128_byte_combine(bf_z_tilde);
      if (bf_x) {
        bf_y_fwd[j_fwd] = bf128_add(bf_y_fwd[j_fwd], bf128_byte_combine(bf_x + 8 * j_fwd));
      }
      
      // Backward final
      for (unsigned int i = 0; i < 8; ++i) {
        if (bf_x) {
          bf_z_tilde[i] = bf128_add(bf_z_tilde[i], bf_x[ird + i]);
        }
      }

      bf128_t bf_y_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        bf_y_tilde[i] = bf128_add(bf128_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                  bf_z_tilde[(i + 2) % 8]);
      }
      bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

      bf_y_bkwd_final[4 * c + r] = bf128_byte_combine(bf_y_tilde);
    }

  } else {
    const bf128_t bf_two   = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);

    unsigned int i_counter = 0;
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_128F_Nwd * round + 32 * c;
      const unsigned int iy = 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        bf128_t bf_vole[8];
        for (unsigned int k = 0; k < 8; k++){
          bf_vole[k] = *get_vole_aes_128(vbb, i + 8 * r + k);
          i_counter++;
        }
        bf_z_hat[r] = bf128_byte_combine(bf_vole);
        if (bf_x) {
          bf_x_hat[r] = bf128_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf128_zero();
        }

        // Backward
        unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
        unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % FAEST_EM_128F_Nwd;
        bf128_t bf_y_tilde[8];
        for (unsigned int k = 0; k < 8; k++) {
          bf_y_tilde[k] = bf128_add(bf128_add(bf_vole[(k + 7) % 8], bf_vole[(k + 5) % 8]),
                                    bf_vole[(k + 2) % 8]);
        }
        bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

        bf_y_bkwd[4 * c_bkwd + r_bkwd] = bf128_byte_combine(bf_y_tilde);
      }
      // Foward
      bf_y_fwd[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_two), bf128_mul(bf_z_hat[1], bf_three));
      bf_y_fwd[iy + 0] = bf128_add(bf_y_fwd[iy + 0], bf_z_hat[2]);
      bf_y_fwd[iy + 0] = bf128_add(bf_y_fwd[iy + 0], bf_z_hat[3]);
      bf_y_fwd[iy + 0] = bf128_add(bf_y_fwd[iy + 0], bf_x_hat[0]);

      bf_y_fwd[iy + 1] = bf128_add(bf_z_hat[0], bf128_mul(bf_z_hat[1], bf_two));
      bf_y_fwd[iy + 1] = bf128_add(bf_y_fwd[iy + 1], bf128_mul(bf_z_hat[2], bf_three));
      bf_y_fwd[iy + 1] = bf128_add(bf_y_fwd[iy + 1], bf_z_hat[3]);
      bf_y_fwd[iy + 1] = bf128_add(bf_y_fwd[iy + 1], bf_x_hat[1]);

      bf_y_fwd[iy + 2] = bf128_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y_fwd[iy + 2] = bf128_add(bf_y_fwd[iy + 2], bf128_mul(bf_z_hat[2], bf_two));
      bf_y_fwd[iy + 2] = bf128_add(bf_y_fwd[iy + 2], bf128_mul(bf_z_hat[3], bf_three));
      bf_y_fwd[iy + 2] = bf128_add(bf_y_fwd[iy + 2], bf_x_hat[2]);

      bf_y_fwd[iy + 3] = bf128_add(bf128_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y_fwd[iy + 3] = bf128_add(bf_y_fwd[iy + 3], bf_z_hat[2]);
      bf_y_fwd[iy + 3] = bf128_add(bf_y_fwd[iy + 3], bf128_mul(bf_z_hat[3], bf_two));
      bf_y_fwd[iy + 3] = bf128_add(bf_y_fwd[iy + 3], bf_x_hat[3]);
    }
  }
}

static bf128_t get_bf_x_128(const uint8_t* x, const uint8_t* delta,  unsigned int i) {
  const bf128_t bf_delta = bf128_load(delta);
  return bf128_mul_bit(bf_delta, ptr_get_bit(x, i));
}

static void em_enc_backward_128_1_round(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf128_t* y_out, unsigned int round) {
  // only called with Mtag == Mkey == 0
  unsigned int j = round;
  for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      const unsigned int icol = (c - r + FAEST_EM_128F_Nwd) % FAEST_EM_128F_Nwd;
      const unsigned int ird =
          FAEST_EM_128F_LAMBDA + 32 * FAEST_EM_128F_Nwd * j + 32 * icol + 8 * r;
      uint8_t z_tilde = 0;
      if (j < (FAEST_EM_128F_R - 1)) {
        z_tilde = z[ird / 8];
      } else {
        z_tilde = z_out[(ird - 32 * FAEST_EM_128F_Nwd * (j + 1)) / 8] ^ x[ird / 8];
      }

      // (bit spliced)
      // delta is always bot
      // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
      const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

      // Step: 18
      y_out[4 * c + r] = bf128_byte_combine_bits(y_tilde);
    }
  }
}

static void em_enc_forward_backward_128_vbb_verify_round(vbb_t* vbb, const bf8_t* x, const uint8_t* out, const uint8_t* delta, bf128_t* bf_y_fwd, bf128_t* bf_y_bkwd, bf128_t* bf_q_out, unsigned int round) {
  // Step: 2
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, 1), bf128_from_bit(1 ^ 1)), 1 ^ 0);
  if(round ==  FAEST_EM_128F_R ){
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
        unsigned int j          = FAEST_EM_128F_R - 1;
        const unsigned int icol = (c - r + FAEST_EM_128F_Nwd) % FAEST_EM_128F_Nwd;
        const unsigned int ird =
            FAEST_EM_128F_LAMBDA + 32 * FAEST_EM_128F_Nwd * j + 32 * icol + 8 * r;

        bf128_t bf_z_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          // Step: 12
          bf_z_tilde[i] = bf_q_out[ird - 32 * FAEST_EM_128F_Nwd * (j + 1) + i];
          bf_z_tilde[i] = bf128_add(bf_z_tilde[i], get_bf_x_128(x, delta, ird + i));
        }

        bf128_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf128_add(bf128_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                    bf_z_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

        // Step: 18
        bf_y_bkwd[4 * c + r] = bf128_byte_combine(bf_y_tilde);
      }
    }
    return;
  }
  
  if(round == 0){
    for (unsigned int j = 0; j < 4 * FAEST_EM_128F_Nwd; j++) {
      // Setup VOLE
      bf128_t bf_vole[8];
      for (unsigned int i = 0; i < 8; i++) {
        bf_vole[i] = *get_vole_aes_128(vbb, 8 * j + i);
        bf_q_out[i+j*8] =
          bf128_add(bf128_mul_bit(bf_delta, ptr_get_bit(out, 8*j+i)), bf_vole[i]);
      }
      
      bf_y_fwd[j] = bf128_byte_combine(bf_vole);
      bf128_t bf_x_arr[8];
      for (unsigned int i = 0; i < 8; i++) {
        bf_x_arr[i] = get_bf_x_128(x, delta, 8 * j + i);
      }
      bf_y_fwd[j] = bf128_add(bf_y_fwd[j], bf128_byte_combine(&bf_x_arr[0]));
    }
  }else {
    const bf128_t bf_two   = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);

    unsigned int i_counter = 0;
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_128F_Nwd * round + 32 * c;
      const unsigned int iy = 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        bf128_t bf_vole[8];
        for (unsigned int k = 0; k < 8; k++){
          bf_vole[k] = *get_vole_aes_128(vbb, i + 8 * r + k);
          i_counter++;
        }
        // Forward 
        bf_z_hat[r] = bf128_byte_combine(bf_vole);
        bf128_t bf_x_tmp[8];
        for(int k = 0; k < 8; k++){
          bf_x_tmp[k] = get_bf_x_128(x, delta, (i +8*r) + k);
        }
        bf_x_hat[r] = bf128_byte_combine(&bf_x_tmp[0]);

        // Backward
        unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
        unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % FAEST_EM_128F_Nwd;
        bf128_t bf_y_tilde[8];
        for (unsigned int k = 0; k < 8; k++) {
          bf_y_tilde[k] = bf128_add(bf128_add(bf_vole[(k + 7) % 8], bf_vole[(k + 5) % 8]),
                                    bf_vole[(k + 2) % 8]);
        }
        bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);
        bf_y_bkwd[4 * c_bkwd + r_bkwd] = bf128_byte_combine(bf_y_tilde);
      }

      bf_y_fwd[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_two), bf128_mul(bf_z_hat[1], bf_three));
      bf_y_fwd[iy + 0] = bf128_add(bf_y_fwd[iy + 0], bf_z_hat[2]);
      bf_y_fwd[iy + 0] = bf128_add(bf_y_fwd[iy + 0], bf_z_hat[3]);
      bf_y_fwd[iy + 0] = bf128_add(bf_y_fwd[iy + 0], bf_x_hat[0]);

      bf_y_fwd[iy + 1] = bf128_add(bf_z_hat[0], bf128_mul(bf_z_hat[1], bf_two));
      bf_y_fwd[iy + 1] = bf128_add(bf_y_fwd[iy + 1], bf128_mul(bf_z_hat[2], bf_three));
      bf_y_fwd[iy + 1] = bf128_add(bf_y_fwd[iy + 1], bf_z_hat[3]);
      bf_y_fwd[iy + 1] = bf128_add(bf_y_fwd[iy + 1], bf_x_hat[1]);

      bf_y_fwd[iy + 2] = bf128_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y_fwd[iy + 2] = bf128_add(bf_y_fwd[iy + 2], bf128_mul(bf_z_hat[2], bf_two));
      bf_y_fwd[iy + 2] = bf128_add(bf_y_fwd[iy + 2], bf128_mul(bf_z_hat[3], bf_three));
      bf_y_fwd[iy + 2] = bf128_add(bf_y_fwd[iy + 2], bf_x_hat[2]);

      bf_y_fwd[iy + 3] = bf128_add(bf128_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y_fwd[iy + 3] = bf128_add(bf_y_fwd[iy + 3], bf_z_hat[2]);
      bf_y_fwd[iy + 3] = bf128_add(bf_y_fwd[iy + 3], bf128_mul(bf_z_hat[3], bf_two));
      bf_y_fwd[iy + 3] = bf128_add(bf_y_fwd[iy + 3], bf_x_hat[3]);
    }
  }
}

static void em_enc_constraints_Mkey_0_128(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          vbb_t* vbb, zk_hash_128_ctx* a0_ctx,
                                          zk_hash_128_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_128F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf128_t bf_s_dash[16];
  bf128_t bf_s[16];
  bf128_t bf_vs[16];
  bf128_t bf_vs_next[16];
  bf128_t bf_vs_dash[16];
  bf128_t bf_vs_dash_final[16];
  em_enc_forward_backward_128_vbb_round(vbb, NULL, bf_vs, bf_vs_dash, bf_vs_dash_final, 1, 0, NULL, 0);
  for (unsigned int i = 0; i < FAEST_EM_128F_Senc/16; i++){
    em_enc_forward_128_1_round(w, x, bf_s, i);
    em_enc_backward_128_1_round(w, x, w_out, bf_s_dash, i);
    if(i != FAEST_EM_128F_R - 1){
      em_enc_forward_backward_128_vbb_round(vbb, NULL, bf_vs_next, bf_vs_dash, NULL, 1, 0, NULL, i+1);
    }else{
      memcpy(bf_vs_dash, bf_vs_dash_final, sizeof(bf128_t) * 16);
    }

    for (unsigned int j = 0; j < 16; j++) {
      const bf128_t tmp = bf128_mul(bf_vs[j], bf_vs_dash[j]);
      zk_hash_128_update(a0_ctx, tmp);
      zk_hash_128_update(a1_ctx,
                        bf128_add(bf128_add(bf128_mul(bf128_add(bf_s[j], bf_vs[j]),
                                                      bf128_add(bf_s_dash[j], bf_vs_dash[j])),
                                            tmp),
                                  bf128_one()));
    }
    memcpy(bf_vs, bf_vs_next, sizeof(bf128_t) * 16);
  }
}

static void em_enc_constraints_Mkey_1_128(const uint8_t* out, const uint8_t* x, vbb_t* vbb,
                                          const uint8_t* delta, zk_hash_128_ctx* b0_ctx) {
  // Step: 18, 19
  // TODO: compute these on demand in em_enc_backward_128
  const bf128_t bf_delta = bf128_load(delta);
  /*
  bf128_t* bf_x = alloca(sizeof(bf128_t) * 128 * (FAEST_EM_128F_R + 1));
  for (unsigned int i = 0; i < 128 * (FAEST_EM_128F_R + 1); i++) {
    bf_x[i] = bf128_mul_bit(bf_delta, ptr_get_bit(x, i));
  }
  */

  // Step 21
  bf128_t* bf_q_out = alloca(sizeof(bf128_t) * FAEST_EM_128F_LAMBDA);
  bf128_t minus_part = bf128_mul(bf_delta, bf_delta);

  bf128_t bf_qs[FAEST_EM_128F_Nwd*4];
  bf128_t bf_qs_next[FAEST_EM_128F_Nwd*4];
  bf128_t bf_qs_dash[FAEST_EM_128F_Nwd*4];

  em_enc_forward_backward_128_vbb_verify_round(vbb, x, out, delta, bf_qs, NULL, bf_q_out,0);
  for(int i = 0; i < FAEST_EM_128F_R; i++){
    em_enc_forward_backward_128_vbb_verify_round(vbb, x, out, delta, bf_qs_next, bf_qs_dash, bf_q_out, i+1);
    for (unsigned int j = 0; j < FAEST_EM_128F_Nwd*4; j++) {
      zk_hash_128_update(b0_ctx, bf128_add(bf128_mul(bf_qs[j], bf_qs_dash[j]), minus_part));
    }
    memcpy(bf_qs, bf_qs_next, sizeof(bf128_t) * FAEST_EM_128F_Nwd*4);
  }
}

static void em_prove_128(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                         const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde) {
  // copy expanded key in to an array
  uint8_t x[FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1)];
  {
    aes_round_keys_t round_keys;
    aes128_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_128F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_128F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_128_ctx a0_ctx;
  zk_hash_128_ctx a1_ctx;

  zk_hash_128_init(&a0_ctx, chall);
  zk_hash_128_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_128(out, x, w, vbb, &a0_ctx, &a1_ctx);

  zk_hash_128_finalize(a_tilde, &a1_ctx, bf128_load(get_vole_u(vbb) + FAEST_EM_128F_Lenc / 8));
  zk_hash_128_finalize(b_tilde, &a0_ctx, bf128_sum_poly_vbb(vbb, FAEST_EM_128F_Lenc));
}

static uint8_t* em_verify_128(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                              const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out, uint8_t* q_tilde) {
  const uint8_t* delta = chall_3;

  // copy expanded key in to an array
  uint8_t x[FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    aes128_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_128F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_128F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_128_ctx b0_ctx;
  zk_hash_128_init(&b0_ctx, chall_2);
  em_enc_constraints_Mkey_1_128(out, x, vbb, delta, &b0_ctx);

  //uint8_t* q_tilde = malloc(FAEST_EM_128F_LAMBDA / 8);
  zk_hash_128_finalize(q_tilde, &b0_ctx, bf128_sum_poly_vbb(vbb, FAEST_EM_128F_Lenc));

  bf128_t bf_qtilde = bf128_load(q_tilde);
  bf128_store(q_tilde, bf128_add(bf_qtilde, bf128_mul(bf128_load(a_tilde), bf128_load(delta))));

  return q_tilde;
}

// EM-192

static void em_enc_forward_192_1_round(const uint8_t* z, const uint8_t* x, bf192_t* bf_y, unsigned int round) { // Step: 2
  if(round == 0){
    for (unsigned int j = 0; j < 4 * FAEST_EM_192F_Nwd; j++) {
      bf_y[j] = bf192_add(bf192_byte_combine_bits(z[j]), bf192_byte_combine_bits(x[j]));
    }
  } else {
    const bf192_t bf_two   = bf192_byte_combine_bits(2);
    const bf192_t bf_three = bf192_byte_combine_bits(3);

    unsigned int j = round;
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_192F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf192_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf192_byte_combine_bits(x[(i + 8 * r) / 8]);
      }

      bf_y[iy + 0] = bf192_add(bf192_mul(bf_z_hat[0], bf_two), bf192_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf192_add(bf_z_hat[0], bf192_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf192_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf192_add(bf192_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[3]);
    }
  }
}

static bf192_t get_bf_x_192(const uint8_t* x, const uint8_t* delta,  unsigned int i) {
  const bf192_t bf_delta = bf192_load(delta);
  return bf192_mul_bit(bf_delta, ptr_get_bit(x, i));
}

static void em_enc_backward_192_1_round(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf192_t* y_out, unsigned int round) {
  // only called with Mtag == Mkey == 0
  unsigned int j = round;
  for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      const unsigned int icol = (c - r + FAEST_EM_192F_Nwd) % FAEST_EM_192F_Nwd;
      const unsigned int ird =
          FAEST_EM_192F_LAMBDA + 32 * FAEST_EM_192F_Nwd * j + 32 * icol + 8 * r;
      uint8_t z_tilde = 0;
      if (j < (FAEST_EM_192F_R - 1)) {
        z_tilde = z[ird / 8];
      } else {
        z_tilde = z_out[(ird - 32 * FAEST_EM_192F_Nwd * (j + 1)) / 8] ^ x[ird / 8];
      }

      // (bit spliced)
      // delta is always bot
      // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
      const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

      // Step: 18
      y_out[4 * c + r] = bf192_byte_combine_bits(y_tilde);
    }
  }
}

static void em_enc_forward_backward_192_vbb_verify_round(vbb_t* vbb, const bf8_t* x, const uint8_t* out, const uint8_t* delta, bf192_t* bf_y_fwd, bf192_t* bf_y_bkwd, bf192_t* bf_q_out, unsigned int round) {
  // Step: 2
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, 1), bf192_from_bit(1 ^ 1)), 1 ^ 0);
  if(round ==  FAEST_EM_192F_R ){
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
        unsigned int j          = FAEST_EM_192F_R - 1;
        const unsigned int icol = (c - r + FAEST_EM_192F_Nwd) % FAEST_EM_192F_Nwd;
        const unsigned int ird =
            FAEST_EM_192F_LAMBDA + 32 * FAEST_EM_192F_Nwd * j + 32 * icol + 8 * r;

        bf192_t bf_z_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          // Step: 12
          bf_z_tilde[i] = bf_q_out[ird - 32 * FAEST_EM_192F_Nwd * (j + 1) + i];
          bf_z_tilde[i] = bf192_add(bf_z_tilde[i], get_bf_x_192(x, delta, ird + i));
        }

        bf192_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf192_add(bf192_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                    bf_z_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

        // Step: 18
        bf_y_bkwd[4 * c + r] = bf192_byte_combine(bf_y_tilde);
      }
    }
    return;
  }
  
  if(round == 0){
    for (unsigned int j = 0; j < 4 * FAEST_EM_192F_Nwd; j++) {
      // Setup VOLE
      bf192_t bf_vole[8];
      for (unsigned int i = 0; i < 8; i++) {
        bf_vole[i] = *get_vole_aes_192(vbb, 8 * j + i);
        bf_q_out[i+j*8] =
          bf192_add(bf192_mul_bit(bf_delta, ptr_get_bit(out, 8*j+i)), bf_vole[i]);
      }
      
      bf_y_fwd[j] = bf192_byte_combine(bf_vole);
      bf192_t bf_x_arr[8];
      for (unsigned int i = 0; i < 8; i++) {
        bf_x_arr[i] = get_bf_x_192(x, delta, 8 * j + i);
      }
      bf_y_fwd[j] = bf192_add(bf_y_fwd[j], bf192_byte_combine(&bf_x_arr[0]));
    }
  }else {
    const bf192_t bf_two   = bf192_byte_combine_bits(2);
    const bf192_t bf_three = bf192_byte_combine_bits(3);

    unsigned int i_counter = 0;
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_192F_Nwd * round + 32 * c;
      const unsigned int iy = 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        bf192_t bf_vole[8];
        for (unsigned int k = 0; k < 8; k++){
          bf_vole[k] = *get_vole_aes_192(vbb, i + 8 * r + k);
          i_counter++;
        }
        // Forward 
        bf_z_hat[r] = bf192_byte_combine(bf_vole);
        bf192_t bf_x_tmp[8];
        for(int k = 0; k < 8; k++){
          bf_x_tmp[k] = get_bf_x_192(x, delta, (i +8*r) + k);
        }
        bf_x_hat[r] = bf192_byte_combine(&bf_x_tmp[0]);

        // Backward
        unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
        unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % FAEST_EM_192F_Nwd;
        bf192_t bf_y_tilde[8];
        for (unsigned int k = 0; k < 8; k++) {
          bf_y_tilde[k] = bf192_add(bf192_add(bf_vole[(k + 7) % 8], bf_vole[(k + 5) % 8]),
                                    bf_vole[(k + 2) % 8]);
        }
        bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);
        bf_y_bkwd[4 * c_bkwd + r_bkwd] = bf192_byte_combine(bf_y_tilde);
      }

      bf_y_fwd[iy + 0] = bf192_add(bf192_mul(bf_z_hat[0], bf_two), bf192_mul(bf_z_hat[1], bf_three));
      bf_y_fwd[iy + 0] = bf192_add(bf_y_fwd[iy + 0], bf_z_hat[2]);
      bf_y_fwd[iy + 0] = bf192_add(bf_y_fwd[iy + 0], bf_z_hat[3]);
      bf_y_fwd[iy + 0] = bf192_add(bf_y_fwd[iy + 0], bf_x_hat[0]);

      bf_y_fwd[iy + 1] = bf192_add(bf_z_hat[0], bf192_mul(bf_z_hat[1], bf_two));
      bf_y_fwd[iy + 1] = bf192_add(bf_y_fwd[iy + 1], bf192_mul(bf_z_hat[2], bf_three));
      bf_y_fwd[iy + 1] = bf192_add(bf_y_fwd[iy + 1], bf_z_hat[3]);
      bf_y_fwd[iy + 1] = bf192_add(bf_y_fwd[iy + 1], bf_x_hat[1]);

      bf_y_fwd[iy + 2] = bf192_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y_fwd[iy + 2] = bf192_add(bf_y_fwd[iy + 2], bf192_mul(bf_z_hat[2], bf_two));
      bf_y_fwd[iy + 2] = bf192_add(bf_y_fwd[iy + 2], bf192_mul(bf_z_hat[3], bf_three));
      bf_y_fwd[iy + 2] = bf192_add(bf_y_fwd[iy + 2], bf_x_hat[2]);

      bf_y_fwd[iy + 3] = bf192_add(bf192_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y_fwd[iy + 3] = bf192_add(bf_y_fwd[iy + 3], bf_z_hat[2]);
      bf_y_fwd[iy + 3] = bf192_add(bf_y_fwd[iy + 3], bf192_mul(bf_z_hat[3], bf_two));
      bf_y_fwd[iy + 3] = bf192_add(bf_y_fwd[iy + 3], bf_x_hat[3]);
    }
  }
}

static void em_enc_forward_backward_192_vbb_round(vbb_t* vbb, const bf192_t* bf_x, bf192_t* bf_y_fwd, bf192_t* bf_y_bkwd, bf192_t* bf_y_bkwd_final, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, unsigned int round) {
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  
  if (round == 0){
    for (unsigned int chunk_idx = 0; chunk_idx < FAEST_EM_192F_LAMBDA; chunk_idx += 8) {
      unsigned int j   = FAEST_EM_192F_R - 1;
      unsigned int ird = chunk_idx + 32 * FAEST_EM_192F_Nwd * (j + 1);
      unsigned int r   = chunk_idx % 32 / 8;
      unsigned int c   = ((chunk_idx - 8 * r) / 32 + r) % FAEST_EM_192F_Nwd;

      // Setup VOLE
      bf192_t bf_z_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        bf_z_tilde[i] = *get_vole_aes_192(vbb, chunk_idx + i);
      }

      // Forward
      unsigned int j_fwd = chunk_idx / 8;
      bf_y_fwd[j_fwd] = bf192_byte_combine(bf_z_tilde);
      if (bf_x) {
        bf_y_fwd[j_fwd] = bf192_add(bf_y_fwd[j_fwd], bf192_byte_combine(bf_x + 8 * j_fwd));
      }
      
      // Backward final
      for (unsigned int i = 0; i < 8; ++i) {
        if (bf_x) {
          bf_z_tilde[i] = bf192_add(bf_z_tilde[i], bf_x[ird + i]);
        }
      }

      bf192_t bf_y_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        bf_y_tilde[i] = bf192_add(bf192_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                  bf_z_tilde[(i + 2) % 8]);
      }
      bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

      bf_y_bkwd_final[4 * c + r] = bf192_byte_combine(bf_y_tilde);
    }

  } else {
    const bf192_t bf_two   = bf192_byte_combine_bits(2);
    const bf192_t bf_three = bf192_byte_combine_bits(3);

    unsigned int i_counter = 0;
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_192F_Nwd * round + 32 * c;
      const unsigned int iy = 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        bf192_t bf_vole[8];
        for (unsigned int k = 0; k < 8; k++){
          bf_vole[k] = *get_vole_aes_192(vbb, i + 8 * r + k);
          i_counter++;
        }
        bf_z_hat[r] = bf192_byte_combine(bf_vole);
        if (bf_x) {
          bf_x_hat[r] = bf192_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf192_zero();
        }

        // Backward
        unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
        unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % FAEST_EM_192F_Nwd;
        bf192_t bf_y_tilde[8];
        for (unsigned int k = 0; k < 8; k++) {
          bf_y_tilde[k] = bf192_add(bf192_add(bf_vole[(k + 7) % 8], bf_vole[(k + 5) % 8]),
                                    bf_vole[(k + 2) % 8]);
        }
        bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

        bf_y_bkwd[4 * c_bkwd + r_bkwd] = bf192_byte_combine(bf_y_tilde);
      }
      // Foward
      bf_y_fwd[iy + 0] = bf192_add(bf192_mul(bf_z_hat[0], bf_two), bf192_mul(bf_z_hat[1], bf_three));
      bf_y_fwd[iy + 0] = bf192_add(bf_y_fwd[iy + 0], bf_z_hat[2]);
      bf_y_fwd[iy + 0] = bf192_add(bf_y_fwd[iy + 0], bf_z_hat[3]);
      bf_y_fwd[iy + 0] = bf192_add(bf_y_fwd[iy + 0], bf_x_hat[0]);

      bf_y_fwd[iy + 1] = bf192_add(bf_z_hat[0], bf192_mul(bf_z_hat[1], bf_two));
      bf_y_fwd[iy + 1] = bf192_add(bf_y_fwd[iy + 1], bf192_mul(bf_z_hat[2], bf_three));
      bf_y_fwd[iy + 1] = bf192_add(bf_y_fwd[iy + 1], bf_z_hat[3]);
      bf_y_fwd[iy + 1] = bf192_add(bf_y_fwd[iy + 1], bf_x_hat[1]);

      bf_y_fwd[iy + 2] = bf192_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y_fwd[iy + 2] = bf192_add(bf_y_fwd[iy + 2], bf192_mul(bf_z_hat[2], bf_two));
      bf_y_fwd[iy + 2] = bf192_add(bf_y_fwd[iy + 2], bf192_mul(bf_z_hat[3], bf_three));
      bf_y_fwd[iy + 2] = bf192_add(bf_y_fwd[iy + 2], bf_x_hat[2]);

      bf_y_fwd[iy + 3] = bf192_add(bf192_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y_fwd[iy + 3] = bf192_add(bf_y_fwd[iy + 3], bf_z_hat[2]);
      bf_y_fwd[iy + 3] = bf192_add(bf_y_fwd[iy + 3], bf192_mul(bf_z_hat[3], bf_two));
      bf_y_fwd[iy + 3] = bf192_add(bf_y_fwd[iy + 3], bf_x_hat[3]);
    }
  }
}

static void em_enc_constraints_Mkey_0_192(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          vbb_t* vbb, zk_hash_192_ctx* a0_ctx,
                                          zk_hash_192_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_192F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf192_t bf_s_dash[4*FAEST_EM_192F_Nwd];
  bf192_t bf_s[4*FAEST_EM_192F_Nwd];
  bf192_t bf_vs[4*FAEST_EM_192F_Nwd];
  bf192_t bf_vs_next[4*FAEST_EM_192F_Nwd];
  bf192_t bf_vs_dash[4*FAEST_EM_192F_Nwd];
  bf192_t bf_vs_dash_final[4*FAEST_EM_192F_Nwd];
  em_enc_forward_backward_192_vbb_round(vbb, NULL, bf_vs, bf_vs_dash, bf_vs_dash_final, 1, 0, NULL, 0);
  for (unsigned int i = 0; i < FAEST_EM_192F_Senc/(4*FAEST_EM_192F_Nwd); i++){
    em_enc_forward_192_1_round(w, x, bf_s, i);
    em_enc_backward_192_1_round(w, x, w_out, bf_s_dash, i);
    if(i != FAEST_EM_192F_R - 1){
      em_enc_forward_backward_192_vbb_round(vbb, NULL, bf_vs_next, bf_vs_dash, NULL, 1, 0, NULL, i+1);
    }else{
      memcpy(bf_vs_dash, bf_vs_dash_final, sizeof(bf192_t) * 4*FAEST_EM_192F_Nwd);
    }

    for (unsigned int j = 0; j < (4*FAEST_EM_192F_Nwd); j++) {
      const bf192_t tmp = bf192_mul(bf_vs[j], bf_vs_dash[j]);
      zk_hash_192_update(a0_ctx, tmp);
      zk_hash_192_update(a1_ctx,
                        bf192_add(bf192_add(bf192_mul(bf192_add(bf_s[j], bf_vs[j]),
                                                      bf192_add(bf_s_dash[j], bf_vs_dash[j])),
                                            tmp),
                                  bf192_one()));
    }
    memcpy(bf_vs, bf_vs_next, sizeof(bf192_t) * 4*FAEST_EM_192F_Nwd);
  }
}

static void em_enc_constraints_Mkey_1_192(const uint8_t* out, const uint8_t* x, vbb_t* vbb,
                                          const uint8_t* delta, zk_hash_192_ctx* b0_ctx) {
  // Step: 18, 19
  // TODO: compute these on demand in em_enc_backward_192
  const bf192_t bf_delta = bf192_load(delta);
  /*
  bf192_t* bf_x = alloca(sizeof(bf192_t) * 192 * (FAEST_EM_192F_R + 1));
  for (unsigned int i = 0; i < 192 * (FAEST_EM_192F_R + 1); i++) {
    bf_x[i] = bf192_mul_bit(bf_delta, ptr_get_bit(x, i));
  }
  */

  // Step 21
  bf192_t* bf_q_out = alloca(sizeof(bf192_t) * FAEST_EM_192F_LAMBDA);
  bf192_t minus_part = bf192_mul(bf_delta, bf_delta);

  bf192_t bf_qs[FAEST_EM_192F_Nwd*4];
  bf192_t bf_qs_next[FAEST_EM_192F_Nwd*4];
  bf192_t bf_qs_dash[FAEST_EM_192F_Nwd*4];

  em_enc_forward_backward_192_vbb_verify_round(vbb, x, out, delta, bf_qs, NULL, bf_q_out,0);
  for(int i = 0; i < FAEST_EM_192F_R; i++){
    em_enc_forward_backward_192_vbb_verify_round(vbb, x, out, delta, bf_qs_next, bf_qs_dash, bf_q_out, i+1);
    for (unsigned int j = 0; j < FAEST_EM_192F_Nwd*4; j++) {
      zk_hash_192_update(b0_ctx, bf192_add(bf192_mul(bf_qs[j], bf_qs_dash[j]), minus_part));
    }
    memcpy(bf_qs, bf_qs_next, sizeof(bf192_t) * FAEST_EM_192F_Nwd*4);
  }
}

static void em_prove_192(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                         const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde) {
  // copy expanded key in to an array
  uint8_t x[FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael192_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_192F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_192F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_192_ctx a0_ctx;
  zk_hash_192_ctx a1_ctx;

  zk_hash_192_init(&a0_ctx, chall);
  zk_hash_192_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_192(out, x, w, vbb, &a0_ctx, &a1_ctx);

  zk_hash_192_finalize(a_tilde, &a1_ctx, bf192_load(get_vole_u(vbb) + FAEST_EM_192F_Lenc / 8));
  zk_hash_192_finalize(b_tilde, &a0_ctx, bf192_sum_poly_vbb(vbb, FAEST_EM_192F_Lenc));
}

static uint8_t* em_verify_192(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                              const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out, uint8_t* q_tilde) {
  const uint8_t* delta = chall_3;

  // copy expanded key in to an array
  uint8_t x[FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael192_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_192F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_192F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_192_ctx b0_ctx;
  zk_hash_192_init(&b0_ctx, chall_2);
  em_enc_constraints_Mkey_1_192(out, x, vbb, delta, &b0_ctx);

  //uint8_t* q_tilde = malloc(FAEST_EM_192F_LAMBDA / 8);
  zk_hash_192_finalize(q_tilde, &b0_ctx, bf192_sum_poly_vbb(vbb, FAEST_EM_192F_Lenc));

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

// EM-256

static void em_enc_forward_256_1_round(const uint8_t* z, const uint8_t* x, bf256_t* bf_y, unsigned int round) { // Step: 2
  if(round == 0){
    for (unsigned int j = 0; j < 4 * FAEST_EM_256F_Nwd; j++) {
      bf_y[j] = bf256_add(bf256_byte_combine_bits(z[j]), bf256_byte_combine_bits(x[j]));
    }
  } else {
    const bf256_t bf_two   = bf256_byte_combine_bits(2);
    const bf256_t bf_three = bf256_byte_combine_bits(3);

    unsigned int j = round;
    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_256F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf256_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf256_byte_combine_bits(x[(i + 8 * r) / 8]);
      }

      bf_y[iy + 0] = bf256_add(bf256_mul(bf_z_hat[0], bf_two), bf256_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf256_add(bf_z_hat[0], bf256_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf256_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf256_add(bf256_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[3]);
    }
  }
}

static bf256_t get_bf_x_256(const uint8_t* x, const uint8_t* delta,  unsigned int i) {
  const bf256_t bf_delta = bf256_load(delta);
  return bf256_mul_bit(bf_delta, ptr_get_bit(x, i));
}

static void em_enc_forward_backward_256_vbb_round(vbb_t* vbb, const bf256_t* bf_x, bf256_t* bf_y_fwd, bf256_t* bf_y_bkwd, bf256_t* bf_y_bkwd_final, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, unsigned int round) {
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  
  if (round == 0){
    for (unsigned int chunk_idx = 0; chunk_idx < FAEST_EM_256F_LAMBDA; chunk_idx += 8) {
      unsigned int j   = FAEST_EM_256F_R - 1;
      unsigned int ird = chunk_idx + 32 * FAEST_EM_256F_Nwd * (j + 1);
      unsigned int r   = chunk_idx % 32 / 8;
      unsigned int c   = ((chunk_idx - 8 * r) / 32 + r) % FAEST_EM_256F_Nwd;

      if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
        c = (c + 1) % FAEST_EM_256F_Nwd;
      }

      // Setup VOLE
      bf256_t bf_z_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        bf_z_tilde[i] = *get_vole_aes_256(vbb, chunk_idx + i);
      }

      // Forward
      unsigned int j_fwd = chunk_idx / 8;
      bf_y_fwd[j_fwd] = bf256_byte_combine(bf_z_tilde);
      if (bf_x) {
        bf_y_fwd[j_fwd] = bf256_add(bf_y_fwd[j_fwd], bf256_byte_combine(bf_x + 8 * j_fwd));
      }
      
      // Backward final
      for (unsigned int i = 0; i < 8; ++i) {
        if (bf_x) {
          bf_z_tilde[i] = bf256_add(bf_z_tilde[i], bf_x[ird + i]);
        }
      }

      bf256_t bf_y_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        bf_y_tilde[i] = bf256_add(bf256_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                  bf_z_tilde[(i + 2) % 8]);
      }
      bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

      bf_y_bkwd_final[4 * c + r] = bf256_byte_combine(bf_y_tilde);
    }

  } else {
    const bf256_t bf_two   = bf256_byte_combine_bits(2);
    const bf256_t bf_three = bf256_byte_combine_bits(3);

    unsigned int i_counter = 0;
    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_256F_Nwd * round + 32 * c;
      const unsigned int iy = 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        bf256_t bf_vole[8];
        for (unsigned int k = 0; k < 8; k++){
          bf_vole[k] = *get_vole_aes_256(vbb, i + 8 * r + k);
          i_counter++;
        }
        bf_z_hat[r] = bf256_byte_combine(bf_vole);
        if (bf_x) {
          bf_x_hat[r] = bf256_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf256_zero();
        }

        // Backward
        unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
        unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % FAEST_EM_256F_Nwd;
        if (FAEST_EM_256F_Nwd == 8 && r_bkwd >= 2) {
          c_bkwd = (c_bkwd + 1) % FAEST_EM_256F_Nwd;
        }
        bf256_t bf_y_tilde[8];
        for (unsigned int k = 0; k < 8; k++) {
          bf_y_tilde[k] = bf256_add(bf256_add(bf_vole[(k + 7) % 8], bf_vole[(k + 5) % 8]),
                                    bf_vole[(k + 2) % 8]);
        }
        bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

        bf_y_bkwd[4 * c_bkwd + r_bkwd] = bf256_byte_combine(bf_y_tilde);
      }
      // Foward
      bf_y_fwd[iy + 0] = bf256_add(bf256_mul(bf_z_hat[0], bf_two), bf256_mul(bf_z_hat[1], bf_three));
      bf_y_fwd[iy + 0] = bf256_add(bf_y_fwd[iy + 0], bf_z_hat[2]);
      bf_y_fwd[iy + 0] = bf256_add(bf_y_fwd[iy + 0], bf_z_hat[3]);
      bf_y_fwd[iy + 0] = bf256_add(bf_y_fwd[iy + 0], bf_x_hat[0]);

      bf_y_fwd[iy + 1] = bf256_add(bf_z_hat[0], bf256_mul(bf_z_hat[1], bf_two));
      bf_y_fwd[iy + 1] = bf256_add(bf_y_fwd[iy + 1], bf256_mul(bf_z_hat[2], bf_three));
      bf_y_fwd[iy + 1] = bf256_add(bf_y_fwd[iy + 1], bf_z_hat[3]);
      bf_y_fwd[iy + 1] = bf256_add(bf_y_fwd[iy + 1], bf_x_hat[1]);

      bf_y_fwd[iy + 2] = bf256_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y_fwd[iy + 2] = bf256_add(bf_y_fwd[iy + 2], bf256_mul(bf_z_hat[2], bf_two));
      bf_y_fwd[iy + 2] = bf256_add(bf_y_fwd[iy + 2], bf256_mul(bf_z_hat[3], bf_three));
      bf_y_fwd[iy + 2] = bf256_add(bf_y_fwd[iy + 2], bf_x_hat[2]);

      bf_y_fwd[iy + 3] = bf256_add(bf256_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y_fwd[iy + 3] = bf256_add(bf_y_fwd[iy + 3], bf_z_hat[2]);
      bf_y_fwd[iy + 3] = bf256_add(bf_y_fwd[iy + 3], bf256_mul(bf_z_hat[3], bf_two));
      bf_y_fwd[iy + 3] = bf256_add(bf_y_fwd[iy + 3], bf_x_hat[3]);
    }
  }
}

static void em_enc_backward_256_1_round(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf256_t* y_out, unsigned int round) {
  // only called with Mtag == Mkey == 0
  unsigned int j = round;
  for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      unsigned int icol = (c - r + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
      if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
        icol = (icol - 1 + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
      }
      
      const unsigned int ird =
          FAEST_EM_256F_LAMBDA + 32 * FAEST_EM_256F_Nwd * j + 32 * icol + 8 * r;
      uint8_t z_tilde = 0;
      if (j < (FAEST_EM_256F_R - 1)) {
        z_tilde = z[ird / 8];
      } else {
        z_tilde = z_out[(ird - 32 * FAEST_EM_256F_Nwd * (j + 1)) / 8] ^ x[ird / 8];
      }

      // (bit spliced)
      // delta is always bot
      // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
      const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

      // Step: 18
      y_out[4 * c + r] = bf256_byte_combine_bits(y_tilde);
    }
  }
}

static void em_enc_forward_backward_256_vbb_verify_round(vbb_t* vbb, const bf8_t* x, const uint8_t* out, const uint8_t* delta, bf256_t* bf_y_fwd, bf256_t* bf_y_bkwd, bf256_t* bf_q_out, unsigned int round) {
  // Step: 2
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, 1), bf256_from_bit(1 ^ 1)), 1 ^ 0);
  if(round ==  FAEST_EM_256F_R ){
    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
        unsigned int j          = FAEST_EM_256F_R - 1;
        unsigned int icol = (c - r + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
        if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
          icol = (icol - 1 + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
        }
        const unsigned int ird =
            FAEST_EM_256F_LAMBDA + 32 * FAEST_EM_256F_Nwd * j + 32 * icol + 8 * r;

        bf256_t bf_z_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          // Step: 12
          bf_z_tilde[i] = bf_q_out[ird - 32 * FAEST_EM_256F_Nwd * (j + 1) + i];
          bf_z_tilde[i] = bf256_add(bf_z_tilde[i], get_bf_x_256(x, delta, ird + i));
        }

        bf256_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf256_add(bf256_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                    bf_z_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

        // Step: 18
        bf_y_bkwd[4 * c + r] = bf256_byte_combine(bf_y_tilde);
      }
    }
    return;
  }
  
  if(round == 0){
    for (unsigned int j = 0; j < 4 * FAEST_EM_256F_Nwd; j++) {
      // Setup VOLE
      bf256_t bf_vole[8];
      for (unsigned int i = 0; i < 8; i++) {
        bf_vole[i] = *get_vole_aes_256(vbb, 8 * j + i);
        bf_q_out[i+j*8] =
          bf256_add(bf256_mul_bit(bf_delta, ptr_get_bit(out, 8*j+i)), bf_vole[i]);
      }
      
      bf_y_fwd[j] = bf256_byte_combine(bf_vole);
      bf256_t bf_x_arr[8];
      for (unsigned int i = 0; i < 8; i++) {
        bf_x_arr[i] = get_bf_x_256(x, delta, 8 * j + i);
      }
      bf_y_fwd[j] = bf256_add(bf_y_fwd[j], bf256_byte_combine(&bf_x_arr[0]));
    }
  }else {
    const bf256_t bf_two   = bf256_byte_combine_bits(2);
    const bf256_t bf_three = bf256_byte_combine_bits(3);

    unsigned int i_counter = 0;
    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_256F_Nwd * round + 32 * c;
      const unsigned int iy = 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        bf256_t bf_vole[8];
        for (unsigned int k = 0; k < 8; k++){
          bf_vole[k] = *get_vole_aes_256(vbb, i + 8 * r + k);
          i_counter++;
        }
        // Forward 
        bf_z_hat[r] = bf256_byte_combine(bf_vole);
        bf256_t bf_x_tmp[8];
        for(int k = 0; k < 8; k++){
          bf_x_tmp[k] = get_bf_x_256(x, delta, (i +8*r) + k);
        }
        bf_x_hat[r] = bf256_byte_combine(&bf_x_tmp[0]);

        // Backward
        unsigned int r_bkwd = ((i_counter-1) % 32) / 8;
        unsigned int c_bkwd = ((i_counter-1 - 8 * r_bkwd) / 32 + r_bkwd) % FAEST_EM_256F_Nwd;
        if (FAEST_EM_256F_Nwd == 8 && r_bkwd >= 2) {
          c_bkwd = (c_bkwd + 1) % FAEST_EM_256F_Nwd;
        }
        bf256_t bf_y_tilde[8];
        for (unsigned int k = 0; k < 8; k++) {
          bf_y_tilde[k] = bf256_add(bf256_add(bf_vole[(k + 7) % 8], bf_vole[(k + 5) % 8]),
                                    bf_vole[(k + 2) % 8]);
        }
        bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);
        bf_y_bkwd[4 * c_bkwd + r_bkwd] = bf256_byte_combine(bf_y_tilde);
      }

      bf_y_fwd[iy + 0] = bf256_add(bf256_mul(bf_z_hat[0], bf_two), bf256_mul(bf_z_hat[1], bf_three));
      bf_y_fwd[iy + 0] = bf256_add(bf_y_fwd[iy + 0], bf_z_hat[2]);
      bf_y_fwd[iy + 0] = bf256_add(bf_y_fwd[iy + 0], bf_z_hat[3]);
      bf_y_fwd[iy + 0] = bf256_add(bf_y_fwd[iy + 0], bf_x_hat[0]);

      bf_y_fwd[iy + 1] = bf256_add(bf_z_hat[0], bf256_mul(bf_z_hat[1], bf_two));
      bf_y_fwd[iy + 1] = bf256_add(bf_y_fwd[iy + 1], bf256_mul(bf_z_hat[2], bf_three));
      bf_y_fwd[iy + 1] = bf256_add(bf_y_fwd[iy + 1], bf_z_hat[3]);
      bf_y_fwd[iy + 1] = bf256_add(bf_y_fwd[iy + 1], bf_x_hat[1]);

      bf_y_fwd[iy + 2] = bf256_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y_fwd[iy + 2] = bf256_add(bf_y_fwd[iy + 2], bf256_mul(bf_z_hat[2], bf_two));
      bf_y_fwd[iy + 2] = bf256_add(bf_y_fwd[iy + 2], bf256_mul(bf_z_hat[3], bf_three));
      bf_y_fwd[iy + 2] = bf256_add(bf_y_fwd[iy + 2], bf_x_hat[2]);

      bf_y_fwd[iy + 3] = bf256_add(bf256_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y_fwd[iy + 3] = bf256_add(bf_y_fwd[iy + 3], bf_z_hat[2]);
      bf_y_fwd[iy + 3] = bf256_add(bf_y_fwd[iy + 3], bf256_mul(bf_z_hat[3], bf_two));
      bf_y_fwd[iy + 3] = bf256_add(bf_y_fwd[iy + 3], bf_x_hat[3]);
    }
  }
}

static void em_enc_constraints_Mkey_0_256(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          vbb_t* vbb, zk_hash_256_ctx* a0_ctx,
                                          zk_hash_256_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_256F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf256_t bf_s_dash[4*FAEST_EM_256F_Nwd];
  bf256_t bf_s[4*FAEST_EM_256F_Nwd];
  bf256_t bf_vs[4*FAEST_EM_256F_Nwd];
  bf256_t bf_vs_next[4*FAEST_EM_256F_Nwd];
  bf256_t bf_vs_dash[4*FAEST_EM_256F_Nwd];
  bf256_t bf_vs_dash_final[4*FAEST_EM_256F_Nwd];
  em_enc_forward_backward_256_vbb_round(vbb, NULL, bf_vs, bf_vs_dash, bf_vs_dash_final, 1, 0, NULL, 0);
  for (unsigned int i = 0; i < FAEST_EM_256F_Senc/(4*FAEST_EM_256F_Nwd); i++){
    em_enc_forward_256_1_round(w, x, bf_s, i);
    em_enc_backward_256_1_round(w, x, w_out, bf_s_dash, i);
    if(i != FAEST_EM_256F_R - 1){
      em_enc_forward_backward_256_vbb_round(vbb, NULL, bf_vs_next, bf_vs_dash, NULL, 1, 0, NULL, i+1);
    }else{
      memcpy(bf_vs_dash, bf_vs_dash_final, sizeof(bf256_t) * 4*FAEST_EM_256F_Nwd);
    }

    for (unsigned int j = 0; j < (4*FAEST_EM_256F_Nwd); j++) {
      const bf256_t tmp = bf256_mul(bf_vs[j], bf_vs_dash[j]);
      zk_hash_256_update(a0_ctx, tmp);
      zk_hash_256_update(a1_ctx,
                        bf256_add(bf256_add(bf256_mul(bf256_add(bf_s[j], bf_vs[j]),
                                                      bf256_add(bf_s_dash[j], bf_vs_dash[j])),
                                            tmp),
                                  bf256_one()));
    }
    memcpy(bf_vs, bf_vs_next, sizeof(bf256_t) * 4*FAEST_EM_256F_Nwd);
  }
}

static void em_enc_constraints_Mkey_1_256(const uint8_t* out, const uint8_t* x, vbb_t* vbb,
                                          const uint8_t* delta, zk_hash_256_ctx* b0_ctx) {
  // Step: 18, 19
  // TODO: compute these on demand in em_enc_backward_256
  const bf256_t bf_delta = bf256_load(delta);
  /*
  bf256_t* bf_x = alloca(sizeof(bf256_t) * 256 * (FAEST_EM_256F_R + 1));
  for (unsigned int i = 0; i < 256 * (FAEST_EM_256F_R + 1); i++) {
    bf_x[i] = bf256_mul_bit(bf_delta, ptr_get_bit(x, i));
  }
  */

  // Step 21
  bf256_t* bf_q_out = alloca(sizeof(bf256_t) * FAEST_EM_256F_LAMBDA);
  bf256_t minus_part = bf256_mul(bf_delta, bf_delta);

  bf256_t bf_qs[FAEST_EM_256F_Nwd*4];
  bf256_t bf_qs_next[FAEST_EM_256F_Nwd*4];
  bf256_t bf_qs_dash[FAEST_EM_256F_Nwd*4];

  em_enc_forward_backward_256_vbb_verify_round(vbb, x, out, delta, bf_qs, NULL, bf_q_out,0);
  for(int i = 0; i < FAEST_EM_256F_R; i++){
    em_enc_forward_backward_256_vbb_verify_round(vbb, x, out, delta, bf_qs_next, bf_qs_dash, bf_q_out, i+1);
    for (unsigned int j = 0; j < FAEST_EM_256F_Nwd*4; j++) {
      zk_hash_256_update(b0_ctx, bf256_add(bf256_mul(bf_qs[j], bf_qs_dash[j]), minus_part));
    }
    memcpy(bf_qs, bf_qs_next, sizeof(bf256_t) * FAEST_EM_256F_Nwd*4);
  }
}

static void em_prove_256(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                         const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde) {
  // copy expanded key in to an array
  uint8_t x[FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael256_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_256F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_256F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_256_ctx a0_ctx;
  zk_hash_256_ctx a1_ctx;

  zk_hash_256_init(&a0_ctx, chall);
  zk_hash_256_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_256(out, x, w, vbb, &a0_ctx, &a1_ctx);

  zk_hash_256_finalize(a_tilde, &a1_ctx, bf256_load(get_vole_u(vbb) + FAEST_EM_256F_Lenc / 8));
  zk_hash_256_finalize(b_tilde, &a0_ctx, bf256_sum_poly_vbb(vbb, FAEST_EM_256F_Lenc));
}

static uint8_t* em_verify_256(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                              const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out, uint8_t* q_tilde) {
  const uint8_t* delta = chall_3;

  // copy expanded key in to an array
  uint8_t x[FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael256_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_256F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_256F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_256_ctx b0_ctx;
  zk_hash_256_init(&b0_ctx, chall_2);
  em_enc_constraints_Mkey_1_256(out, x, vbb, delta, &b0_ctx);

  //uint8_t* q_tilde = malloc(FAEST_EM_256F_LAMBDA / 8);
  zk_hash_256_finalize(q_tilde, &b0_ctx, bf256_sum_poly_vbb(vbb, FAEST_EM_256F_Lenc));

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}

// dispatchers

void aes_prove(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
               const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      aes_prove_256(w, vbb, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_256(w, vbb, in, out, chall, a_tilde, b_tilde);
    }
    break;
  case 192:
    if (params->faest_param.Lke) {
      aes_prove_192(w, vbb, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_192(w, vbb, in, out, chall, a_tilde, b_tilde);
    }
    break;
  default:
    if (params->faest_param.Lke) {
      aes_prove_128(w, vbb, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_128(w, vbb, in, out, chall, a_tilde, b_tilde);
    }
  }
}

uint8_t* aes_verify(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params, uint8_t* q_tilde) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      return aes_verify_256(vbb, chall_2, chall_3, a_tilde, in, out, q_tilde);
    } else {
      return em_verify_256(vbb, chall_2, chall_3, a_tilde, in, out, q_tilde);
    }
  case 192:
    if (params->faest_param.Lke) {
      return aes_verify_192(vbb, chall_2, chall_3, a_tilde, in, out, q_tilde);
    } else {
      return em_verify_192(vbb, chall_2, chall_3, a_tilde, in, out, q_tilde);
    }
  default:
    if (params->faest_param.Lke) {
      return aes_verify_128(vbb, chall_2, chall_3, a_tilde, in, out, q_tilde);
    } else {
      return em_verify_128(vbb, chall_2, chall_3, a_tilde, in, out, q_tilde);
    }
  }
}
