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

// dispatchers

void aes_prove(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
               const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params) {
  aes_prove_128(w, vbb, in, out, chall, a_tilde, b_tilde, params);
}

uint8_t* aes_verify(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params, uint8_t* q_tilde) {
  return aes_verify_128(vbb, chall_2, chall_3, a_tilde, in, out, q_tilde);
}
