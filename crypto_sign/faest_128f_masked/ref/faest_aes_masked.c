#include "faest.h"
#include "fields.h"
#include "vole.h"
#include "universal_hashing.h"
#include "utils.h"
#include "parameters.h"
#include "shuffle.h"

static const bf8_t Rcon[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

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

static void aes_key_schedule_backward_1_round_share(const uint8_t* x, const uint8_t* xk, uint8_t* out,
                                              unsigned int j, const faest_paramset_t* params, bool addRcon) {
  // Step: 1 skipped (sanity check)
  unsigned int lambda = params->faest_param.lambda;

  // Step: 2
  unsigned int iwd;
  if (lambda == 192){
    iwd   = 192/8 * j/4;
  } else {
    iwd   = 128/8 * j/4;
  }
  bool rmvRcon = addRcon;
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
    if(addRcon == true){
      out[c] = y_tilde ^ 0x5;
    }else{
      out[c] = y_tilde;
    }
  }
}

static void aes_key_schedule_backward_128_vbb_vk_round_share(vbb_t* vbb, uint8_t Mtag, uint8_t Mkey,
                                                 const uint8_t* delta, bf128_t* bf_out, unsigned int j, unsigned int share) {
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
      bf128_t VOLE = *get_vole_aes_128_share(vbb, (8 * (j+c) + i) + FAEST_128F_LAMBDA, share);
      add_vole_to_vk_cache_share(vbb, (8 * (j+c) + i) + FAEST_128F_LAMBDA, &VOLE, share);
      bf_x_tilde[i] = VOLE;
    }
    for(unsigned int i = 0; i < 8; i++){
      bf_x_tilde[i] = bf128_add(bf_x_tilde[i], *get_vk_128_share(vbb, iwd + 8 * c + i, share));
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

static void aes_key_schedule_128_masked(const uint8_t* w_share, vbb_t* vbb,
                                                    zk_hash_128_ctx* a0_ctx,
                                                    zk_hash_128_ctx* a1_ctx, uint8_t* k,
                                                    const faest_paramset_t* params) {
  uint8_t w_dash[2][4] = {0};
  bf128_t v_w_dash[2][4 * 8] = {0};
  aes_key_schedule_forward_1(w_share, k, params);
  aes_key_schedule_forward_1(w_share + FAEST_128F_L/8, k + (FAEST_128F_R + 1) * 128 / 8, params);

  const unsigned int Ske    = params->faest_param.Ske;

  unsigned int iwd   = 0;
  for (unsigned int j = 0; j < Ske; j+=4) {
    bf128_t bf_k_hat_share[2][4] = {0};
    bf128_t bf_v_k_hat_share[2][4] = {0};
    bf128_t bf_w_dash_hat_share[2][4] = {0};
    bf128_t bf_v_w_dash_hat_share[2][4] = {0};


    // Share 0
    aes_key_schedule_backward_1_round_share(w_share + FAEST_128F_LAMBDA / 8, k, &w_dash[0][0], j, params, false);
    aes_key_schedule_backward_128_vbb_vk_round_share(vbb, 1, 0, NULL, &v_w_dash[0][0], j, 0);
    
    uint8_t decoy_in[9] = {0};
    bf128_t decoy_out[9] = {0};
    for (unsigned int i = 0; i < 9; i++){
      decoy_in[i] = bf8_rand();
      decoy_out[i] = bf128_zero();
    }

    uint8_t permutation[4];
    for(unsigned int i = 0; i < 4; i++){
      permutation[i] = i;
    }
    shuffle_4(permutation);

    for (unsigned int i = 0; i <= 3; i++) {
      uint8_t r = permutation[i];
      
      bf_v_k_hat_share[0][(r + 3) % 4] = bf128_byte_combine_vk_share(vbb, (96 + iwd*8 + 8 * r), 0);
      bf_v_w_dash_hat_share[0][r]      = bf128_byte_combine(v_w_dash[0] + (8 * r));

      byte_combine_bits_shuffle_11(decoy_in, decoy_out, k + (96 + iwd*8 + 8 * r) / 8, w_dash[0] + (8 * r) / 8, bf_k_hat_share[0] + (r + 3) % 4, bf_w_dash_hat_share[0] + r);
    }

    // share 1
    aes_key_schedule_backward_1_round_share(w_share + FAEST_128F_LAMBDA / 8 + FAEST_128F_L/8, k + (FAEST_128F_R + 1) * 128 / 8, &w_dash[1][0], j, params, true);
    aes_key_schedule_backward_128_vbb_vk_round_share(vbb, 1, 0, NULL, &v_w_dash[1][0], j, 1);

    for (unsigned int i = 0; i < 9; i++){
      decoy_in[i] = bf8_rand();
      decoy_out[i] = bf128_zero();
    }
    shuffle_4(permutation);
    for (unsigned int i = 0; i <= 3; i++) {
      uint8_t r = permutation[i];
      bf_v_k_hat_share[1][(r + 3) % 4] = bf128_byte_combine_vk_share(vbb, (96 + iwd*8 + 8 * r), 1);
      bf_v_w_dash_hat_share[1][r]      = bf128_byte_combine(v_w_dash[1] + (8 * r));

      byte_combine_bits_shuffle_11(decoy_in, decoy_out, (k + (FAEST_128F_R + 1) * 128 / 8) +(96 + iwd*8 + 8 * r) / 8, w_dash[1]+ (8 * r) / 8, bf_k_hat_share[1]+(r + 3) % 4, bf_w_dash_hat_share[1] + r);
    }

    // hash the shares
    for (unsigned int r = 0; r <= 3; r++) {
      const bf128_t part_a = bf128_add(bf_v_k_hat_share[0][r], bf_k_hat_share[0][r]);
      const bf128_t part_b = bf128_add(bf_w_dash_hat_share[0][r], bf_v_w_dash_hat_share[0][r]);
      const bf128_t part_d = bf128_add(bf_k_hat_share[1][r], bf_v_k_hat_share[1][r]);
      const bf128_t part_c = bf128_add(bf_v_w_dash_hat_share[1][r], bf_w_dash_hat_share[1][r]);

      // instead of storing in A0, A1, hash it
      bf128_t mask1 = bf128_rand();
      bf128_t mask2 = bf128_rand();
      const bf128_t tmp_0 = bf128_add(bf128_mul(bf_v_k_hat_share[1][r], bf_v_w_dash_hat_share[0][r]), bf128_add(bf128_mul(bf_v_k_hat_share[0][r], bf_v_w_dash_hat_share[0][r]), bf128_add(bf128_mul(bf_v_k_hat_share[0][r], bf_v_w_dash_hat_share[1][r]), mask1)));
      zk_hash_128_update(a0_ctx, tmp_0);
      const bf128_t share_0 = bf128_add(bf128_mul(part_d, part_b), bf128_add(bf128_add(bf128_mul(part_a, part_b), bf128_add(bf128_mul(part_a, part_c), mask2)), tmp_0));
      zk_hash_128_update(a1_ctx, share_0);
      
      const bf128_t tmp_1 = bf128_add(bf128_mul(bf_v_k_hat_share[1][r], bf_v_w_dash_hat_share[1][r]), mask1);
      zk_hash_128_update(a0_ctx + 1, tmp_1);
      const bf128_t share_1 = bf128_add(bf128_add(bf128_add(bf128_mul(part_d, part_c), mask2), tmp_1), bf128_one());
      zk_hash_128_update(a1_ctx+ 1, share_1); 
    }
    iwd += 128 / 8;
  }
}

static void aes_enc_forward_128_1_round(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf128_t* bf_y, int round) {
  if (round == 0){
    uint8_t decoy_in[4] = {0};
    bf128_t decoy_out[4] = {0};
    for (unsigned int i = 0; i < 4; i++){
      decoy_in[i] = bf8_rand();
      decoy_out[i] = bf128_zero();
    }

    uint8_t permutation[16];
    for(unsigned int i = 0; i < 16; i++){
      permutation[i] = i;
    }
    shuffle_16(permutation);

    for (unsigned int i = 0; i < 16; i++) {
      bf128_t xin_field;
      bf128_t xk_field;
      byte_combine_bits_shuffle_6(decoy_in, decoy_out, in + permutation[i], xk + permutation[i], &xin_field, &xk_field);

      bf_y[permutation[i]] = bf128_add(xin_field, xk_field);
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
      
      
      uint8_t decoy_in[9] = {0};
      bf128_t decoy_out[9] = {0};
      for (unsigned int i = 0; i < 9; i++){
        decoy_in[i] = bf8_rand();
        decoy_out[i] = bf128_zero();
      }

      uint8_t permutation[4];
      for(unsigned int i = 0; i < 4; i++){
        permutation[i] = i;
      }
      shuffle_4(permutation);
      for (unsigned int i = 0; i <= 3; i++) {
        uint8_t r = permutation[i];

        byte_combine_bits_shuffle_11(decoy_in, decoy_out, x + (ix + 8 * r) / 8, xk + (ik + 8 * r) / 8, bf_x_hat + r, bf_xk_hat + r);
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

static void aes_enc_backward_128_1_round_share(const uint8_t* x, const uint8_t* xk, const uint8_t* out,
                                   bf128_t* y_out, unsigned int round, unsigned int share) {
  // called only with Mtag == Mkey == 0

  uint8_t xtilde;
  // Step:2..4
  unsigned int j = round;

  uint8_t decoy_in[4] = {0};
  bf128_t decoy_out[4] = {0};
  for (unsigned int i = 0; i < 4; i++){
    decoy_in[i] = bf8_rand();
    decoy_out[i] = bf128_zero();
  }

  uint8_t permutation[16];
  for(unsigned int i = 0; i < 16; i++){
    permutation[i] = i;
  }
  shuffle_16(permutation);

  for (unsigned int i = 0; i < 16; i++){
    unsigned int p = permutation[i];
    unsigned int c = p / 4;
    unsigned int r = p % 4;
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
    const uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2) ^ (share * 0x5);

    uint8_t single_decoy_in = bf8_rand();
    bf128_t single_decoy_out = bf128_zero();
    byte_combine_bits_shuffle_6(decoy_in, decoy_out, &ytilde, &single_decoy_in, y_out+4*c+r, &single_decoy_out);

  }
}

static void aes_enc_forward_backward_128_share(vbb_t* vbb, unsigned int offset, const uint8_t* in, const uint8_t* out,
                                         uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                         bf128_t* vs, bf128_t* vs_old, bf128_t* vs_dash, unsigned int round, unsigned int share) {
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
            bf_x_tilde[i] = bf128_add(bf_xout, *get_vk_128_share(vbb, 128 + ird + i, share));
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
      vs[i] = bf128_add(bf128_byte_combine(bf_xin), bf128_byte_combine_vk_share(vbb, (8 * i), share));
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
        bf128_t t = *get_vole_aes_128_share(vbb, offset + ix + 8 * r + i, share);
        memcpy(bf_x_tilde + i, &t, sizeof(bf128_t));
        i_counter++;

        // forward part
      }
      bf_x_hat[r]  =  bf128_byte_combine(bf_x_tilde);
      bf_xk_hat[r] = bf128_byte_combine_vk_share(vbb, (ik + 8 * r), share);

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

static void aes_enc_constraints_128_masked(const uint8_t* in_share, const uint8_t* out_share, const uint8_t* w_share,
                                           vbb_t* vbb, unsigned int offset, const uint8_t* k_share,
                                           zk_hash_128_ctx* a0_ctx, zk_hash_128_ctx* a1_ctx) {
  unsigned int w_offset = offset / 8;
  w_share += w_offset;

  bf128_t s_share[2][16];
  bf128_t vs_share[2][16];
  bf128_t vs_share_old[2][16];
  bf128_t s_dash_share[2][16];
  bf128_t vs_dash_share[2][16];

  for (unsigned int i = 0; i < FAEST_128F_R; i++){
    if (i != 0){
      memcpy(vs_share, vs_share_old, sizeof(vs_share));
    }

    aes_enc_forward_128_1_round(w_share, k_share, in_share, s_share[0], i);
    aes_enc_backward_128_1_round_share(w_share, k_share, out_share, s_dash_share[0], i, 0);
    aes_enc_forward_backward_128_share(vbb, offset, in_share, out_share, 1, 0, NULL, vs_share[0], vs_share_old[0], vs_dash_share[0], i, 0);
    
    aes_enc_forward_128_1_round(w_share + + FAEST_128F_L/8, k_share + (FAEST_128F_R + 1) * 128 / 8, in_share + MAX_LAMBDA_BYTES, s_share[1], i);
    aes_enc_backward_128_1_round_share(w_share + + FAEST_128F_L/8, k_share + (FAEST_128F_R + 1) * 128 / 8, out_share + 16, s_dash_share[1], i, 1);
    aes_enc_forward_backward_128_share(vbb, offset, in_share + MAX_LAMBDA_BYTES, out_share + MAX_LAMBDA_BYTES, 1, 0, NULL, vs_share[1], vs_share_old[1], vs_dash_share[1], i, 1);

    for (unsigned int j = 0; j < 16; j++){
      const bf128_t part_a = bf128_add(vs_share[0][j], s_share[0][j]);
      const bf128_t part_b = bf128_add(s_dash_share[0][j], vs_dash_share[0][j]);
      const bf128_t part_d = bf128_add(s_share[1][j], vs_share[1][j]);
      const bf128_t part_c = bf128_add(vs_dash_share[1][j], s_dash_share[1][j]);

      bf128_t mask1 = bf128_rand();
      bf128_t mask2 = bf128_rand();
      const bf128_t tmp_0 = bf128_add(bf128_mul(vs_share[1][j], vs_dash_share[0][j]), bf128_add(bf128_mul(vs_share[0][j], vs_dash_share[0][j]), bf128_add(bf128_mul(vs_share[0][j], vs_dash_share[1][j]), mask1)));
      zk_hash_128_update(a0_ctx, tmp_0);
      const bf128_t share_0 = bf128_add(bf128_mul(part_d, part_b), bf128_add(bf128_add(bf128_mul(part_a, part_b), bf128_add(bf128_mul(part_a, part_c), mask2)), tmp_0));
      zk_hash_128_update(a1_ctx, share_0);

      const bf128_t tmp_1 = bf128_add(bf128_mul(vs_share[1][j], vs_dash_share[1][j]), mask1);
      zk_hash_128_update(a0_ctx + 1, tmp_1);
      const bf128_t share_1 = bf128_add(bf128_add(bf128_add(bf128_mul(part_d, part_c), mask2), tmp_1), bf128_one());
      zk_hash_128_update(a1_ctx + 1, share_1);
    }
  }
}

static void aes_prove_128_masked(const uint8_t* w_share, vbb_t* vbb, const uint8_t* in_share, const uint8_t* out_share,
                          const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
                          const faest_paramset_t* params) {
  uint8_t* k_share = alloca(2 * (FAEST_128F_R + 1) * 128 / 8);

  zk_hash_128_ctx a0_ctx_share[2];
  zk_hash_128_ctx a1_ctx_share[2];

  zk_hash_128_init(&a0_ctx_share[0], chall);
  zk_hash_128_init(&a1_ctx_share[0], chall);
  zk_hash_128_init(&a0_ctx_share[1], chall);
  zk_hash_128_init(&a1_ctx_share[1], chall);

  aes_key_schedule_128_masked(w_share, vbb, &a0_ctx_share[0], &a1_ctx_share[0], k_share, params);

  unsigned int offset = FAEST_128F_Lke;
  aes_enc_constraints_128_masked(in_share, out_share, w_share, vbb, offset, k_share, &a0_ctx_share[0], &a1_ctx_share[0]);

  uint8_t a_tilde_share[2][FAEST_128F_LAMBDA / 8];
  uint8_t b_tilde_share[2][FAEST_128F_LAMBDA / 8];

  zk_hash_128_finalize(a_tilde_share[0], &a1_ctx_share[0], bf128_load(get_vole_u_share(vbb, 0) + FAEST_128F_L / 8));
  zk_hash_128_finalize(b_tilde_share[0], &a0_ctx_share[0], bf128_sum_poly_vbb_share(vbb, FAEST_128F_L, 0));

  zk_hash_128_finalize(a_tilde_share[1], &a1_ctx_share[1], bf128_load(get_vole_u_share(vbb, 1) + FAEST_128F_L / 8));
  zk_hash_128_finalize(b_tilde_share[1], &a0_ctx_share[1], bf128_sum_poly_vbb_share(vbb, FAEST_128F_L, 1));

  // LEAKS PUBLIC VALUES...
  for(int i = 0; i < FAEST_128F_LAMBDA / 8; i++){
    a_tilde[i] = a_tilde_share[0][i] ^ a_tilde_share[1][i];
    b_tilde[i] = b_tilde_share[0][i] ^ b_tilde_share[1][i];
  }
}

void aes_prove_masked(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
               const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params) {
  aes_prove_128_masked(w, vbb, in, out, chall, a_tilde, b_tilde, params);
}