#include "faest.h"
#include "fields.h"
#include "vole.h"
#include "universal_hashing.h"
#include "utils.h"
#include "parameters.h"


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

void aes_key_schedule_backward_128_vbb_vk_round_share(vbb_t* vbb, bf128_t* bf_out, unsigned int j, unsigned int share);

static void __attribute__ ((noinline)) aes_key_schedule_128_masked(const uint8_t* w_share, vbb_t* vbb,
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

    aes_key_schedule_backward_1_round_share(w_share + FAEST_128F_LAMBDA / 8, k, &w_dash[0][0], j, params, false);
    aes_key_schedule_backward_128_vbb_vk_round_share(vbb, &v_w_dash[0][0], j, 0);
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 10..11
      bf_k_hat_share[0][(r + 3) % 4]   = bf128_byte_combine_bits_sclf(k[(96 + iwd*8 + 8 * r) / 8]);
      bf_v_k_hat_share[0][(r + 3) % 4] = bf128_byte_combine_vk_share(vbb, (96 + iwd*8 + 8 * r), 0);
      bf_w_dash_hat_share[0][r]        = bf128_byte_combine_bits_sclf(w_dash[0][(8 * r) / 8]);
      bf_v_w_dash_hat_share[0][r]      = bf128_byte_combine_sclf(v_w_dash[0] + (8 * r));
    }

    aes_key_schedule_backward_1_round_share(w_share + FAEST_128F_LAMBDA / 8 + FAEST_128F_L/8, k + (FAEST_128F_R + 1) * 128 / 8, &w_dash[1][0], j, params, true);
    aes_key_schedule_backward_128_vbb_vk_round_share(vbb, &v_w_dash[1][0], j, 1);
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 10..11
      bf_k_hat_share[1][(r + 3) % 4]   = bf128_byte_combine_bits_sclf((k + (FAEST_128F_R + 1) * 128 / 8)[(96 + iwd*8 + 8 * r) / 8]);
      bf_v_k_hat_share[1][(r + 3) % 4] = bf128_byte_combine_vk_share(vbb, (96 + iwd*8 + 8 * r), 1);
      bf_w_dash_hat_share[1][r]        = bf128_byte_combine_bits_sclf(w_dash[1][(8 * r) / 8]);
      bf_v_w_dash_hat_share[1][r]      = bf128_byte_combine_sclf(v_w_dash[1] + (8 * r));
    }
    for (unsigned int r = 0; r <= 3; r++) {
      bf128_t part_a, part_b, part_c, part_d;
      bf128_add_wrapper(&part_a, bf_v_k_hat_share[0] + r, bf_k_hat_share[0] + r);
      bf128_add_wrapper(&part_b, bf_w_dash_hat_share[0] + r, bf_v_w_dash_hat_share[0] + r);
      bf128_add_wrapper(&part_d, bf_k_hat_share[1] + r, bf_v_k_hat_share[1] + r);
      bf128_add_wrapper(&part_c, bf_v_w_dash_hat_share[1] + r, bf_w_dash_hat_share[1] + r);

      bf128_t tmp_x;
      // instead of storing in A0, A1, hash it
      bf128_t mask1 = bf128_rand();
      bf128_t mask2 = bf128_rand();

      bf128_t tmp_0;
      bf128_mul_wrapper(&tmp_x, bf_v_w_dash_hat_share[1] + r, bf_v_k_hat_share[0] + r);
      bf128_add_wrapper(&tmp_0, &tmp_x, &mask1);
      bf128_mul_wrapper(&tmp_x, bf_v_k_hat_share[0] + r, bf_v_w_dash_hat_share[0] + r);
      bf128_add_wrapper(&tmp_0, &tmp_0, &tmp_x);
      bf128_mul_wrapper(&tmp_x, bf_v_w_dash_hat_share[0] + r, bf_v_k_hat_share[1] + r);
      bf128_add_wrapper(&tmp_0, &tmp_0, &tmp_x);
      zk_hash_128_update(a0_ctx, tmp_0);

      bf128_t share_0;
      bf128_mul_wrapper(&tmp_x, &part_a, &part_c);
      bf128_add_wrapper(&share_0, &tmp_x, &mask2);
      bf128_mul_wrapper(&tmp_x, &part_a, &part_b);
      bf128_add_wrapper(&share_0, &share_0, &tmp_x);
      bf128_add_wrapper(&share_0, &share_0, &tmp_0);
      bf128_mul_wrapper(&tmp_x, &part_d, &part_b);
      bf128_add_wrapper(&share_0, &share_0, &tmp_x);
      zk_hash_128_update(a1_ctx, share_0);

      bf128_t tmp_1;
      bf128_mul_wrapper(&tmp_x, bf_v_k_hat_share[1] + r, bf_v_w_dash_hat_share[1] + r);
      bf128_add_wrapper(&tmp_1, &tmp_x, &mask1);
      zk_hash_128_update(a0_ctx + 1, tmp_1);

      bf128_t share_1;
      bf128_mul_wrapper(&tmp_x, &part_d, &part_c);
      bf128_add_wrapper(&share_1, &tmp_x, &mask2);
      bf128_add_wrapper(&share_1, &share_1, &tmp_1);
      tmp_x = bf128_one();
      bf128_add_wrapper(&share_1, &share_1, &tmp_x);
      zk_hash_128_update(a1_ctx + 1, share_1); 
    }
    iwd += 128 / 8;
  }
}

static void aes_enc_forward_128_1_round(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf128_t* bf_y, int round) {
  if (round == 0){
    for (unsigned int i = 0; i < 16; i++) {
      const uint8_t xin = in[i];
      bf_y[i] = bf128_add(bf128_byte_combine_bits_sclf(xin), bf128_byte_combine_bits_sclf(xk[i]));
    }
  }

  if (round > 0){
    const bf128_t bf_two   = bf128_byte_combine_bits_sclf(2);
    const bf128_t bf_three = bf128_byte_combine_bits_sclf(3);
    unsigned int j = round;
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf128_byte_combine_bits_sclf(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf128_byte_combine_bits_sclf(xk[(ik + 8 * r) / 8]);
      }

      // Step : 14
      bf_y[iy + 0] = bf128_add(bf_xk_hat[0], bf128_mul_sclf(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul_sclf(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul_sclf(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul_sclf(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul_sclf(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul_sclf(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf128_add(bf_xk_hat[3], bf128_mul_sclf(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul_sclf(bf_x_hat[3], bf_two));
    }
  }
}

static void aes_enc_backward_128_1_round_share(const uint8_t* x, const uint8_t* xk, const uint8_t* out,
                                   bf128_t* y_out, unsigned int round, unsigned int share) {
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
      const uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2) ^ (share * 0x5);

      // Step: 18
      y_out[4 * c + r] = bf128_byte_combine_bits_sclf(ytilde);
    }
  }
}

void aes_enc_forward_backward_128_share(vbb_t* vbb, unsigned int offset, const uint8_t* in, const uint8_t* out,
                                         bf128_t* vs, bf128_t* vs_old, bf128_t* vs_dash, unsigned int round, unsigned int share);


void __attribute__ ((noinline)) aes_enc_hash_loop( bf128_t s_share[2][16], bf128_t vs_share[2][16], bf128_t s_dash_share[2][16], bf128_t vs_dash_share[2][16], zk_hash_128_ctx* a0_ctx, zk_hash_128_ctx* a1_ctx, unsigned int j){
  bf128_t part_a, part_b, part_c, part_d;
  bf128_add_wrapper(&part_a, vs_share[0] + j, s_share[0] + j);
  bf128_add_wrapper(&part_b, s_dash_share[0] + j, vs_dash_share[0] + j);
  bf128_add_wrapper(&part_d, s_share[1] + j, vs_share[1] + j);
  bf128_add_wrapper(&part_c, vs_dash_share[1] + j, s_dash_share[1] + j);

  bf128_t tmp_x;
  bf128_t mask1 = bf128_rand();

  bf128_t tmp_0;
  bf128_mul_wrapper(&tmp_x, vs_dash_share[1] + j, vs_share[0] + j);
  bf128_add_wrapper(&tmp_0, &tmp_x, &mask1);
  bf128_mul_wrapper(&tmp_x, vs_share[0] + j, vs_dash_share[0] + j);
  bf128_add_wrapper(&tmp_0, &tmp_0, &tmp_x);
  bf128_mul_wrapper(&tmp_x, vs_dash_share[0] + j, vs_share[1] + j);
  bf128_add_wrapper(&tmp_0, &tmp_0, &tmp_x);
  zk_hash_128_update(a0_ctx, tmp_0);

  bf128_t mask2 = bf128_rand();
  bf128_t share_0;
  bf128_mul_wrapper(&tmp_x, &part_a, &part_c);
  bf128_add_wrapper(&share_0, &tmp_x, &mask2);
  bf128_mul_wrapper(&tmp_x, &part_a, &part_b);
  bf128_add_wrapper(&share_0, &share_0, &tmp_x);
  bf128_add_wrapper(&share_0, &share_0, &tmp_0);
  bf128_mul_wrapper(&tmp_x, &part_d, &part_b);
  bf128_add_wrapper(&share_0, &share_0, &tmp_x);
  zk_hash_128_update(a1_ctx, share_0);

  bf128_t tmp_1;
  bf128_mul_wrapper(&tmp_x, vs_share[1] + j, vs_dash_share[1] + j);
  bf128_add_wrapper(&tmp_1, &tmp_x, &mask1);
  zk_hash_128_update(a0_ctx + 1, tmp_1);

  bf128_t share_1;
  bf128_mul_wrapper(&tmp_x, &part_d, &part_c);
  bf128_add_wrapper(&share_1, &tmp_x, &mask2);
  bf128_add_wrapper(&share_1, &share_1, &tmp_1);
  tmp_x = bf128_one();
  bf128_add_wrapper(&share_1, &share_1, &tmp_x);
  zk_hash_128_update(a1_ctx + 1, share_1);
}

static void  __attribute__ ((noinline)) aes_enc_constraints_128_masked(const uint8_t* in_share, const uint8_t* out_share, const uint8_t* w_share,
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
    aes_enc_forward_backward_128_share(vbb, offset, in_share, out_share, vs_share[0], vs_share_old[0], vs_dash_share[0], i, 0);
    
    aes_enc_forward_128_1_round(w_share + + FAEST_128F_L/8, k_share + (FAEST_128F_R + 1) * 128 / 8, in_share + MAX_LAMBDA_BYTES, s_share[1], i);
    aes_enc_backward_128_1_round_share(w_share + + FAEST_128F_L/8, k_share + (FAEST_128F_R + 1) * 128 / 8, out_share + 16, s_dash_share[1], i, 1);
    aes_enc_forward_backward_128_share(vbb, offset, in_share + MAX_LAMBDA_BYTES, out_share + MAX_LAMBDA_BYTES, vs_share[1], vs_share_old[1], vs_dash_share[1], i, 1);

    for (unsigned int j = 0; j < 16; j++){
      aes_enc_hash_loop(s_share, vs_share, s_dash_share, vs_dash_share, a0_ctx, a1_ctx, j);
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