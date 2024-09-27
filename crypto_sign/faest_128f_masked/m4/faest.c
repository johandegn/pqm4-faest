/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest.h"
#include "aes.h"
#include "faest_aes.h"
#include "randomness.h"
#include "random_oracle.h"
#include "utils.h"
#include "vole.h"
#include "universal_hashing.h"
#include "vbb.h"

// helpers to compute position in signature (sign)

ATTR_PURE static inline uint8_t* signature_c(uint8_t* base_ptr, unsigned int index,
                                             const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + index * ell_hat_bytes;
}

ATTR_PURE static inline uint8_t* signature_u_tilde(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes;
}

ATTR_PURE static inline uint8_t* signature_d(uint8_t* base_ptr, const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes;
}

ATTR_PURE static inline uint8_t* signature_a_tilde(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes;
}

ATTR_PURE static inline uint8_t* signature_pdec(uint8_t* base_ptr, unsigned int index,
                                                const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.t0;
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr + index * (params->faest_param.k0 + 2) * lambda_bytes;
  } else {
    return base_ptr +
           ((index - tau0) * (params->faest_param.k1 + 2) + tau0 * (params->faest_param.k0 + 2)) *
               lambda_bytes;
  }
}

ATTR_PURE static inline uint8_t* signature_com(uint8_t* base_ptr, unsigned int index,
                                               const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.t0;
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr +
           (index * (params->faest_param.k0 + 2) + params->faest_param.k0) * lambda_bytes;
  } else {
    return base_ptr + ((index - tau0) * (params->faest_param.k1 + 2) + params->faest_param.k1 +
                       tau0 * (params->faest_param.k0 + 2)) *
                          lambda_bytes;
  }
}

ATTR_PURE static inline uint8_t* signature_chall_3(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const size_t lambda_bytes = params->faest_param.lambda / 8;
  return base_ptr + params->faest_param.sigSize - IV_SIZE - lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_iv(uint8_t* base_ptr, const faest_paramset_t* params) {
  return base_ptr + params->faest_param.sigSize - IV_SIZE;
}

// helpers to compute position in signature (verify)

ATTR_PURE inline const uint8_t* dsignature_c(const uint8_t* base_ptr, unsigned int index,
                                             const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + index * ell_hat_bytes;
}

ATTR_PURE inline const uint8_t* dsignature_u_tilde(const uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes;
}

ATTR_PURE inline const uint8_t* dsignature_d(const uint8_t* base_ptr,
                                             const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_a_tilde(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes;
}

ATTR_PURE inline const uint8_t* dsignature_pdec(const uint8_t* base_ptr, unsigned int index,
                                                const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.t0;
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr + index * (params->faest_param.k0 + 2) * lambda_bytes;
  } else {
    return base_ptr +
           ((index - tau0) * (params->faest_param.k1 + 2) + tau0 * (params->faest_param.k0 + 2)) *
               lambda_bytes;
  }
}

ATTR_PURE inline const uint8_t* dsignature_com(const uint8_t* base_ptr, unsigned int index,
                                               const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.t0;
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr +
           (index * (params->faest_param.k0 + 2) + params->faest_param.k0) * lambda_bytes;
  } else {
    return base_ptr + ((index - tau0) * (params->faest_param.k1 + 2) + params->faest_param.k1 +
                       tau0 * (params->faest_param.k0 + 2)) *
                          lambda_bytes;
  }
}

ATTR_PURE inline const uint8_t* dsignature_chall_3(const uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const size_t lambda_bytes = params->faest_param.lambda / 8;
  return base_ptr + params->faest_param.sigSize - IV_SIZE - lambda_bytes;
}

ATTR_PURE inline const uint8_t* dsignature_iv(const uint8_t* base_ptr,
                                              const faest_paramset_t* params) {
  return base_ptr + params->faest_param.sigSize - IV_SIZE;
}

static void hash_mu(uint8_t* mu, const uint8_t* owf_input, const uint8_t* owf_output,
                    size_t owf_size, const uint8_t* msg, size_t msglen, unsigned int lambda) {
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  H1_update(&h1_ctx, owf_input, owf_size);
  H1_update(&h1_ctx, owf_output, owf_size);
  H1_update(&h1_ctx, msg, msglen);
  H1_final(&h1_ctx, mu, 2 * lambda / 8);
}

static void hash_challenge_1(uint8_t* chall_1, const uint8_t* mu, const uint8_t* hcom,
                             const uint8_t* c, const uint8_t* iv, unsigned int lambda,
                             unsigned int ell, unsigned int tau) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int ell_hat_bytes = ell / 8 + lambda_bytes * 2 + UNIVERSAL_HASH_B;

  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);
  H2_update(&h2_ctx, mu, lambda_bytes * 2);
  H2_update(&h2_ctx, hcom, lambda_bytes * 2);
  H2_update(&h2_ctx, c, ell_hat_bytes * (tau - 1));
  H2_update(&h2_ctx, iv, IV_SIZE);
  H2_final(&h2_ctx, chall_1, 5 * lambda_bytes + 8);
}

static void hash_challenge_2(uint8_t* chall_2, const uint8_t* chall_1, const uint8_t* u_tilde,
                             const uint8_t* h_v, const uint8_t* d, unsigned int lambda,
                             unsigned int ell) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int ell_bytes     = ell / 8;
  const unsigned int u_tilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;

  H2_context_t h2_ctx_1;
  H2_init(&h2_ctx_1, lambda);
  H2_update(&h2_ctx_1, chall_1, 5 * lambda_bytes + 8);
  H2_update(&h2_ctx_1, u_tilde, u_tilde_bytes);
  H2_update(&h2_ctx_1, h_v, 2 * lambda_bytes);
  H2_update(&h2_ctx_1, d, ell_bytes);
  H2_final(&h2_ctx_1, chall_2, 3 * lambda_bytes + 8);
}

static void hash_challenge_3(uint8_t* chall_3, const uint8_t* chall_2, const uint8_t* a_tilde,
                             const uint8_t* b_tilde, unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  H2_context_t h2_ctx_2;
  H2_init(&h2_ctx_2, lambda);
  H2_update(&h2_ctx_2, chall_2, 3 * lambda_bytes + 8);
  H2_update(&h2_ctx_2, a_tilde, lambda_bytes);
  H2_update(&h2_ctx_2, b_tilde, lambda_bytes);
  H2_final(&h2_ctx_2, chall_3, lambda_bytes);
}

void faest_sign(uint8_t* sig, const uint8_t* msg, size_t msglen, const uint8_t* owf_key,
                const uint8_t* owf_input, const uint8_t* owf_output, const uint8_t* rho,
                size_t rholen, const faest_paramset_t* params) {
  const unsigned int l           = params->faest_param.l;
  const unsigned int ell_bytes   = l / 8;
  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int lambdaBytes = lambda / 8;
  const unsigned int tau         = params->faest_param.tau;
  const unsigned int tau0        = params->faest_param.t0;
  const unsigned int ell_hat     = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  // const unsigned int ell_hat_bytes = ell_hat / 8;

  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, owf_output, params->faest_param.pkSize / 2, msg, msglen, lambda);

  uint8_t rootkey[MAX_LAMBDA_BYTES];
  uint8_t key_share[2][MAX_LAMBDA_BYTES]       = {0};
  uint8_t owf_input_share[2][MAX_LAMBDA_BYTES] = {0};
  uint8_t owf_output_share[2][16];

  {
#ifdef KECCAK_MASK_NONE
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, lambda);
    H3_update(&h3_ctx, owf_key, lambdaBytes);
    H3_update(&h3_ctx, mu, lambdaBytes * 2);
    if (rho && rholen) {
      H3_update(&h3_ctx, rho, rholen);
    }
    H3_final(&h3_ctx, rootkey, lambdaBytes, signature_iv(sig, params));
#else
    uint8_t rootkey_share[MAX_LAMBDA_BYTES]  = {0};
    uint8_t owf_key_share0[MAX_LAMBDA_BYTES] = {0};
    uint8_t owf_key_share1[MAX_LAMBDA_BYTES] = {0};
    uint8_t mu_share[2 * MAX_LAMBDA_BYTES]   = {0};
    uint8_t rho_share[MAX_LAMBDA_BYTES]      = {0};
    uint8_t iv_share[16]                     = {0};
    uint8_t* rootkey_shares[2]               = {rootkey, rootkey_share};
    uint8_t* iv_shares[2]                    = {signature_iv(sig, params), iv_share};
    const uint8_t* owf_key_shares[2]         = {owf_key_share0, owf_key_share1};
    const uint8_t* mu_shares[2]              = {mu, mu_share};
    const uint8_t* rho_shares[2]             = {rho, rho_share};
    H3_context_t h3_ctx;
    rand_mask(owf_key_share0, lambdaBytes);
    for (size_t i = 0; i < lambdaBytes; i++) {
      owf_key_share1[i] = owf_key[i] ^ owf_key_share0[i];
    }

    for (int i = 0; i < MAX_LAMBDA_BYTES; i++) {
      rand_mask(&key_share[0][i], 1);
      key_share[1][i] = owf_key[i] ^ key_share[0][i];
    }
    for (int i = 0; i < MAX_LAMBDA_BYTES; i++) {
      rand_mask(&owf_input_share[0][i], 1);
      owf_input_share[1][i] = owf_input[i] ^ owf_input_share[0][i];
    }
    for (int i = 0; i < 16; i++) {
      rand_mask(&owf_output_share[0][i], 1);
      owf_output_share[1][i] = owf_output[i] ^ owf_output_share[0][i];
    }

    H3_init(&h3_ctx, lambda);
    H3_update(&h3_ctx, owf_key_shares, lambdaBytes);
    H3_update(&h3_ctx, mu_shares, lambdaBytes * 2);
    if (rho && rholen) {
      H3_update(&h3_ctx, rho_shares, rholen);
    }
    H3_final(&h3_ctx, rootkey_shares, lambdaBytes, iv_shares);
    for (size_t i = 0; i < lambdaBytes; i++) {
      rootkey[i] ^= rootkey_share[i];
    }
    for (size_t i = 0; i < 16; i++) {
      signature_iv(sig, params)[i] ^= iv_share[i];
    }
#endif
  }

  vbb_t vbb;
  // TODO: find a solution for setting argument (dynamic or static)?
  const unsigned int len = ell_hat;
  uint8_t* hcom          = alloca(MAX_LAMBDA_BYTES * 2);
  uint8_t* u             = alloca(ell_hat / 8);
  uint8_t* v_cache       = alloca(len * lambdaBytes);
  uint8_t* v_buf         = alloca(lambdaBytes);
  uint8_t* vk_buf        = NULL;
  uint8_t* vk_cache      = NULL;
  if (!(params->faest_paramid > 6)) {
    vk_buf   = alloca(lambdaBytes);
    vk_cache = alloca(params->faest_param.Lke * lambdaBytes);
  }
  init_stack_allocations_sign(&vbb, hcom, u, v_cache, v_buf, vk_buf, vk_cache);
  init_vbb_sign(&vbb, len, rootkey, signature_iv(sig, params), signature_c(sig, 0, params), params);

  uint8_t chall_1[(5 * MAX_LAMBDA_BYTES) + 8];
  hash_challenge_1(chall_1, mu, get_com_hash(&vbb), signature_c(sig, 0, params),
                   signature_iv(sig, params), lambda, l, tau);

  vole_hash(signature_u_tilde(sig, params), chall_1, get_vole_u(&vbb), l, lambda);

  prepare_hash_sign(&vbb);
  uint8_t h_v[MAX_LAMBDA_BYTES * 2];
  {
    H1_context_t h1_ctx_1;
    H1_init(&h1_ctx_1, lambda);

    uint8_t V_tilde[MAX_LAMBDA_BYTES + UNIVERSAL_HASH_B];
    for (unsigned int i = 0; i != lambda; ++i) {
      vole_hash(V_tilde, chall_1, get_vole_v_hash(&vbb, i), l, lambda);
      H1_update(&h1_ctx_1, V_tilde, lambdaBytes + UNIVERSAL_HASH_B);
    }
    H1_final(&h1_ctx_1, h_v, lambdaBytes * 2);
  }

#define WITNESS_MASKING
#ifdef WITNESS_MASKING
  // secret sharing the key and then computing the extended witness
  uint8_t* w_share = alloca(2 * (l + 7) / 8);

  // secret sharing the input, (Not needed only to remove false positive leakage)

  w_share = aes_extend_witness_masked(&key_share[0][0], &owf_input_share[0][0], params, w_share);

  xor_u8_array(w_share, get_vole_u(&vbb), signature_d(sig, params), ell_bytes);
  xor_u8_array(w_share + (l + 7) / 8, signature_d(sig, params), signature_d(sig, params),
               ell_bytes);
#else
  uint8_t* w = alloca((l + 7) / 8);
  w          = aes_extend_witness(owf_key, owf_input, params, w);
  xor_u8_array(w, get_vole_u(&vbb), signature_d(sig, params), ell_bytes);
#endif

  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_2(chall_2, chall_1, signature_u_tilde(sig, params), h_v, signature_d(sig, params),
                   lambda, l);

  prepare_aes_sign(&vbb);
  uint8_t b_tilde[MAX_LAMBDA_BYTES];

#ifdef WITNESS_MASKING
  if (params->faest_paramid == 1 || params->faest_paramid == 2) {
    uint8_t* vk_mask = alloca(params->faest_param.Lke * lambdaBytes);
    uint8_t* v_mask = alloca(ell_hat * lambdaBytes);
    uint8_t* u_mask = alloca(ell_hat/8);
    setup_mask_storage(&vbb, vk_mask, v_mask, u_mask);
    aes_prove_masked(w_share, &vbb, &owf_input_share[0][0], &owf_output_share[0][0], chall_2,
                     signature_a_tilde(sig, params), b_tilde, params);
  } else {
    uint8_t* w = alloca((l + 7) / 8);
    for (unsigned int i = 0; i < (l + 7) / 8; i++) {
      w[i] = w_share[i] ^ w_share[i + (l + 7) / 8];
    }
    aes_prove(w, &vbb, owf_input, owf_output, chall_2, signature_a_tilde(sig, params), b_tilde,
              params);
  }
#else
  aes_prove(w, &vbb, owf_input, owf_output, chall_2, signature_a_tilde(sig, params), b_tilde,
            params);
#endif

  // free(w);
  // w = NULL;

  hash_challenge_3(signature_chall_3(sig, params), chall_2, signature_a_tilde(sig, params), b_tilde,
                   lambda);

  for (unsigned int i = 0; i < tau; i++) {
    uint8_t s_[MAX_DEPTH];
    ChalDec(signature_chall_3(sig, params), i, params->faest_param.k0, params->faest_param.t0,
            params->faest_param.k1, params->faest_param.t1, s_);
    const unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    vector_open_ondemand(&vbb, i, s_, signature_pdec(sig, i, params), signature_com(sig, i, params),
                         depth);
  }
  // clean_vbb(&vbb);
}

int faest_verify(const uint8_t* msg, size_t msglen, const uint8_t* sig, const uint8_t* owf_input,
                 const uint8_t* owf_output, const faest_paramset_t* params) {
  const unsigned int l           = params->faest_param.l;
  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int lambdaBytes = lambda / 8;
  const unsigned int tau         = params->faest_param.tau;
  const unsigned int ell_hat     = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;

  vbb_t vbb;
  const unsigned int len = ell_hat;
  uint8_t* hcom          = alloca(MAX_LAMBDA_BYTES * 2);
  uint8_t* q_cache       = alloca(len * lambdaBytes);
  uint8_t* Dtilde_buf    = alloca(lambdaBytes + UNIVERSAL_HASH_B);
  uint8_t* v_buf         = alloca(lambdaBytes);
  uint8_t* vk_buf        = NULL;
  uint8_t* vk_cache      = NULL;
  if (!(params->faest_paramid > 6)) {
    vk_buf   = alloca(lambdaBytes);
    vk_cache = alloca(params->faest_param.Lke * lambdaBytes);
  }
  init_stack_allocations_verify(&vbb, hcom, q_cache, Dtilde_buf, v_buf, vk_buf, vk_cache);
  init_vbb_verify(&vbb, len, params, sig);

  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, owf_output, params->faest_param.pkSize / 2, msg, msglen, lambda);

  uint8_t chall_1[(5 * MAX_LAMBDA_BYTES) + 8];
  hash_challenge_1(chall_1, mu, vbb.com_hash, dsignature_c(sig, 0, params),
                   dsignature_iv(sig, params), lambda, l, tau);

  prepare_hash_verify(&vbb);
  uint8_t h_v[MAX_LAMBDA_BYTES * 2];
  {
    H1_context_t h1_ctx_1;
    H1_init(&h1_ctx_1, lambda);

    uint8_t Q_tilde[MAX_LAMBDA_BYTES + UNIVERSAL_HASH_B];
    for (unsigned int i = 0; i != lambda; ++i) {
      vole_hash(Q_tilde, chall_1, get_vole_q_hash(&vbb, i), l, lambda);
      xor_u8_array(Q_tilde, get_dtilde(&vbb, i), Q_tilde, lambdaBytes + UNIVERSAL_HASH_B);
      H1_update(&h1_ctx_1, Q_tilde, lambdaBytes + UNIVERSAL_HASH_B);
    }
    H1_final(&h1_ctx_1, h_v, lambdaBytes * 2);
  }

  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_2(chall_2, chall_1, dsignature_u_tilde(sig, params), h_v,
                   dsignature_d(sig, params), lambda, l);

  prepare_aes_verify(&vbb);
  uint8_t* q_tilde = alloca(lambdaBytes);
  uint8_t* b_tilde =
      aes_verify(&vbb, chall_2, dsignature_chall_3(sig, params), dsignature_a_tilde(sig, params),
                 owf_input, owf_output, params, q_tilde);

  uint8_t chall_3[MAX_LAMBDA_BYTES];
  hash_challenge_3(chall_3, chall_2, dsignature_a_tilde(sig, params), b_tilde, lambda);
  // free(b_tilde);
  // b_tilde = NULL;
  // clean_vbb(&vbb);

  return memcmp(chall_3, dsignature_chall_3(sig, params), lambdaBytes) == 0 ? 0 : -1;
}
