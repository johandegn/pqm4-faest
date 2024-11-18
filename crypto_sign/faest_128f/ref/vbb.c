#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>

#include "macros.h"
#include "vbb.h"
#include "vole.h"
#include "vc.h"
#include "instances.h"
#include "faest.h"
#include "faest_aes.h"
#include "fields.h"
#include "parameters.h"

static void setup_vk_cache(vbb_t* vbb);

ATTR_CONST ATTR_ALWAYS_INLINE static inline bool is_em_variant(faest_paramid_t id) {
  return id > 6;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bool is_column_cached(vbb_t* vbb, unsigned int index) {
  bool above_cache_start = index >= vbb->cache_idx;
  bool below_cache_end   = index < vbb->cache_idx + vbb->column_count;
  return above_cache_start && below_cache_end;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bool is_row_cached(vbb_t* vbb, unsigned int index) {
  bool above_cache_start = index >= vbb->cache_idx;
  bool below_cache_end   = index < vbb->cache_idx + vbb->row_count;
  return above_cache_start && below_cache_end;
}

static void recompute_hash_sign(vbb_t* vbb, unsigned int start, unsigned int end) {
  const unsigned int lambda = vbb->params->faest_param.lambda;
  const unsigned int ell    = vbb->params->faest_param.l;
  const unsigned int ellhat = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int capped_end   = MIN(end, lambda);

  partial_vole_commit_column(vbb->root_key, vbb->iv, ellhat, start, capped_end,
                             vole_mode_v(vbb->vole_cache), vbb->params);
  vbb->cache_idx = start;
}

static void recompute_vole_row(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda = vbb->params->faest_param.lambda;
  const unsigned int ell    = vbb->params->faest_param.l;
  const unsigned int ellhat = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;

  if (start >= ell) {
    start = start - len + 1;
  }
  if (len >= ell + lambda) {
    start = 0;
  } else if (start + len > ell + lambda) {
    start = ell + lambda - len;
  }

  partial_vole_commit_row(vbb->root_key, vbb->iv, ellhat, start, start + len, vbb->params,
                          vbb->vole_cache);
  vbb->cache_idx = start;
}

// len is the number of OLE v's that is allowed to be stored in memory.
// Hence we store (at most) len*lambda in memory.
void init_vbb_sign(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv,
                   uint8_t* c, const faest_paramset_t* params) {
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat       = params->faest_param.l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int row_count    = len;
  const unsigned int column_count = (size_t)row_count * (size_t)lambda_bytes / (size_t)ellhat_bytes;
  assert(column_count >= 1);

  vbb->party        = SIGNER;
  vbb->iv           = iv;
  vbb->params       = params;
  vbb->root_key     = root_key;
  vbb->full_size    = len >= ellhat;
  vbb->row_count    = row_count;
  vbb->column_count = column_count;

  sign_vole_mode_ctx_t mode =
      vbb->full_size ? vole_mode_all_sign(vbb->vole_cache, vbb->vole_U, vbb->com_hash, c)
                     : vole_mode_u_hcom_c(vbb->vole_U, vbb->com_hash, c);

  partial_vole_commit_column(vbb->root_key, vbb->iv, ellhat, 0, lambda, mode, vbb->params);
}

void init_stack_allocations_sign(vbb_t* vbb, uint8_t* hcom, uint8_t* u, uint8_t* v,
                                 uint8_t* v_buffer, uint8_t* vk_buffer, uint8_t* vk_cache) {
  vbb->com_hash   = hcom;
  vbb->vole_U     = u;
  vbb->vole_cache = v;
  vbb->v_buf      = v_buffer;
  vbb->vk_buf     = vk_buffer;
  vbb->vk_cache   = vk_cache;
}

void prepare_hash_sign(vbb_t* vbb) {
  if (vbb->full_size) {
    vbb->cache_idx = 0;
    return;
  }
  recompute_hash_sign(vbb, 0, vbb->column_count);
}

void prepare_aes_sign(vbb_t* vbb) {
  if (vbb->full_size) {
    vbb->cache_idx = 0;
  } else {
    recompute_vole_row(vbb, 0, vbb->row_count);
  }
  if (!is_em_variant(vbb->params->faest_paramid)) {
    setup_vk_cache(vbb);
  }
}

void vector_open_ondemand(vbb_t* vbb, unsigned int idx, const uint8_t* s_, uint8_t* sig_pdec,
                          uint8_t* sig_com, unsigned int depth) {
  const unsigned int lambda       = vbb->params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int tau          = vbb->params->faest_param.tau;
  uint8_t* expanded_keys          = alloca(tau * lambda_bytes);
  prg(vbb->root_key, vbb->iv, expanded_keys, lambda, lambda_bytes * tau);

  vec_com_t vec_com;
  vector_commitment(expanded_keys + lambda_bytes * idx, lambda, depth, NULL, &vec_com);
  vector_open(&vec_com, s_, sig_pdec, sig_com, depth, vbb->iv, lambda);
}

static inline void apply_correction_values_column(vbb_t* vbb, unsigned int start,
                                                  unsigned int len) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int tau           = vbb->params->faest_param.tau;
  const unsigned int tau0          = vbb->params->faest_param.t0;
  const unsigned int tau1          = vbb->params->faest_param.t1;
  const unsigned int k0            = vbb->params->faest_param.k0;
  const unsigned int k1            = vbb->params->faest_param.k1;

  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* c      = dsignature_c(vbb->sig, 0, vbb->params);
  unsigned int col_idx  = k0;
  vbb->cache_idx = start;
  for (unsigned int i = 1; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;
    if (start >= col_idx + depth) {
      col_idx += depth;
      continue;
    }
    if (col_idx >= start + len) {
      break;
    }
    uint8_t delta[MAX_DEPTH];
    ChalDec(chall3, i, k0, tau0, k1, tau1, delta);

    for (unsigned int d = 0; d < depth; d++, col_idx++) {
      if (start > col_idx) {
        continue;
      }
      masked_xor_u8_array(
          vbb->vole_cache + (col_idx - start) * ell_hat_bytes, c + (i - 1) * ell_hat_bytes,
          vbb->vole_cache + (col_idx - start) * ell_hat_bytes, delta[d], ell_hat_bytes);

      if (col_idx + 1 >= start + len) {
        return;
      }
    }
  }
}

static void setup_pdec_com(vbb_t* vbb, const uint8_t** pdec, const uint8_t** com) {
  const unsigned int tau = vbb->params->faest_param.tau;
  for (unsigned int i = 0; i < tau; ++i) {
    pdec[i] = dsignature_pdec(vbb->sig, i, vbb->params);
    com[i]  = dsignature_com(vbb->sig, i, vbb->params);
  }
}

// Verifier implementation
static void recompute_hash_verify(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda  = vbb->params->faest_param.lambda;
  const unsigned int ell     = vbb->params->faest_param.l;
  const unsigned int ell_hat = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int amount  = MIN(len, lambda - start);
  const uint8_t* chall3      = dsignature_chall_3(vbb->sig, vbb->params);

  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  setup_pdec_com(vbb, pdec, com);

  partial_vole_reconstruct_column(vbb->iv, chall3, pdec, com, ell_hat, start, amount,
                                  vole_mode_q(vbb->vole_cache), vbb->params);
  apply_correction_values_column(vbb, start, amount);
  vbb->cache_idx = start;
}

static void apply_correction_values_row(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int tau           = vbb->params->faest_param.tau;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int tau0          = vbb->params->faest_param.t0;
  const unsigned int tau1          = vbb->params->faest_param.t1;
  const unsigned int k0            = vbb->params->faest_param.k0;
  const unsigned int k1            = vbb->params->faest_param.k1;
  const uint8_t* chall3            = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* c                 = dsignature_c(vbb->sig, 0, vbb->params);
  const unsigned int len_bytes = (len + 7) / 8;

  uint8_t* c_trunc      = alloca(len_bytes);
  uint8_t* c_trunc_mask = alloca(len_bytes);
  
  unsigned long col_idx  = k0;
  for (unsigned int i = 1; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;

    uint8_t delta[MAX_DEPTH];
    ChalDec(chall3, i, k0, tau0, k1, tau1, delta);
    // Make a c_trunc
    const uint8_t* c_ptr = c + (i - 1) * ell_hat_bytes;

    unsigned int bit_offset = start % 8;
    unsigned int start_byte = start / 8;
    // If aligned, copy over
    if (bit_offset == 0) {
      memcpy(c_trunc, c_ptr + start_byte, len_bytes);
    } else { // If not aligned
      for (unsigned int j = 0; j < len_bytes; j++) {
        c_trunc[j] =
            (c_ptr[start_byte + j] >> bit_offset) | (c_ptr[start_byte + j + 1] << (8 - bit_offset));
      }
      // Get last part
      c_trunc[len_bytes - 1] = (c_ptr[start_byte + len_bytes - 1] >> bit_offset);
      unsigned int rest      = len - (len_bytes - 1) * 8;
      if (rest > 8 - bit_offset) {
        // Get extra part
        c_trunc[len_bytes - 1] |= c_ptr[start_byte + len_bytes] << (8 - bit_offset);
      }
    }
    // Clear final bits
    unsigned int bits_to_clear = (8 - (len % 8)) % 8;
    c_trunc[len_bytes - 1] &= (uint8_t)0xFF >> bits_to_clear;

    for (unsigned int d = 0; d < depth; d++) {
      // TODO: if delta[d] == 0 continue; else mask = 0xFF
      uint8_t mask = -(delta[d] & 1);
      for(unsigned int k = 0; k < len_bytes; k++) {
        c_trunc_mask[k] = c_trunc[k] & mask;
      }

      uint8_t* q = vbb->vole_cache;
      
      // Shift offset and XOR into v
      unsigned long row_bit_index   = (unsigned long)len * (col_idx);
      unsigned long row_bit_offset  = row_bit_index % 8;
      unsigned long row_byte_offset = row_bit_index / 8;
      // Apply first byte
      q[row_byte_offset] ^= c_trunc_mask[0] << row_bit_offset;
      // Apply the remaining
      for (unsigned int k = 1; k < len_bytes; k++) {
        q[row_byte_offset + k] ^=
            (c_trunc_mask[k - 1] >> (8 - row_bit_offset)) | (c_trunc_mask[k] << row_bit_offset);
      }

      if (row_bit_offset != 0) {
        q[row_byte_offset + len_bytes] ^= c_trunc_mask[len_bytes - 1] >> (8 - row_bit_offset);
      }
      col_idx++;
    }

  }
  free(c_trunc);
  free(c_trunc_mask);

  vbb->cache_idx = start;
}

static void apply_witness_values_row(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int ell          = vbb->params->faest_param.l;
  const unsigned int tau          = vbb->params->faest_param.tau;
  const unsigned int tau0         = vbb->params->faest_param.t0;
  const unsigned int tau1         = vbb->params->faest_param.t1;
  const unsigned int k0           = vbb->params->faest_param.k0;
  const unsigned int k1           = vbb->params->faest_param.k1;
  const uint8_t* d_ptr            = dsignature_d(vbb->sig, vbb->params);
  const uint8_t* chall3           = dsignature_chall_3(vbb->sig, vbb->params);
  unsigned int col_idx            = 0;
  unsigned int effective_len      = MIN(len, ell - start);
  uint8_t delta[MAX_DEPTH];

  if(start >= ell){
    return;
  }
  
  // Make d_trunc
  const unsigned int len_bytes  = (effective_len + 7) / 8;
  uint8_t* d_trunc              = alloca(len_bytes);
  uint8_t* d_trunc_mask         = alloca(len_bytes);

  unsigned int bit_offset = start % 8;
  unsigned int start_byte = start / 8;
  // If aligned, copy over
  if (bit_offset == 0) {
    memcpy(d_trunc, d_ptr + start_byte, len_bytes);
  } else { // If not aligned
    for (unsigned int j = 0; j < len_bytes; j++) {
      d_trunc[j] =
          (d_ptr[start_byte + j] >> bit_offset) | (d_ptr[start_byte + j + 1] << (8 - bit_offset));
    }
    // Get last part
    d_trunc[len_bytes - 1] = (d_ptr[start_byte + len_bytes - 1] >> bit_offset);
    unsigned int rest      = effective_len - (len_bytes - 1) * 8;
    if (rest > 8 - bit_offset) {
      // Get extra part
      d_trunc[len_bytes - 1] |= d_ptr[start_byte + len_bytes] << (8 - bit_offset);
    }
  }
  // Clear final bits
  unsigned int bits_to_clear = (8 - (effective_len % 8)) % 8;
  d_trunc[len_bytes - 1] &= (uint8_t)0xFF >> bits_to_clear;

  for (unsigned int t = 0; t < tau; t++) {
    unsigned int depth = t < tau0 ? k0 : k1;
    ChalDec(chall3, t, k0, tau0, k1, tau1, delta);

    for (unsigned int d = 0; d < depth; d++) {
      // Mul d_truct with delta[d] store in d_trunc_mask
      uint8_t mask = -(delta[d] & 1);
      for(unsigned int k = 0; k < len_bytes; k++) {
        d_trunc_mask[k] = d_trunc[k] & mask;
      }

      // Apply the d_trunc_mask to the vole_cache
      uint8_t* q = vbb->vole_cache;
      
      // Shift offset and XOR into v
      unsigned long row_bit_index   = (unsigned long)len * (col_idx);
      unsigned long row_bit_offset  = row_bit_index % 8;
      unsigned long row_byte_offset = row_bit_index / 8;
      // Apply first byte
      q[row_byte_offset] ^= d_trunc_mask[0] << row_bit_offset;
      // Apply the remaining
      for (unsigned int k = 1; k < len_bytes; k++) {
        q[row_byte_offset + k] ^=
            (d_trunc_mask[k - 1] >> (8 - row_bit_offset)) | (d_trunc_mask[k] << row_bit_offset);
      }

      if (row_bit_offset != 0) {
        q[row_byte_offset + len_bytes] ^= d_trunc_mask[len_bytes - 1] >> (8 - row_bit_offset);
      }
      col_idx++;
    }
  }
  free(d_trunc);
  free(d_trunc_mask);
}

static void apply_witness_values_column(vbb_t* vbb) {
  const unsigned int tau           = vbb->params->faest_param.tau;
  const unsigned int t0            = vbb->params->faest_param.t0;
  const unsigned int k0            = vbb->params->faest_param.k0;
  const unsigned int t1            = vbb->params->faest_param.t1;
  const unsigned int k1            = vbb->params->faest_param.k1;
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int l             = vbb->params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;

  // Apply withness to CMO cache
  unsigned int size = vbb->params->faest_param.l;
  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(dsignature_chall_3(vbb->sig, vbb->params), i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(dsignature_d(vbb->sig, vbb->params), vbb->vole_cache + col * ell_hat_bytes,
                     vbb->vole_cache + col * ell_hat_bytes, (size + 7) / 8);
      }
    }
  }
}

static void recompute_vole_row_reconstruct(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda  = vbb->params->faest_param.lambda;
  const unsigned int ell     = vbb->params->faest_param.l;
  const unsigned int ell_hat = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;

  if (start >= ell) {
    start = start - len + 1;
  }
  if (len >= ell + lambda) {
    start = 0;
  } else if (start + len > ell + lambda) {
    start = ell + lambda - len;
  }

  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  setup_pdec_com(vbb, pdec, com);
  partial_vole_reconstruct_row(vbb->iv, chall3, pdec, com, vbb->vole_cache, ell_hat, vbb->params,
                               start, len);
  apply_correction_values_row(vbb, start, len);
  apply_witness_values_row(vbb, start, len);
  
  vbb->cache_idx = start;
}

void init_vbb_verify(vbb_t* vbb, unsigned int len, const faest_paramset_t* params,
                     const uint8_t* sig) {
  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int l             = params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int row_count     = len;
  const unsigned int column_count =
      (size_t)row_count * (size_t)lambda_bytes / (size_t)ell_hat_bytes;
  assert(column_count >= 1);

  vbb->party        = VERIFIER;
  vbb->params       = params;
  vbb->iv           = dsignature_iv(sig, params);
  vbb->full_size    = len >= ell_hat;
  vbb->sig          = sig;
  vbb->row_count    = row_count;
  vbb->column_count = column_count;

  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  setup_pdec_com(vbb, pdec, com);

  verify_vole_mode_ctx_t vole_mode = (vbb->full_size)
                                         ? vole_mode_all_verify(vbb->vole_cache, vbb->com_hash)
                                         : vole_mode_hcom(vbb->com_hash);
  partial_vole_reconstruct_column(vbb->iv, chall3, pdec, com, ell_hat, 0, lambda, vole_mode,
                                  vbb->params);
  if (vbb->full_size) {
    apply_correction_values_column(vbb, 0, lambda);
  }
}

void init_stack_allocations_verify(vbb_t* vbb, uint8_t* hcom, uint8_t* q, uint8_t* dtilde,
                                   uint8_t* v_buffer, uint8_t* vk_buffer, uint8_t* vk_cache) {
  vbb->com_hash   = hcom;
  vbb->vole_cache = q;
  vbb->Dtilde_buf = dtilde;
  vbb->v_buf      = v_buffer;
  vbb->vk_buf     = vk_buffer;
  vbb->vk_cache   = vk_cache;
}

void prepare_hash_verify(vbb_t* vbb) {
  if (vbb->full_size) {
    vbb->cache_idx = 0;
    return;
  }
  recompute_hash_verify(vbb, 0, vbb->column_count);
}

const uint8_t* get_dtilde(vbb_t* vbb, unsigned int idx) {
  const unsigned int tau0         = vbb->params->faest_param.t0;
  const unsigned int tau1         = vbb->params->faest_param.t1;
  const unsigned int lambda       = vbb->params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int utilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0           = vbb->params->faest_param.k0;
  const unsigned int k1           = vbb->params->faest_param.k1;

  unsigned int t = 0;
  unsigned int j = 0;
  if (idx < k0 * tau0) {
    t = idx / k0;
    j = idx % k0;
  } else {
    t = tau0 + (idx - k0 * tau0) / k1;
    j = (idx - k0 * tau0) % k1;
  }

  uint8_t delta[MAX_DEPTH];
  ChalDec(dsignature_chall_3(vbb->sig, vbb->params), t, k0, tau0, k1, tau1, delta);
  memset(vbb->Dtilde_buf, 0, utilde_bytes);
  masked_xor_u8_array(vbb->Dtilde_buf, dsignature_u_tilde(vbb->sig, vbb->params), vbb->Dtilde_buf,
                      delta[j], utilde_bytes);
  return vbb->Dtilde_buf;
}

void prepare_aes_verify(vbb_t* vbb) {
  if (vbb->full_size) {
    apply_witness_values_column(vbb);
    vbb->cache_idx = 0;
  } else {
    recompute_vole_row_reconstruct(vbb, 0, vbb->row_count);
  }
  if (!is_em_variant(vbb->params->faest_paramid)) {
    setup_vk_cache(vbb);
  }
}

// Get voles for hashing
const uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;

  if (!is_column_cached(vbb, idx)) {
    unsigned int cmo_budget = vbb->column_count;
    recompute_hash_sign(vbb, idx, idx + cmo_budget);
  }
  
  const unsigned int offset = idx - vbb->cache_idx;
  return vbb->vole_cache + offset * ell_hat_bytes;
}

const uint8_t* get_vole_q_hash(vbb_t* vbb, unsigned int idx) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;

  if (!is_column_cached(vbb, idx)) {
    unsigned int cmo_budget = vbb->column_count;
    recompute_hash_verify(vbb, idx, cmo_budget);
  }

  const unsigned int offset = idx - vbb->cache_idx;
  return vbb->vole_cache + offset * ell_hat_bytes;
}

static inline uint8_t* get_vole_row(vbb_t* vbb, unsigned int idx) {
  unsigned int lambda       = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;

  // Check if the idx is within the cache
  if (!is_row_cached(vbb, idx)) {
    unsigned int rmo_budget = vbb->row_count;
    if (vbb->party == VERIFIER) {
      recompute_vole_row_reconstruct(vbb, idx, rmo_budget);
    } else {
      recompute_vole_row(vbb, idx, rmo_budget);
    }
  }

  // Always transpose the VOLE access
  // Compute the new idx based on the starting position of the cache
  unsigned int idx_relative = idx - vbb->cache_idx;
  memset(vbb->v_buf, 0, lambda_bytes);
  // Transpose the VOLE into the buffer
  for (unsigned int column = 0; column != lambda; ++column) {
    ptr_set_bit(vbb->v_buf, ptr_get_bit(vbb->vole_cache, idx_relative + vbb->row_count * column),
                column);
  }
  return vbb->v_buf;
}

const bf256_t* get_vole_aes_256(vbb_t* vbb, unsigned int idx) {
  return (bf256_t*)get_vole_row(vbb, idx);
}

const bf192_t* get_vole_aes_192(vbb_t* vbb, unsigned int idx) {
  return (bf192_t*)get_vole_row(vbb, idx);
}

const bf128_t* get_vole_aes_128(vbb_t* vbb, unsigned int idx) {
  return (bf128_t*)get_vole_row(vbb, idx);
}

const uint8_t* get_vole_u(vbb_t* vbb) {
  return vbb->vole_U;
}

const uint8_t* get_com_hash(vbb_t* vbb) {
  return vbb->com_hash;
}

// V_k cache

void add_vole_to_vk_cache(vbb_t* vbb, unsigned int idx, bf128_t* vole){
  unsigned int lambda_bytes = vbb->params->faest_param.lambda / 8;
  unsigned int offset = idx * lambda_bytes;
  memcpy(vbb->vk_cache + offset, vole, lambda_bytes);
}

void add_vole_to_vk_cache_192(vbb_t* vbb, unsigned int idx, bf192_t* vole){
  unsigned int lambda_bytes = vbb->params->faest_param.lambda / 8;
  unsigned int offset = idx * lambda_bytes;
  memcpy(vbb->vk_cache + offset, vole, lambda_bytes);
}

void add_vole_to_vk_cache_256(vbb_t* vbb, unsigned int idx, bf256_t* vole){
  unsigned int lambda_bytes = vbb->params->faest_param.lambda / 8;
  unsigned int offset = idx * lambda_bytes;
  memcpy(vbb->vk_cache + offset, vole, lambda_bytes);
}

static void setup_vk_cache(vbb_t* vbb) {
  unsigned int lambda_bytes = vbb->params->faest_param.lambda / 8;

  for (unsigned int i = 0; i < vbb->params->faest_param.lambda; i++) {
    unsigned int offset = i * lambda_bytes;
    memcpy(vbb->vk_cache + offset, get_vole_row(vbb, i), lambda_bytes);
  }
}

static inline uint8_t* get_vk(vbb_t* vbb, unsigned int idx) {
  unsigned int lambda_bytes = vbb->params->faest_param.lambda / 8;
  assert(idx < vbb->params->faest_param.Lke);

  unsigned int offset = idx * lambda_bytes;
  return (vbb->vk_cache + offset);
}

const bf128_t* get_vk_128(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_128F_LAMBDA) {
    const bf128_t* vk = (bf128_t*)get_vk(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf128_t));
    return (bf128_t*)vbb->vk_buf;
  }

  unsigned int j = idx / 32 + FAEST_128F_Nwd;
  if ((j % FAEST_128F_Nwd) == 0 || (FAEST_128F_Nwd > 6 && (j % FAEST_128F_Nwd) == 4)) {
    unsigned int i_wd       = FAEST_128F_LAMBDA;
    unsigned int factor_128 = (idx / 128) - 1;
    unsigned int offset_128 = idx % 128;
    unsigned int index      = i_wd + factor_128 * 32 + offset_128;
    const bf128_t* vk       = (bf128_t*)get_vk(vbb, index);
    memcpy(vbb->vk_buf, vk, sizeof(bf128_t));
    return (bf128_t*)vbb->vk_buf;
  }

  // Lhs recursive call
  const bf128_t* lhs_ptr = get_vk_128(vbb, idx - FAEST_128F_Nwd * 32);
  bf128_t lhs            = *lhs_ptr;
  // Rhs recursive call
  const bf128_t* rhs_ptr = get_vk_128(vbb, idx - 32);
  bf128_t rhs            = *rhs_ptr;

  bf128_t vk = bf128_add(lhs, rhs);
  memcpy(vbb->vk_buf, &vk, sizeof(bf128_t));
  return (bf128_t*)vbb->vk_buf;
}

const bf192_t* get_vk_192(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_192F_LAMBDA) {
    const bf192_t* vk = (bf192_t*)get_vk(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf192_t));
    return (bf192_t*)vbb->vk_buf;
  }

  unsigned int j = idx / 32 + FAEST_192F_Nwd;
  if ((j % FAEST_192F_Nwd) == 0 || (FAEST_192F_Nwd > 6 && (j % FAEST_192F_Nwd) == 4)) {
    unsigned int i_wd       = FAEST_192F_LAMBDA;
    unsigned int factor_192 = (idx / 192) - 1;
    unsigned int offset_192 = idx % 192;
    unsigned int index      = i_wd + factor_192 * 32 + offset_192;
    const bf192_t* vk       = (bf192_t*)get_vk(vbb, index);
    memcpy(vbb->vk_buf, vk, sizeof(bf192_t));
    return (bf192_t*)vbb->vk_buf;
  }

  // Lhs recursive call
  const bf192_t* lhs_ptr = get_vk_192(vbb, idx - FAEST_192F_Nwd * 32);
  bf192_t lhs            = *lhs_ptr;
  // Rhs recursive call
  const bf192_t* rhs_ptr = get_vk_192(vbb, idx - 32);
  bf192_t rhs            = *rhs_ptr;

  bf192_t vk = bf192_add(lhs, rhs);
  memcpy(vbb->vk_buf, &vk, 3 * sizeof(uint64_t));
  return (bf192_t*)vbb->vk_buf;
}

const bf256_t* get_vk_256(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_256F_LAMBDA) {
    const bf256_t* vk = (bf256_t*)get_vk(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf256_t));
    return (bf256_t*)vbb->vk_buf;
  }

  // We go from j=N_k to j=4(R+1)
  // In our case this is j=8 to j=4(14+1)=60
  // Based on each j we have 32 tags in vk
  unsigned int j = idx / 32 + FAEST_256F_Nwd;
  if ((j % FAEST_256F_Nwd) == 0 || (FAEST_256F_Nwd > 6 && (j % FAEST_256F_Nwd) == 4)) {
    unsigned int i_wd       = FAEST_256F_LAMBDA;
    unsigned int factor_128 = (idx / 128) - 2;
    unsigned int offset_128 = idx % 128;
    unsigned int index      = i_wd + factor_128 * 32 + offset_128;
    const bf256_t* vk       = (bf256_t*)get_vk(vbb, index);
    memcpy(vbb->vk_buf, vk, sizeof(bf256_t));
    return (bf256_t*)vbb->vk_buf;
  }

  // Lhs recursive call
  const bf256_t* lhs_ptr = get_vk_256(vbb, idx - FAEST_256F_Nwd * 32);
  bf256_t lhs            = *lhs_ptr;
  // Rhs recursive call
  const bf256_t* rhs_ptr = get_vk_256(vbb, idx - 32);
  bf256_t rhs            = *rhs_ptr;

  bf256_t vk = bf256_add(lhs, rhs);
  memcpy(vbb->vk_buf, &vk, sizeof(bf256_t));
  return (bf256_t*)vbb->vk_buf;
}
