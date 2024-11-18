/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"

#include <stdbool.h>
#include <string.h>

int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout) {
  if (i >= t0 + t1) {
    return 0;
  }

  unsigned int lo;
  unsigned int hi;
  if (i < t0) {
    lo = i * k0;
    hi = ((i + 1) * k0);
  } else {
    unsigned int t = i - t0;
    lo             = (t0 * k0) + (t * k1);
    hi             = (t0 * k0) + ((t + 1) * k1);
  }

  assert(hi - lo == k0 || hi - lo == k1);
  for (unsigned int j = lo; j < hi; ++j) {
    // set_bit(chalout, i - lo, get_bit(chal, i));
    chalout[j - lo] = ptr_get_bit(chal, j);
  }
  return 1;
}

void partial_vole_commit_column(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                                unsigned int start, unsigned int end,
                                sign_vole_mode_ctx_t vole_mode, const faest_paramset_t* params) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth    = MAX(k0, k1);

  uint8_t* expanded_keys = alloca(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);
  uint8_t* path = alloca(lambda_bytes * max_depth * 2);
  uint8_t* r    = alloca(ellhat_bytes);
  uint8_t* sd  = alloca(lambda_bytes);
  uint8_t* com = alloca(lambda_bytes * 2);
  H1_context_t hcom_ctx;
  H1_context_t com_ctx;
  uint8_t* h = NULL;
  if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
    h = alloca(lambda_bytes * 2);
    H1_init(&hcom_ctx, lambda);
  }

  // STEP 1: To commit to [start,end] we first compute which trees we need to consider
  unsigned int depth_tau_0 = tau0 * k0;
  unsigned int k0_trees_begin = (start < depth_tau_0) ? start / k0 : tau0;
  unsigned int k1_trees_begin = (start < depth_tau_0) ? 0 : (start - depth_tau_0) / k1;
  unsigned int k0_trees_end   = (end < depth_tau_0) ? (end + (k0 - 1)) / k0 : tau0; // ceiled
  unsigned int k1_trees_end =
      (end < depth_tau_0) ? 0 : (end - depth_tau_0 + (k1 - 1)) / k1; // ceiled

  unsigned int tree_start = k0_trees_begin + k1_trees_begin;
  unsigned int tree_end   = k0_trees_end + k1_trees_end;

  // Compute the cummulative sum of the tree depths until the requested start
  unsigned int v_progress = k0 * k0_trees_begin + k1 * k1_trees_begin;

  for (unsigned int t = tree_start; t < tree_end; t++) {
    bool is_first_tree      = (t == 0);
    unsigned int tree_depth = t < tau0 ? k0 : k1;

    // v_cache_offset is used to compute the index we should write v to relative to our cache
    unsigned int v_cache_offset =
        (v_progress > start) ? v_progress - start : 0; // (i.e. MAX(v_progress-start, 0))
    // [v_start, v_end] is the v's that t provides (capped by requested start/end)
    unsigned int v_start = MAX(v_progress, start);
    unsigned int v_end   = MIN(end, v_progress + tree_depth);

    // (Setup for STEP 2)
    const unsigned int num_seeds = 1 << tree_depth;

    uint8_t* u_ptr = NULL;
    if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
      H1_init(&com_ctx, lambda);
      u_ptr = is_first_tree ? vole_mode.u : vole_mode.c + (t - 1) * ellhat_bytes;
      memset(u_ptr, 0, ellhat_bytes);
    }
    if (vole_mode.mode != EXCLUDE_V) {
      unsigned int v_count = v_end - v_start;
      memset(vole_mode.v + (v_cache_offset * ellhat_bytes), 0, v_count * ellhat_bytes);
    }

    vec_com_t vec_com;
    vector_commitment(expanded_keys + t * lambda_bytes, lambda, tree_depth, path, &vec_com);

    // STEP 2: For this tree, extract all seeds and commitments and compute according to the
    // VOLE-mode
    for (unsigned int i = 0; i < num_seeds; i++) {
      extract_sd_com(&vec_com, iv, lambda, i, sd, com);
      prg(sd, iv, r, lambda, ellhat_bytes); // Seed expansion

      if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
        int factor_32 = ellhat_bytes / 4;
        H1_update(&com_ctx, com, lambda_bytes * 2);
        xor_u32_array((uint32_t*)u_ptr, (uint32_t*)r, (uint32_t*)u_ptr, factor_32);
        xor_u8_array(u_ptr + factor_32 * 4, r + factor_32 * 4, u_ptr + factor_32 * 4,
                     ellhat_bytes - factor_32 * 4);
      }
      if (vole_mode.mode != EXCLUDE_V) {
        for (unsigned int j = v_start; j < v_end; j++) {
          // Instead of writing v_j at V[j], use the v_cache_offset
          uint8_t* write_idx = (vole_mode.v + (j - v_start + v_cache_offset) * ellhat_bytes);
          unsigned int t_v   = j - v_progress; // That is; t provides depth num of v's where t_v
                                               // reflects the current v \in [0, depth]
          // Apply r if the i/2^t_v is odd
          if ((i >> t_v) & 1) {
            int factor_32 = ellhat_bytes / 4;
            xor_u32_array((uint32_t*)write_idx, (uint32_t*)r, (uint32_t*)write_idx, factor_32);
            xor_u8_array(write_idx + factor_32 * 4, r + factor_32 * 4, write_idx + factor_32 * 4,
                         ellhat_bytes - factor_32 * 4);
          }
        }
      }
    }

    if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
      if (!is_first_tree) {
        xor_u8_array(vole_mode.u, u_ptr, u_ptr, ellhat_bytes); // Correction values
      }
      H1_final(&com_ctx, h, lambda_bytes * 2);
      H1_update(&hcom_ctx, h, lambda_bytes * 2);
    }

    v_progress += tree_depth;
  }

  if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
    H1_final(&hcom_ctx, vole_mode.hcom, lambda_bytes * 2);
    // free(h);
  }

  /*
  free(expanded_keys);
  free(path);
  */
}

void partial_vole_commit_row(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             unsigned int start, unsigned int end, const faest_paramset_t* params,
                             uint8_t* v) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth    = MAX(k0, k1);

  unsigned int len       = end - start;
  unsigned int len_bytes = (len + 7) / 8;

  uint8_t* expanded_keys = alloca(tau * lambda_bytes);
  uint8_t* sd            = alloca(lambda_bytes);
  uint8_t* com           = alloca(lambda_bytes * 2);
  uint8_t* r             = alloca(ellhat_bytes);
  uint8_t* path          = alloca(lambda_bytes * max_depth * 2);
  uint8_t* r_trunc       = alloca(len_bytes);

  vec_com_t vec_com;
  memset(v, 0, ((size_t)len) * (size_t)lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

  unsigned int col_idx = 0;
  // Iterate over each tree
  for (unsigned int t = 0; t < tau; t++) {
    unsigned int depth         = t < tau0 ? k0 : k1;
    unsigned int num_instances = 1 << depth;

    vector_commitment(expanded_keys + t * lambda_bytes, lambda, depth, path, &vec_com);

    // Iterate each seed emmited from the tree
    for (unsigned int i = 0; i < num_instances; i++) {
      extract_sd_com(&vec_com, iv, lambda, i, sd, com);
      prg(sd, iv, r, lambda, ellhat_bytes);

      // Extract and align the requested part of r
      unsigned int bit_offset = start % 8;
      unsigned int start_byte = start / 8;
      // If aligned, copy over
      if (bit_offset == 0) {
        memcpy(r_trunc, r + start_byte, len_bytes);
      } else { // If not aligned
        for (unsigned int j = 0; j < len_bytes; j++) {
          r_trunc[j] =
              (r[start_byte + j] >> bit_offset) | (r[start_byte + j + 1] << (8 - bit_offset));
        }
        // Get last part
        r_trunc[len_bytes - 1] = (r[start_byte + len_bytes - 1] >> bit_offset);
        unsigned int rest      = len - (len_bytes - 1) * 8;
        if (rest > 8 - bit_offset) {
          // Get extra part
          r_trunc[len_bytes - 1] |= r[start_byte + len_bytes] << (8 - bit_offset);
        }
      }
      // Clear final bits
      unsigned int bit_to_clear = (8 - (len % 8)) % 8;
      r_trunc[len_bytes - 1] &= (uint8_t)0xFF >> bit_to_clear;

      // XOR directly into v instead of maintaining a stack to save memory
      for (unsigned int j = 0; j < depth; j++) {
        // Only apply to correct entries
        if ((i >> j) & 1) {
          // Shift offset and XOR into v
          unsigned long row_bit_index   = (unsigned long)len * (unsigned long)(col_idx + j);
          unsigned long row_bit_offset  = row_bit_index % 8;
          unsigned long row_byte_offset = row_bit_index / 8;
          // Apply first byte
          v[row_byte_offset] ^= r_trunc[0] << row_bit_offset;
          // Apply the remaining
          for (unsigned int k = 1; k < len_bytes; k++) {
            v[row_byte_offset + k] ^=
                (r_trunc[k - 1] >> (8 - row_bit_offset)) | (r_trunc[k] << row_bit_offset);
          }

          if (row_bit_offset != 0) {
            v[row_byte_offset + len_bytes] ^= r_trunc[len_bytes - 1] >> (8 - row_bit_offset);
          }
        }
      }
    }
    col_idx += depth;
  }
  /*
  free(r);
  free(expanded_keys);
  free(path);
  free(r_trunc);
  */
}

void partial_vole_reconstruct_column(const uint8_t* iv, const uint8_t* chall,
                                     const uint8_t* const* pdec, const uint8_t* const* com_j,
                                     unsigned int ellhat, unsigned int start, unsigned int len,
                                     verify_vole_mode_ctx_t vole_mode,
                                     const faest_paramset_t* params) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int tau1         = params->faest_param.t1;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  H1_context_t hcom_ctx;
  H1_context_t com_ctx;
  uint8_t* h = NULL;
  if (vole_mode.mode != EXCLUDE_HCOM) {
    H1_init(&hcom_ctx, lambda);
    h = alloca(lambda_bytes * 2);
  }

  unsigned int max_depth = MAX(k0, k1);
  vec_com_rec_t vec_com_rec;
  vec_com_rec.b     = alloca(max_depth * sizeof(uint8_t));
  vec_com_rec.nodes = alloca(max_depth * lambda_bytes);
  memset(vec_com_rec.nodes, 0, max_depth * lambda_bytes);
  vec_com_rec.com_j   = alloca(lambda_bytes * 2);
  uint8_t* tree_nodes = alloca(lambda_bytes * (max_depth - 1) * 2);

  uint8_t* sd                  = alloca(lambda_bytes);
  uint8_t* com                 = alloca(lambda_bytes * 2);
  uint8_t* r                   = alloca(ellhat_bytes);

  unsigned int end = start + len;
  // STEP 1: To commit to [start,end] we first compute which trees we need to consider
  unsigned int depth_tau_0    = tau0 * k0;
  unsigned int k0_trees_begin = (start < depth_tau_0) ? start / k0 : tau0;
  unsigned int k1_trees_begin = (start < depth_tau_0) ? 0 : (start - depth_tau_0) / k1;
  unsigned int k0_trees_end   = (end < depth_tau_0) ? (end + (k0 - 1)) / k0 : tau0; // ceiled
  unsigned int k1_trees_end =
      (end < depth_tau_0) ? 0 : (end - depth_tau_0 + (k1 - 1)) / k1; // ceiled

  unsigned int tree_start = k0_trees_begin + k1_trees_begin;
  unsigned int tree_end   = k0_trees_end + k1_trees_end;

  // Compute the cummulative sum of the tree depths until the requested start
  unsigned int q_progress = k0 * k0_trees_begin + k1 * k1_trees_begin;

  for (unsigned int t = tree_start; t < tree_end; t++) {
    unsigned int tree_depth = t < tau0 ? k0 : k1;

    // q_cache_offset is used to compute the index we should write q to relative to our cache
    unsigned int q_cache_offset =
        (q_progress > start) ? q_progress - start : 0; // (i.e. MAX(q_progress-start, 0))
    // [q_begin, q_end] is the q's that t provides (capped by requested start/end)
    unsigned int q_begin = MAX(q_progress, start);
    unsigned int q_end   = MIN(end, q_progress + tree_depth);

    uint8_t chalout[MAX_DEPTH];
    ChalDec(chall, t, k0, tau0, k1, tau1, chalout);
    vector_reconstruction(pdec[t], com_j[t], chalout, lambda, tree_depth, tree_nodes, &vec_com_rec);
    unsigned int offset = NumRec(tree_depth, vec_com_rec.b);

    const unsigned int num_seeds = 1 << tree_depth;
    unsigned int q_count         = q_end - q_begin;

    if (vole_mode.mode != EXCLUDE_HCOM) {
      H1_init(&com_ctx, lambda);
    }
    if (vole_mode.mode != EXCLUDE_Q) {
      memset(vole_mode.q + (q_cache_offset * ellhat_bytes), 0, q_count * ellhat_bytes);
    }

    for (unsigned int i = 0; i < num_seeds; i++) {
      unsigned int offset_index = i ^ offset;
      if (offset_index == 0) {
        // As a verifier, we do not have the first seed (i.e. seed i with offset)
        if (vole_mode.mode != EXCLUDE_HCOM) {
          H1_update(&com_ctx, vec_com_rec.com_j, lambda_bytes * 2);
        }
        continue; // Skip the first seed
      }

      extract_sd_com_rec(&vec_com_rec, iv, lambda, i, sd, com);

      if (vole_mode.mode != EXCLUDE_HCOM) {
        H1_update(&com_ctx, com, lambda_bytes * 2);
      }
      if (vole_mode.mode != EXCLUDE_Q) {
        prg(sd, iv, r, lambda, ellhat_bytes);
        for (unsigned int j = q_begin; j < q_end; j++) {
          uint8_t* write_idx = (vole_mode.q + (j - q_begin + q_cache_offset) * ellhat_bytes);
          unsigned int q_v   = j - q_progress;

          // Apply r if i/2^q_v is odd
          if ((offset_index >> q_v) & 1) {
            xor_u8_array(write_idx, r, write_idx, ellhat_bytes);
          }
        }
      }
    }
    if (vole_mode.mode != EXCLUDE_HCOM) {
      H1_final(&com_ctx, h, lambda_bytes * 2);
      H1_update(&hcom_ctx, h, lambda_bytes * 2);
    }

    q_progress += tree_depth;
  }
  if (vole_mode.mode != EXCLUDE_HCOM) {
    H1_final(&hcom_ctx, vole_mode.hcom, lambda_bytes * 2);
  }
}

void partial_vole_reconstruct_row(const uint8_t* iv, const uint8_t* chall,
                                  const uint8_t* const* pdec, const uint8_t* const* com_j,
                                  uint8_t* q, unsigned int ellhat, const faest_paramset_t* params,
                                  unsigned int start, unsigned int len) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int tau1         = params->faest_param.t1;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  unsigned int len_bytes = (len + 7) / 8;

  unsigned int max_depth = MAX(k0, k1);
  vec_com_rec_t vec_com_rec;
  vec_com_rec.b     = alloca(max_depth * sizeof(uint8_t));
  vec_com_rec.nodes = alloca(max_depth * lambda_bytes);
  memset(vec_com_rec.nodes, 0, max_depth * lambda_bytes);
  vec_com_rec.com_j   = alloca(lambda_bytes * 2);
  uint8_t* tree_nodes = alloca(lambda_bytes * (max_depth - 1) * 2);

  uint8_t* sd  = alloca(lambda_bytes);
  uint8_t* com = alloca(lambda_bytes * 2);
  uint8_t* r   = alloca(ellhat_bytes);
  uint8_t* r_trunc = alloca(len_bytes);
  memset(q, 0, len * lambda_bytes);

  unsigned int col_idx = 0;
  for (unsigned int t = 0; t < tau; t++) {
    unsigned int depth = t < tau0 ? k0 : k1;
    uint8_t chalout[MAX_DEPTH];
    ChalDec(chall, t, k0, tau0, k1, tau1, chalout);
    vector_reconstruction(pdec[t], com_j[t], chalout, lambda, depth, tree_nodes, &vec_com_rec);
    unsigned int offset = NumRec(depth, vec_com_rec.b);

    const unsigned int num_instances = 1 << depth;
    for (unsigned int i = 0; i < num_instances; i++) {
      extract_sd_com_rec(&vec_com_rec, iv, lambda, i, sd, com);
      prg(sd, iv, r, lambda, ellhat_bytes);
      unsigned int offset_index = i ^ offset;

      // Extract and align the requested part of r
      unsigned int bit_offset = start % 8;
      unsigned int start_byte = start / 8;
      // If aligned, copy over
      if (bit_offset == 0) {
        memcpy(r_trunc, r + start_byte, len_bytes);
      } else { // If not aligned
        for (unsigned int j = 0; j < len_bytes; j++) {
          r_trunc[j] =
              (r[start_byte + j] >> bit_offset) | (r[start_byte + j + 1] << (8 - bit_offset));
        }
        // Get last part
        r_trunc[len_bytes - 1] = (r[start_byte + len_bytes - 1] >> bit_offset);
        unsigned int rest      = len - (len_bytes - 1) * 8;
        if (rest > 8 - bit_offset) {
          // Get extra part
          r_trunc[len_bytes - 1] |= r[start_byte + len_bytes] << (8 - bit_offset);
        }
      }
      // Clear final bits
      unsigned int bit_to_clear = (8 - (len % 8)) % 8;
      r_trunc[len_bytes - 1] &= (uint8_t)0xFF >> bit_to_clear;

      // XOR directly into v instead of maintaining a stack to save memory
      for (unsigned int j = 0; j < depth; j++) {
        // Only apply to correct entries
        if ((offset_index >> j) & 1) {
          // Shift offset and XOR into v
          unsigned long row_bit_index   = (unsigned long)len * (unsigned long)(col_idx + j);
          unsigned long row_bit_offset  = row_bit_index % 8;
          unsigned long row_byte_offset = row_bit_index / 8;
          // Apply first byte
          q[row_byte_offset] ^= r_trunc[0] << row_bit_offset;
          // Apply the remaining
          for (unsigned int k = 1; k < len_bytes; k++) {
            q[row_byte_offset + k] ^=
                (r_trunc[k - 1] >> (8 - row_bit_offset)) | (r_trunc[k] << row_bit_offset);
          }

          if (row_bit_offset != 0) {
            q[row_byte_offset + len_bytes] ^= r_trunc[len_bytes - 1] >> (8 - row_bit_offset);
          }
        }
      }
    }

    col_idx += depth;
  }
}
