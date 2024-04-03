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

void StreamConstructVole(const uint8_t* iv, stream_vec_com_t* sVecCom, unsigned int lambda, unsigned int outLenBytes, uint8_t* u, uint8_t* v, uint8_t* h) {
  unsigned int depth = sVecCom->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;

  uint8_t* stack = alloca((depth+1) * outLenBytes);
  memset(stack, 0, (depth+1) * outLenBytes);
  memset(v, 0, depth * outLenBytes);
  unsigned int stack_index = 1;
#define V(idx) (v + (idx)*outLenBytes)
#define STACK_PEAK(depth) (stack + (stack_index - 1 - depth) * outLenBytes)


  uint8_t sd[MAX_LAMBDA_BYTES];
  uint8_t com[MAX_LAMBDA_BYTES * 2];

  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  get_sd_com(sVecCom, iv, lambda, 0, sd, com);
  H1_update(&h1_ctx, com, lambda_bytes * 2);
  prg(sd, iv, stack, lambda, outLenBytes);

  unsigned int j;
  // Step: 3..4
  for (unsigned int i = 1; i < num_instances; i++) {
    get_sd_com(sVecCom, iv, lambda, i, sd, com);
    H1_update(&h1_ctx, com, lambda_bytes * 2);

    prg(sd, iv, STACK_PEAK(-1), lambda, outLenBytes);
    stack_index++;

    j = 0;
    for (unsigned int x = i + 1; !(x & 1); x >>= 1) {
      xor_u8_array(V(j), STACK_PEAK(0), V(j), outLenBytes);
      ++j;
      xor_u8_array(STACK_PEAK(0), STACK_PEAK(1), STACK_PEAK(1), outLenBytes);
      stack_index--;
    }
  }
  
  // Step: 10
  if (u != NULL) {
    memcpy(u, stack, outLenBytes);
  }

  H1_final(&h1_ctx, h, lambda_bytes * 2);
}

void StreamReconstructVole(const uint8_t* iv, stream_vec_com_rec_t* sVecComRec, unsigned int lambda, unsigned int outLenBytes, uint8_t* q, uint8_t* h) {
  unsigned int depth = sVecComRec->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;

  uint8_t* stack = alloca((depth+1) * outLenBytes);
  memset(stack, 0, (depth+1) * outLenBytes);
  memset(q, 0, depth * outLenBytes);
  unsigned int stack_index = 1;
#define Q(idx) (q + (idx)*outLenBytes)
#define STACK_PEAK(depth) (stack + (stack_index - 1 - depth) * outLenBytes)


  uint8_t sd[MAX_LAMBDA_BYTES];
  uint8_t com[MAX_LAMBDA_BYTES * 2];

  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  unsigned int j;
  // Step: 3..4
  for (unsigned int i = 1; i < num_instances; i++) {
    get_sd_com_rec(sVecComRec, iv, lambda, i, sd, com);

    prg(sd, iv, STACK_PEAK(-1), lambda, outLenBytes);
    stack_index++;

    j = 0;
    for (unsigned int x = i + 1; !(x & 1); x >>= 1) {
      xor_u8_array(Q(j), STACK_PEAK(0), Q(j), outLenBytes);
      ++j;
      xor_u8_array(STACK_PEAK(0), STACK_PEAK(1), STACK_PEAK(1), outLenBytes);
      stack_index--;
    }
  }


  // Run through to compute h
  unsigned int offset = NumRec(depth, sVecComRec->b);
  for (unsigned int i = 0; i < num_instances; i++) {
    unsigned int index = i ^ offset;
    if (index == 0) {
      H1_update(&h1_ctx, sVecComRec->com_j, lambda_bytes * 2);
      continue;
    }
    get_sd_com_rec(sVecComRec, iv, lambda, index, sd, com);
    H1_update(&h1_ctx, com, lambda_bytes * 2);
  }

  H1_final(&h1_ctx, h, lambda_bytes * 2);
}

void stream_vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                        const faest_paramset_t* params, uint8_t* hcom, stream_vec_com_t* sVecCom, uint8_t* c,
                        uint8_t* u, uint8_t** v) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  uint8_t* ui = alloca(tau * ellhat_bytes);

  // Step 1
  uint8_t* expanded_keys = alloca(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

  // for Step 12
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  uint8_t h[MAX_LAMBDA_BYTES * 2];

  uint8_t *path = alloca(MAX(k0, k1) * lambda_bytes);

  unsigned int v_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    // Step 4
    unsigned int depth = i < tau0 ? k0 : k1;
    // Step 5
    stream_vector_commitment(expanded_keys + i * lambda_bytes, lambda, &sVecCom[i], depth);
    // Step 6
    sVecCom[i].path = path;
    sVecCom[i].index = sVecCom[i].depth; // Signals no path yet
    StreamConstructVole(iv, &sVecCom[i], lambda, ellhat_bytes, ui + i * ellhat_bytes, v[v_idx], h);
    sVecCom[i].path = NULL;
    // Step 7 (and parts of 8)
    v_idx += depth;
    // Step 12 (part)
    H1_update(&h1_ctx, h, lambda_bytes * 2);
  }
  // Step 9
  memcpy(u, ui, ellhat_bytes);
  for (unsigned int i = 1; i < tau; i++) {
    // Step 11
    xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ellhat_bytes, ellhat_bytes);
  }

  // Step 12: Generating final commitment from all the com commitments
  H1_final(&h1_ctx, hcom, lambda_bytes * 2);
}

void stream_vole_reconstruct(const uint8_t* iv, const uint8_t* chall, const uint8_t* const* pdec,
                             const uint8_t* const* com_j, uint8_t* hcom, uint8_t** q, unsigned int ellhat,
                             const faest_paramset_t* params) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int tau1         = params->faest_param.t1;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  // Step 9
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  stream_vec_com_rec_t sVecComRec;
  unsigned int max_depth = MAX(k0, k1);
  sVecComRec.b = alloca(max_depth * sizeof(uint8_t));
  sVecComRec.nodes = alloca(max_depth * lambda_bytes);
  memset(sVecComRec.nodes, 0, max_depth * lambda_bytes);
  sVecComRec.com_j = alloca(lambda_bytes * 2);
  sVecComRec.path = alloca(lambda_bytes * (max_depth - 1));

  uint8_t h[MAX_LAMBDA_BYTES * 2];
  // Step: 1
  unsigned int q_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    // Step: 2
    unsigned int depth = i < tau0 ? k0 : k1;
    // Step 3
    uint8_t chalout[MAX_DEPTH];
    ChalDec(chall, i, k0, tau0, k1, tau1, chalout);

    // Step 5
    stream_vector_reconstruction(pdec[i], com_j[i], chalout, lambda, depth, &sVecComRec);

    // Step: 7..8
    StreamReconstructVole(iv, &sVecComRec, lambda, ellhat_bytes, q[q_idx], h);
    q_idx += depth;

    // Step 9
    H1_update(&h1_ctx, h, lambda_bytes * 2);
  }

  // Step: 9
  H1_final(&h1_ctx, hcom, lambda_bytes * 2);
}
