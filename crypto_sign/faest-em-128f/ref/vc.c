/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vc.h"
#include "random_oracle.h"
#include "compat.h"
#include "aes.h"
#include "instances.h"

#include <assert.h>
#include <string.h>

unsigned int NumRec(unsigned int depth, const uint8_t* bi) {
  unsigned int out = 0;
  for (unsigned int i = 0; i < depth; i++) {
    out += ((unsigned int)bi[i]) << i;
  }
  return out;
}

// index is the index i for (sd_i, com_i)
void get_sd_com(stream_vec_com_t* sVecCom, const uint8_t* iv, uint32_t lambda, unsigned int index, uint8_t* sd, uint8_t* com) {
  const unsigned int lambdaBytes = lambda / 8;

  uint8_t children[MAX_LAMBDA_BYTES * 2];
  uint8_t* l_child = children;
  uint8_t* r_child = l_child + lambdaBytes;

  size_t lo = 0;
  size_t leaf_count = (1 << sVecCom->depth);
  size_t hi = leaf_count - 1;
  size_t center;

  // Find starting point from path memory
  size_t i = 0;
  if (sVecCom->path != NULL && sVecCom->index != sVecCom->depth) {
    for (; i < sVecCom->depth; i++) {
      center = (hi - lo) / 2 + lo;
      if (index <= center) { // Left
        if (sVecCom->index > center)
          break;
        hi = center;
      }
      else { // Right
        if (sVecCom->index < center + 1)
          break;
        lo = center + 1;
      }
    }
  }

  // Set starting node
  uint8_t* node;
  if (i > 0)
    node = sVecCom->path + (i - 1) * lambdaBytes;
  else
    node = sVecCom->rootKey;


  // Continue computing until leaf is reached
  for (; i < sVecCom->depth; i++) {
    prg(node, iv, children, lambda, lambdaBytes * 2);

    center = (hi - lo) / 2 + lo;
    if (index <= center) { // Left
      node = l_child;
      hi = center;
    }
    else { // Right
      node = r_child;
      lo = center + 1;
    }
    if (sVecCom->path != NULL)
      memcpy(sVecCom->path + i * lambdaBytes, node, lambdaBytes);
  }

  sVecCom->index = index;

  H0_context_t h0_ctx;
  H0_init(&h0_ctx, lambda);
  H0_update(&h0_ctx, node, lambdaBytes);
  H0_update(&h0_ctx, iv, IV_SIZE);
  H0_final(&h0_ctx, sd, lambdaBytes, com, (lambdaBytes * 2));
}

void get_sd_com_rec(stream_vec_com_rec_t* sVecComRec, const uint8_t* iv, uint32_t lambda, unsigned int index, uint8_t* sd, uint8_t* com) {
  const unsigned int lambdaBytes = lambda / 8;
  const unsigned int depth = sVecComRec->depth;
  uint8_t children[MAX_LAMBDA_BYTES * 2];
  uint8_t* l_child = children;
  uint8_t* r_child = l_child + lambdaBytes;

  size_t lo = 0;
  size_t leaf_count = (1 << depth);
  size_t hi = leaf_count - 1;
  size_t center;
  uint8_t* b = sVecComRec->b;

  size_t idx = NumRec(depth, b);
  size_t real_index = index ^ idx;

  // Find first known node on path
  size_t i = 0;
  for (; i < depth; i++) {
    uint8_t left_node_known = b[depth - i - 1];

    center = (hi - lo) / 2 + lo;
    if (real_index <= center) { // Left
      hi = center;
      if (left_node_known) {
        break;
      }
    }
    else { // Right
      lo = center + 1;
      if (!left_node_known) {
        break;
      }
    }
  }
  // NOTE: we could have reaced the unknown leaf if malformed signature
  uint8_t* node = sVecComRec->nodes + lambdaBytes * i;
  ++i;
  
  size_t j = 0;
  // Continue in path memory
  if (sVecComRec->path != NULL && sVecComRec->index != depth 
      && sVecComRec->index <= hi && sVecComRec->index >= lo) {
    for (; i < depth; i++, j++) {
      center = (hi - lo) / 2 + lo;
      if (real_index <= center) { // Left
        if (sVecComRec->index > center)
          break;
        hi = center;
      }
      else { // Right
        if (sVecComRec->index < center + 1)
          break;
        lo = center + 1;
      }
    }
  }
  
  // Set starting node
  if (j > 0)
    node = sVecComRec->path + (j - 1) * lambdaBytes;

  // Continue computing until leaf is reached
  for (; i < depth; i++, j++) {
    prg(node, iv, children, lambda, lambdaBytes * 2);

    center = (hi - lo) / 2 + lo;
    if (real_index <= center) { // Left
      node = l_child;
      hi = center;
    }
    else { // Right
      node = r_child;
      lo = center + 1;
    }
    if (sVecComRec->path != NULL)
      memcpy(sVecComRec->path + j * lambdaBytes, node, lambdaBytes);
  }

  sVecComRec->index = real_index;
  
  H0_context_t h0_ctx;
  H0_init(&h0_ctx, lambda);
  H0_update(&h0_ctx, node, lambdaBytes);
  H0_update(&h0_ctx, iv, IV_SIZE);
  H0_final(&h0_ctx, sd, lambdaBytes, com, (lambdaBytes * 2));
}

void stream_vector_commitment(const uint8_t* rootKey, uint32_t lambda, stream_vec_com_t* sVecCom, uint32_t depth) {
  const unsigned int lambdaBytes = lambda / 8;
  memcpy(sVecCom->rootKey, rootKey, lambdaBytes);
  sVecCom->depth = depth;
  sVecCom->path = NULL;
}

void stream_vector_open(stream_vec_com_t* sVecCom, const uint8_t* b, uint8_t* cop,
                 uint8_t* com_j, uint32_t depth,  const uint8_t* iv, uint32_t lambda) {
  // Step: 1
  const unsigned int lambdaBytes = lambda / 8;
  uint8_t children[MAX_LAMBDA_BYTES * 2];
  uint8_t* l_child = children;
  uint8_t* r_child = l_child + lambdaBytes;
  uint8_t* node = sVecCom->rootKey;

  // Step: 3..6
  uint8_t save_left;
  for (uint32_t i = 0; i < depth; i++) {
    // b = 0 => Right
    // b = 1 => Left
    prg(node, iv, children, lambda, lambdaBytes * 2);
    save_left = b[depth - 1 - i];
    if (save_left) {
      memcpy(cop + (lambdaBytes * i), l_child, lambdaBytes);
      node = r_child;
    }
    else {
      memcpy(cop + (lambdaBytes * i), r_child, lambdaBytes);
      node = l_child;
    }
  }

  // Step: 7
  uint64_t leafIndex = NumRec(depth, b);
  uint8_t sd[MAX_LAMBDA_BYTES]; // Byproduct
  get_sd_com(sVecCom, iv, lambda, leafIndex, sd, com_j);
}

void stream_vector_reconstruction(const uint8_t* cop, const uint8_t* com_j, const uint8_t* b, uint32_t lambda, uint32_t depth, stream_vec_com_rec_t* sVecComRec) {
  const unsigned int lambdaBytes = lambda / 8;
  sVecComRec->depth = depth;
  sVecComRec->index = depth; // Signals no path yet
  memcpy(sVecComRec->b, b, depth);
  memcpy(sVecComRec->nodes, cop, lambdaBytes * depth);
  memcpy(sVecComRec->com_j, com_j, lambdaBytes * 2);
}
