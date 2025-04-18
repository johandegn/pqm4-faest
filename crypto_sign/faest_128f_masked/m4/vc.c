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

#include <string.h>

unsigned int NumRec(unsigned int depth, const uint8_t* bi) {
  unsigned int out = 0;
  for (unsigned int i = 0; i < depth; i++) {
    out += ((unsigned int)bi[i]) << i;
  }
  return out;
}

static void H0(const uint8_t* node, uint32_t lambda, const uint8_t* iv, uint8_t* sd, uint8_t* com) {
  const unsigned int lambda_bytes = lambda / 8;
  H0_context_t h0_ctx;
  H0_init(&h0_ctx, lambda);
  H0_update(&h0_ctx, node, lambda_bytes);
  H0_update(&h0_ctx, iv, IV_SIZE);
  H0_final(&h0_ctx, sd, lambda_bytes, com, (lambda_bytes * 2));
}

// index is the index i for (sd_i, com_i)
void get_sd_com(vec_com_t* vec_com, const uint8_t* iv, uint32_t lambda, unsigned int index,
                uint8_t* sd, uint8_t* com) {
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int depth        = vec_com->depth;

  uint8_t* children = alloca(lambda_bytes * 2);
  uint8_t* l_child  = children;
  uint8_t* r_child  = l_child + lambda_bytes;

  size_t lo         = 0;
  size_t leaf_count = (1 << depth);
  size_t hi         = leaf_count - 1;
  size_t center;

  // Find starting point from path memory
  size_t i                = 0;
  unsigned int path_index = vec_com->path.index;
  uint8_t* path_nodes     = vec_com->path.nodes;
  if (path_nodes != NULL && !vec_com->path.empty) {
    for (; i < depth;) {
      i++;
      center = (hi - lo) / 2 + lo;
      if (index <= center) { // Left
        hi = center;
        if (path_index > center)
          break;
      } else { // Right
        lo = center + 1;
        if (path_index < center + 1)
          break;
      }
    }
  }

  // Set starting node
  uint8_t* node;
  if (i > 0) {
    node = path_nodes + (i - 1) * lambda_bytes * 2;
    // if last node was right child
    if (hi != center) {
      node += lambda_bytes;
    }
  } else {
    node = vec_com->rootKey;
  }

  // Continue computing until leaf is reached
  for (; i < depth; i++) {
    prg(node, iv, children, lambda, lambda_bytes * 2);

    center = (hi - lo) / 2 + lo;
    if (index <= center) { // Left
      node = l_child;
      hi   = center;
    } else { // Right
      node = r_child;
      lo   = center + 1;
    }
    if (path_nodes != NULL) {
      memcpy(path_nodes + i * lambda_bytes * 2, node, lambda_bytes * 2);
    }
  }

  if (path_nodes != NULL) {
    vec_com->path.index = index;
    vec_com->path.empty = false;
  }

  H0(node, lambda, iv, sd, com);
  // free(children);
}

void vector_commitment(const uint8_t* rootKey, uint32_t lambda, uint32_t depth, uint8_t* path_nodes,
                       vec_com_t* vec_com) {
  const unsigned int lambda_bytes = lambda / 8;
  memcpy(vec_com->rootKey, rootKey, lambda_bytes);
  vec_com->depth      = depth;
  vec_com->path.empty = true;
  vec_com->path.nodes = path_nodes;
}

void vector_open(vec_com_t* vec_com, const uint8_t* b, uint8_t* cop, uint8_t* com_j, uint32_t depth,
                 const uint8_t* iv, uint32_t lambda) {
  // Step: 1
  const unsigned int lambda_bytes = lambda / 8;
  uint8_t* children               = alloca(lambda_bytes * 2);
  uint8_t* node                   = vec_com->rootKey;

  // Step: 3..6
  uint8_t save_left;
  for (uint32_t i = 0; i < depth; i++) {
    // b = 0 => Right
    // b = 1 => Left
    prg(node, iv, children, lambda, lambda_bytes * 2);
    save_left = b[depth - 1 - i];
    uint8_t* dst_child = children + (lambda_bytes * !save_left);
    node = children + (lambda_bytes * save_left);
    memcpy(cop + (lambda_bytes * i), dst_child, lambda_bytes);
  }

  // Step: 7
  uint64_t leaf_index = NumRec(depth, b);
  uint8_t* sd         = alloca(lambda_bytes); // Byproduct
  extract_sd_com(vec_com, iv, lambda, leaf_index, sd, com_j);
}

// Reconstruction
void vector_reconstruction(const uint8_t* cop, const uint8_t* com_j, const uint8_t* b,
                           uint32_t lambda, uint32_t depth, uint8_t* tree_nodes,
                           vec_com_rec_t* vec_com_rec) {
  const unsigned int lambda_bytes = lambda / 8;
  memcpy(vec_com_rec->nodes, cop, lambda_bytes * depth);
  memcpy(vec_com_rec->b, b, depth);
  memcpy(vec_com_rec->com_j, com_j, lambda_bytes * 2);
  vec_com_rec->depth      = depth;
  vec_com_rec->path.empty = true;
  vec_com_rec->path.nodes = tree_nodes;
}

// index is the index i for (sd_i, com_i)
void extract_sd_com(vec_com_t* vec_com, const uint8_t* iv, uint32_t lambda, unsigned int index,
                uint8_t* sd, uint8_t* com) {
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int depth        = vec_com->depth;

  uint8_t* children = alloca(lambda_bytes * 2);
  uint8_t* l_child  = children;
  uint8_t* r_child  = l_child + lambda_bytes;

  size_t lo         = 0;
  size_t leaf_count = (1 << depth);
  size_t hi         = leaf_count - 1;
  size_t center;

  // Find starting point from path memory
  size_t i                = 0;
  unsigned int path_index = vec_com->path.index;
  uint8_t* path_nodes     = vec_com->path.nodes;
  if (path_nodes != NULL && !vec_com->path.empty) {
    for (; i < depth;) {
      i++;
      center = (hi - lo) / 2 + lo;
      if (index <= center) { // Left
        hi = center;
        if (path_index > center)
          break;
      } else { // Right
        lo = center + 1;
        if (path_index < center + 1)
          break;
      }
    }
  }

  // Set starting node
  uint8_t* node;
  if (i > 0) {
    node = path_nodes + (i - 1) * lambda_bytes * 2;
    // if last node was right child
    if (hi != center) {
      node += lambda_bytes;
    }
  } else {
    node = vec_com->rootKey;
  }

  // Continue computing until leaf is reached
  for (; i < depth; i++) {
    prg(node, iv, children, lambda, lambda_bytes * 2);

    center = (hi - lo) / 2 + lo;
    if (index <= center) { // Left
      node = l_child;
      hi   = center;
    } else { // Right
      node = r_child;
      lo   = center + 1;
    }
    if (path_nodes != NULL) {
      memcpy(path_nodes + i * lambda_bytes * 2, node, lambda_bytes * 2);
    }
  }

  if (path_nodes != NULL) {
    vec_com->path.index = index;
    vec_com->path.empty = false;
  }

  H0(node, lambda, iv, sd, com);
  //free(children);
}

void extract_sd_com_rec(vec_com_rec_t* vec_com_rec, const uint8_t* iv, uint32_t lambda,
                    unsigned int index, uint8_t* sd, uint8_t* com) {
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int depth        = vec_com_rec->depth;
  uint8_t* children               = alloca(lambda_bytes * 2);
  uint8_t* l_child                = children;
  uint8_t* r_child                = l_child + lambda_bytes;

  // Assert we do not request the unknown leaf. (alternatively return zeroes)
  // assert(index != NumRec(depth, vec_com_rec->b));

  size_t lo         = 0;
  size_t leaf_count = (1 << depth);
  size_t hi         = leaf_count - 1;
  size_t center;
  uint8_t* b = vec_com_rec->b;

  // Find first known node on path
  size_t i = 0;
  for (; i < depth; i++) {
    uint8_t left_node_known = b[depth - i - 1];

    center = (hi - lo) / 2 + lo;
    if (index <= center) { // Left
      hi = center;
      if (left_node_known) {
        break;
      }
    } else { // Right
      lo = center + 1;
      if (!left_node_known) {
        break;
      }
    }
  }
  uint8_t* node = vec_com_rec->nodes + lambda_bytes * i;
  ++i;

  size_t j = 0;
  // Continue in path memory
  unsigned int path_index = vec_com_rec->path.index;
  uint8_t* path_nodes     = vec_com_rec->path.nodes;
  if (path_nodes != NULL && !vec_com_rec->path.empty && path_index <= hi && path_index >= lo) {
    for (; i < depth;) {
      i++;
      j++;
      center = (hi - lo) / 2 + lo;
      if (index <= center) { // Left
        hi = center;
        if (path_index > center)
          break;
      } else { // Right
        lo = center + 1;
        if (path_index < center + 1)
          break;
      }
    }
  }

  // Set starting node
  if (j > 0) {
    node = path_nodes + (j - 1) * lambda_bytes * 2;
    // if last node was right child
    if (hi != center) {
      node += lambda_bytes;
    }
  }

  // Continue computing until leaf is reached
  for (; i < depth; i++, j++) {
    prg(node, iv, children, lambda, lambda_bytes * 2);

    center = (hi - lo) / 2 + lo;
    if (index <= center) { // Left
      node = l_child;
      hi   = center;
    } else { // Right
      node = r_child;
      lo   = center + 1;
    }
    if (path_nodes != NULL) {
      memcpy(path_nodes + j * lambda_bytes * 2, children, lambda_bytes * 2);
    }
  }

  if (path_nodes != NULL) {
    vec_com_rec->path.index = index;
    vec_com_rec->path.empty = false;
  }

  H0(node, lambda, iv, sd, com);
}