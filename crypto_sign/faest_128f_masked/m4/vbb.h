#ifndef VBB_H
#define VBB_H

#include <stdint.h>
#include <stdbool.h>

#include "vole.h"

typedef struct vbb_t vbb_t;
#include "fields.h"
#include "faest_aes.h"

typedef enum party_t {
  VERIFIER,
  SIGNER,
} party_t;

struct vbb_t {
  // Signer
  uint8_t* vole_U;
  // Verifier
  const uint8_t* sig;
  uint8_t* Dtilde_buf;
  // Common
  party_t party;
  unsigned int row_count;    // Budget of RMO computation
  unsigned int column_count; // Budget of CMO computation
  unsigned int cache_idx;
  const uint8_t* root_key;
  const faest_paramset_t* params;
  const uint8_t* iv;
  uint8_t* vole_cache;
  uint8_t* com_hash;
  // Optimizing parameters
  bool full_size;
  uint8_t* v_buf;
  // Vk_cache
  uint8_t* vk_buf;
  uint8_t* vk_cache;
  // Masking
  uint8_t* v_mask_cache;
  uint8_t* vk_mask_cache;
  uint8_t* u_mask_cache;
};

// Signer
void init_vbb_sign(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv,
                   uint8_t* c, const faest_paramset_t* params);
void init_stack_allocations_sign(vbb_t* vbb, uint8_t* hcom, uint8_t* u, uint8_t* v,
                                 uint8_t* v_buffer, uint8_t* vk_buffer, uint8_t* vk_cache);
void clean_vbb(vbb_t* vbb);
void prepare_hash_sign(vbb_t* vbb);
void prepare_aes_sign(vbb_t* vbb);
const uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx);
const bf256_t* get_vole_aes_256(vbb_t* vbb, unsigned int idx);
const bf192_t* get_vole_aes_192(vbb_t* vbb, unsigned int idx);
const bf128_t* get_vole_aes_128(vbb_t* vbb, unsigned int idx);
const uint8_t* get_vole_u(vbb_t* vbb);
const uint8_t* get_com_hash(vbb_t* vbb);
void vector_open_ondemand(vbb_t* vbb, unsigned int idx, const uint8_t* s_, uint8_t* sig_pdec,
                          uint8_t* sig_com, unsigned int depth);

// Verifier
void init_vbb_verify(vbb_t* vbb, unsigned int len, const faest_paramset_t* params,
                     const uint8_t* sig);
void init_stack_allocations_verify(vbb_t* vbb, uint8_t* hcom, uint8_t* q, uint8_t* dtilde,
                                   uint8_t* v_buffer, uint8_t* vk_buffer, uint8_t* vk_cache);
void prepare_hash_verify(vbb_t* vbb);
const uint8_t* get_vole_q_hash(vbb_t* vbb, unsigned int idx);
void prepare_aes_verify(vbb_t* vbb);
const uint8_t* get_dtilde(vbb_t* vbb, unsigned int idx);

// Vk_box
const bf128_t* get_vk_128(vbb_t* vbb, unsigned int idx);
const bf192_t* get_vk_192(vbb_t* vbb, unsigned int idx);
const bf256_t* get_vk_256(vbb_t* vbb, unsigned int idx);

// Masked
void setup_mask_storage(vbb_t* vbb, uint8_t* vk_mask, uint8_t* v_mask, uint8_t* u_mask);
void reconstruct_vole(vbb_t* vbb);
const bf128_t* get_vole_aes_128_share(vbb_t* vbb, unsigned int idx, unsigned int share);
const bf128_t* get_vk_128_share(vbb_t* vbb, unsigned int idx, unsigned int share);
const uint8_t* get_vole_u_share(vbb_t* vbb, unsigned int share);

void add_vole_to_vk_cache_share(vbb_t* vbb, unsigned int idx, bf128_t* VOLE, unsigned int share);

#endif