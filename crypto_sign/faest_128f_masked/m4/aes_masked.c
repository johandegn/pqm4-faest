
#include "aes.h"
#include "aes_masked.h"
#include "fields.h"
#include "compat.h"
#include "utils.h"
#include "randomness.h"

#define ROUNDS_128 10
#define ROUNDS_192 12
#define ROUNDS_256 14

#define KEY_WORDS_128 4
#define KEY_WORDS_192 6
#define KEY_WORDS_256 8

#define AES_BLOCK_WORDS 4
#define RIJNDAEL_BLOCK_WORDS_192 6
#define RIJNDAEL_BLOCK_WORDS_256 8

void bf8_inv_masked(bf8_t* a, bf8_t* b);

void compute_sbox_masked(bf8_t* a, bf8_t* b);

void __attribute__ ((noinline)) sub_bytes_masked(aes_block_t* state_0, aes_block_t* state_1, unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    for (unsigned int r = 0; r < AES_NR; r++) {
      bf8_t* share_0 = (*state_0)[c] + r;
      bf8_t* share_1 = (*state_1)[c] + r;
      compute_sbox_masked(share_0, share_1);
    }
  }
}

void __attribute__ ((noinline)) sub_words_masked(bf8_t* words) {
  for (int i = 0; i < 4; i++) {
    compute_sbox_masked(words + i, words + i + AES_NR);
  }
}

void __attribute__ ((noinline)) copy_first_round_key(aes_round_keys_t* round_keys_share, const uint8_t* key_share,
                          unsigned int key_words, unsigned int block_words);


void __attribute__ ((noinline)) setup_round_key_tmp(bf8_t* tmp_share, aes_round_keys_t* round_keys_share, unsigned int k, 
                         unsigned int key_words, unsigned int block_words);


void __attribute__ ((noinline)) finalize_round_key(aes_round_keys_t* round_keys_share, unsigned int k, unsigned int key_words,
                        bf8_t* tmp_share, unsigned int block_words);


void __attribute__ ((noinline)) expand_128key_masked(aes_round_keys_t* round_keys_share, const uint8_t* key_share,
                          unsigned int key_words, unsigned int block_words,
                          unsigned int num_rounds);


void aes128_init_round_keys_masked(aes_round_keys_t* round_key_share, const uint8_t* key) {
  expand_128key_masked(round_key_share, key, KEY_WORDS_128, AES_BLOCK_WORDS, ROUNDS_128);
}

uint8_t* init_round_0_key(uint8_t* w_share[2], uint8_t* w, uint8_t* w_out,
                          const faest_paramset_t* params,
                          const aes_round_keys_t round_keys_share[2]) {
  const unsigned int S_ke   = params->faest_param.Ske;
  const unsigned int lambda = params->faest_param.lambda;
  uint8_t* wtmp = w;
  // Key schedule constraints only needed for normal AES, not EM variant.
  for (unsigned int i = 0; i != params->faest_param.Nwd; ++i) {
    memcpy(w_share[0] + (wtmp - w_out), round_keys_share[0].round_keys[i / 4][i % 4],
           sizeof(aes_word_t));
    wtmp += sizeof(aes_word_t);
  }
  for (unsigned int j = 0, ik = params->faest_param.Nwd; j < S_ke / 4; ++j) {
    memcpy(w_share[0] + (wtmp - w_out), round_keys_share[0].round_keys[ik / 4][ik % 4],
           sizeof(aes_word_t));
    wtmp += sizeof(aes_word_t);
    ik += lambda == 192 ? 6 : 4;
  }


  for (unsigned int i = 0; i != params->faest_param.Nwd; ++i) {
    memcpy(w_share[1] + (w - w_out), round_keys_share[1].round_keys[i / 4][i % 4],
           sizeof(aes_word_t));
    w += sizeof(aes_word_t);
  }
  for (unsigned int j = 0, ik = params->faest_param.Nwd; j < S_ke / 4; ++j) {
    memcpy(w_share[1] + (w - w_out), round_keys_share[1].round_keys[ik / 4][ik % 4],
           sizeof(aes_word_t));
    w += sizeof(aes_word_t);
    ik += lambda == 192 ? 6 : 4;
  }
  return w;
}

void __attribute__ ((noinline)) aes_encrypt_round_masked_inner(aes_block_t* state_share, unsigned int block_words,  aes_round_keys_t* round_key, uint8_t* w_share, uint8_t** w, uint8_t* w_out, unsigned int round) {
  shift_row(state_share[0], block_words);
  store_state(w_share + (*w - w_out), state_share[0], block_words);
  mix_column(state_share[0], block_words);
  add_round_key(round, state_share[0], round_key, block_words);
}


void __attribute__ ((noinline)) aes_encrypt_round_masked(aes_block_t state_share[2], unsigned int block_words,  aes_round_keys_t round_keys_share[2], uint8_t* w_share[2], uint8_t** w, uint8_t* w_out, unsigned int round);

uint8_t* aes_extend_witness_masked_output(uint8_t* w_out, uint8_t * const w_share[2], unsigned int l);

uint8_t* aes_extend_witness_masked(const uint8_t* key_share, const uint8_t* in_share,
                                   const faest_paramset_t* params, uint8_t* w) {
  const unsigned int l          = params->faest_param.l;
  const unsigned int num_rounds = params->faest_param.R;

  // uint8_t* w           = malloc((l + 7) / 8);
  uint8_t* const w_out = w;

  uint8_t* w_share[2] = {0};
  w_share[0]          = alloca((l + 7) / 8);
  w_share[1]          = alloca((l + 7) / 8);
  memset(w_share[0], 0, (l + 7) / 8);
  memset(w_share[1], 0, (l + 7) / 8);

  unsigned int block_words = AES_BLOCK_WORDS;
  aes_round_keys_t round_keys_share[2];

  aes128_init_round_keys_masked(&round_keys_share[0], key_share);

  w = init_round_0_key(w_share, w, w_out, params, round_keys_share);

  aes_block_t state_share[2] = {0};
  load_state(state_share[0], in_share, block_words);
  load_state(state_share[1], in_share + MAX_LAMBDA_BYTES, block_words);
  add_round_key(0, state_share[0], &round_keys_share[0], block_words);
  add_round_key(0, state_share[1], &round_keys_share[1], block_words);

  for (unsigned int round = 1; round < num_rounds; ++round) {
    aes_encrypt_round_masked(state_share, block_words, round_keys_share, w_share, &w, w_out, round);
  }
  return aes_extend_witness_masked_output(w_out, w_share, l);
}
