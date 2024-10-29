
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

void bf8_inv_masked(bf8_t* in_0, bf8_t* in_1, bf8_t* out_0, bf8_t* out_1);
/*
void bf8_inv_masked(bf8_t* in_0, bf8_t* in_1, bf8_t* out_0, bf8_t* out_1) {
  bf8_t r = 0;

  //do {
  //  rand_mask(&r, 1);
  //} while(r == 0);
  
  //rand_mask(&r, 1);
  // map 0 -> 1, and otherwize no change
  //r += (!__builtin_popcount(r));
  
  bf64_t big_r = bf64_rand();
  // map 2^64-1 -> 1
  r = (bf8_t) (big_r % 255 + 1);

  bf8_t x0r = bf8_mul(*in_0, r);
  bf8_t x1r = bf8_mul(r, *in_1);
  bf8_t xr = bf8_add(x0r, x1r);
  bf8_t xr_inv = bf8_inv(xr);
  bf8_t r1 = 0;
  rand_mask(&r1, 1);
  bf8_t y0 = bf8_add(xr_inv, r1);
  bf8_t y1 = r1;
  *out_0 = bf8_mul(y0, r);
  *out_1 = bf8_mul(r, y1);
}
*/

void compute_sbox_masked(bf8_t* in_0, bf8_t* in_1, bf8_t* out_0, bf8_t* out_1);

/*
void __attribute__ ((noinline)) compute_sbox_masked(bf8_t* in_0, bf8_t* in_1, bf8_t* out_0, bf8_t* out_1) {
  bf8_inv_masked(in_0, in_1, out_0, out_1);

  *out_0 = affine_incomplete(*out_0);
  *out_1 = affine(*out_1);
}
*/

void __attribute__ ((noinline)) sub_bytes_masked(aes_block_t* state_0, aes_block_t* state_1, unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    for (unsigned int r = 0; r < AES_NR; r++) {
      bf8_t* share_0 = (*state_0)[c] + r;
      bf8_t* share_1 = (*state_1)[c] + r;
      compute_sbox_masked(share_0, share_1, share_0, share_1);
    }
  }
}

void __attribute__ ((noinline)) sub_words_masked(bf8_t* words) {
  for (int i = 0; i < 4; i++) {
    compute_sbox_masked(words + i, words + i + AES_NR, words + i, words + i + AES_NR);
  }
}

void expand_128key_masked(aes_round_keys_t* round_keys_share, const uint8_t* key_share,
                          unsigned int key_words, unsigned int block_words,
                          unsigned int num_rounds);
/*
void expand_128key_masked(aes_round_keys_t* round_keys_share, const uint8_t* key_share,
                          unsigned int key_words, unsigned int block_words,
                          unsigned int num_rounds) {
  for (unsigned int k = 0; k < key_words; k++) {
    round_keys_share[0].round_keys[k / block_words][k % block_words][0] =
        bf8_load(&key_share[4 * k]);
    round_keys_share[0].round_keys[k / block_words][k % block_words][1] =
        bf8_load(&key_share[(4 * k) + 1]);
    round_keys_share[0].round_keys[k / block_words][k % block_words][2] =
        bf8_load(&key_share[(4 * k) + 2]);
    round_keys_share[0].round_keys[k / block_words][k % block_words][3] =
        bf8_load(&key_share[(4 * k) + 3]);
  }
  for (unsigned int k = 0; k < key_words; k++) {
    round_keys_share[1].round_keys[k / block_words][k % block_words][0] =
        bf8_load(&key_share[4 * k + MAX_LAMBDA_BYTES]);
    round_keys_share[1].round_keys[k / block_words][k % block_words][1] =
        bf8_load(&key_share[(4 * k) + 1 + MAX_LAMBDA_BYTES]);
    round_keys_share[1].round_keys[k / block_words][k % block_words][2] =
        bf8_load(&key_share[(4 * k) + 2 + MAX_LAMBDA_BYTES]);
    round_keys_share[1].round_keys[k / block_words][k % block_words][3] =
        bf8_load(&key_share[(4 * k) + 3 + MAX_LAMBDA_BYTES]);
  }

  for (unsigned int k = key_words; k < block_words * (num_rounds + 1); ++k) {
    bf8_t tmp_share[2][AES_NR];
    for (int i = 0; i < AES_NR; i++){
      tmp_share[0][i] = round_keys_share[0].round_keys[(k - 1) / block_words][(k - 1) % block_words][i];
    }
    for (int i = 0; i < AES_NR; i++){
      tmp_share[1][i] = round_keys_share[1].round_keys[(k - 1) / block_words][(k - 1) % block_words][i];
    }
    if (k % key_words == 0) {
      rot_word(tmp_share[0]);
      rot_word(tmp_share[1]);
      sub_words_masked(&tmp_share[0][0]);
      tmp_share[0][0] ^= round_constants((k / key_words) - 1);
    }

    unsigned int m = k - key_words;
    round_keys_share[0].round_keys[k / block_words][k % block_words][0] =
        round_keys_share[0].round_keys[m / block_words][m % block_words][0] ^ tmp_share[0][0];
    round_keys_share[0].round_keys[k / block_words][k % block_words][1] =
        round_keys_share[0].round_keys[m / block_words][m % block_words][1] ^ tmp_share[0][1];
    round_keys_share[0].round_keys[k / block_words][k % block_words][2] =
        round_keys_share[0].round_keys[m / block_words][m % block_words][2] ^ tmp_share[0][2];
    round_keys_share[0].round_keys[k / block_words][k % block_words][3] =
        round_keys_share[0].round_keys[m / block_words][m % block_words][3] ^ tmp_share[0][3];
    round_keys_share[1].round_keys[k / block_words][k % block_words][0] =
        round_keys_share[1].round_keys[m / block_words][m % block_words][0] ^ tmp_share[1][0];
    round_keys_share[1].round_keys[k / block_words][k % block_words][1] =
        round_keys_share[1].round_keys[m / block_words][m % block_words][1] ^ tmp_share[1][1];
    round_keys_share[1].round_keys[k / block_words][k % block_words][2] =
        round_keys_share[1].round_keys[m / block_words][m % block_words][2] ^ tmp_share[1][2];
    round_keys_share[1].round_keys[k / block_words][k % block_words][3] =
        round_keys_share[1].round_keys[m / block_words][m % block_words][3] ^ tmp_share[1][3];
  }
}
*/

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
/*
void __attribute__ ((noinline)) aes_encrypt_round_masked(aes_block_t state_share[2], unsigned int block_words,  aes_round_keys_t round_keys_share[2], uint8_t* w_share[2], uint8_t** w, uint8_t* w_out, unsigned int round) {
  sub_bytes_masked(&state_share[0], &state_share[1], block_words);

  aes_encrypt_round_masked_inner(&state_share[0], block_words, &round_keys_share[0], w_share[0], w, w_out, round);
  aes_encrypt_round_masked_inner(&state_share[1], block_words, &round_keys_share[1], w_share[1], w, w_out, round);

  *w += sizeof(aes_word_t) * block_words;
}
*/

uint8_t* aes_extend_witness_masked_output(uint8_t* w_out, uint8_t * const w_share[2], unsigned int l);
/*
uint8_t* aes_extend_witness_masked_output(uint8_t* w_out, uint8_t * const w_share[2], unsigned int l) {
  // Setting up output to return
  for (unsigned int i = 0; i < (l + 7) / 8; i++) {
    w_out[i]               = w_share[0][i];
    w_out[i + (l + 7) / 8] = w_share[1][i];
  }
  return w_out;
}
*/

uint8_t* aes_extend_witness_masked(const uint8_t* key_share, const uint8_t* in_share,
                                   const faest_paramset_t* params, uint8_t* w) {
  const unsigned int lambda     = params->faest_param.lambda;
  const unsigned int l          = params->faest_param.l;
  const unsigned int L_ke       = params->faest_param.Lke;
  const unsigned int num_rounds = params->faest_param.R;

  // uint8_t* w           = malloc((l + 7) / 8);
  uint8_t* const w_out = w;

  uint8_t* w_share[2] = {0};
  w_share[0]          = alloca((l + 7) / 8);
  w_share[1]          = alloca((l + 7) / 8);
  memset(w_share[0], 0, (l + 7) / 8);
  memset(w_share[1], 0, (l + 7) / 8);

  unsigned int block_words = AES_BLOCK_WORDS;
  unsigned int beta        = 1;
  switch (params->faest_paramid) {
  case FAEST_192F:
  case FAEST_192S:
  case FAEST_256F:
  case FAEST_256S:
    beta = 2;
    break;
  case FAEST_EM_192F:
  case FAEST_EM_192S:
    block_words = RIJNDAEL_BLOCK_WORDS_192;
    break;
  case FAEST_EM_256F:
  case FAEST_EM_256S:
    block_words = RIJNDAEL_BLOCK_WORDS_256;
    break;
  default:
    break;
  }

  // NOTE - Reconstruct key if not running 128 variant
  // NOTE - All other variants are not masked yet!
  uint8_t* key = NULL;
  uint8_t* in = NULL;
  if (!L_ke || lambda != 128) {
    key = alloca(MAX_LAMBDA_BYTES);
    for (int i = 0; i < MAX_LAMBDA_BYTES; i++) {
      key[i] = key_share[i] ^ key_share[i + MAX_LAMBDA_BYTES];
    }
    in = alloca(MAX_LAMBDA_BYTES);
    for (int i = 0; i < MAX_LAMBDA_BYTES; i++) {
      in[i] = in_share[i] ^ in_share[i + MAX_LAMBDA_BYTES];
    }
  }
  if (!L_ke) {
    // switch input and key for EM
    uint8_t* tmp = key;
    key          = in;
    in           = tmp;
  }
  aes_round_keys_t round_keys_share[2];

  // Step 3
  aes_round_keys_t round_keys;
  switch (lambda) {
  case 256:
    if (block_words == RIJNDAEL_BLOCK_WORDS_256) {
      rijndael256_init_round_keys(&round_keys, key);
    } else {
      aes256_init_round_keys(&round_keys, key);
    }
    for (int i = 0; i < AES_MAX_ROUNDS + 1; i++) {
      for (int j = 0; j < 8; j++) {
        for (int k = 0; k < 4; k++) {
          rand_mask(&round_keys_share[0].round_keys[i][j][k], 1);
          round_keys_share[1].round_keys[i][j][k] =
              round_keys.round_keys[i][j][k] ^ round_keys_share[0].round_keys[i][j][k];
        }
      }
    }
    break;
  case 192:
    if (block_words == RIJNDAEL_BLOCK_WORDS_192) {
      rijndael192_init_round_keys(&round_keys, key);
    } else {
      aes192_init_round_keys(&round_keys, key);
    }
    for (int i = 0; i < AES_MAX_ROUNDS + 1; i++) {
      for (int j = 0; j < 8; j++) {
        for (int k = 0; k < 4; k++) {
          rand_mask(&round_keys_share[0].round_keys[i][j][k], 1);
          round_keys_share[1].round_keys[i][j][k] =
              round_keys.round_keys[i][j][k] ^ round_keys_share[0].round_keys[i][j][k];
        }
      }
    }
    break;
  default:
    if (!L_ke) {
      aes128_init_round_keys(&round_keys, key);
      for (int i = 0; i < AES_MAX_ROUNDS + 1; i++) {
        for (int j = 0; j < 8; j++) {
          for (int k = 0; k < 4; k++) {
            rand_mask(&round_keys_share[0].round_keys[i][j][k], 1);
            round_keys_share[1].round_keys[i][j][k] =
                round_keys.round_keys[i][j][k] ^ round_keys_share[0].round_keys[i][j][k];
          }
        }
      }
      break;
    }
    aes128_init_round_keys_masked(&round_keys_share[0], key_share);
    break;
  }

  // Saving the expanded key parts to the extended witness
  // Step 4
  if (L_ke > 0) {
    w = init_round_0_key(w_share, w, w_out, params, round_keys_share);
  } else {
    // saving the OWF key to the extended witness
    memcpy(w_share[0] + (w - w_out), in, lambda / 8);
    memset(w_share[1] + (w - w_out), 0, lambda / 8);
    w += lambda / 8;
  }

  // Step 10
  for (unsigned b = 0; b < beta; ++b, in += sizeof(aes_word_t) * block_words) {
    // Step 12

    // Masking the input state for encryption
    // The state we mask are taken from the pk, hence it is public, the all states after are secret.
    aes_block_t state_share[2] = {0};
    if (!L_ke || lambda != 128) {
      aes_block_t state;
      load_state(state, in, block_words);
      for (unsigned int c = 0; c < block_words; c++) {
        for (unsigned int r = 0; r < AES_NR; r++) {
          rand_mask(&state_share[0][c][r], 1);
          state_share[1][c][r] = state[c][r] ^ state_share[0][c][r];
        }
      }
    } else {
      load_state(state_share[0], in_share, block_words);
      load_state(state_share[1], in_share + MAX_LAMBDA_BYTES, block_words);
    }
    // Step 13
    add_round_key(0, state_share[0], &round_keys_share[0], block_words);
    add_round_key(0, state_share[1], &round_keys_share[1], block_words);

    for (unsigned int round = 1; round < num_rounds; ++round) {
      aes_encrypt_round_masked(state_share, block_words, round_keys_share, w_share, &w, w_out, round);
    }
    // last round is not commited to, so not computed
  }

  // Setting up output to return
  return aes_extend_witness_masked_output(w_out, w_share, l);
}
