#include "config.h"
/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "aes.h"

#include "fields.h"
#include "compat.h"
#include "utils.h"
#include "randomness.h"

#if defined(PQCLEAN)
#include "aes-publicinputs.h"
#endif

#if defined(HAVE_OPENSSL)
#include <openssl/evp.h>
#endif
#include <string.h>

#define ROUNDS_128 10
#define ROUNDS_192 12
#define ROUNDS_256 14

#define KEY_WORDS_128 4
#define KEY_WORDS_192 6
#define KEY_WORDS_256 8

#define AES_BLOCK_WORDS 4
#define RIJNDAEL_BLOCK_WORDS_192 6
#define RIJNDAEL_BLOCK_WORDS_256 8

static const bf8_t round_constants_table[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

bf8_t __attribute__ ((noinline)) round_constants(bf8_t b) {
  return round_constants_table[b];
}

static int contains_zero(const bf8_t* block) {
  return !block[0] | !block[1] | !block[2] | !block[3];
}
static const bf8_t affine_table[256] = {
    0x63, 0x7c, 0x5d, 0x42, 0x1f, 0x00, 0x21, 0x3e, 0x9b, 0x84, 0xa5, 0xba, 0xe7, 0xf8, 0xd9, 0xc6,
    0x92, 0x8d, 0xac, 0xb3, 0xee, 0xf1, 0xd0, 0xcf, 0x6a, 0x75, 0x54, 0x4b, 0x16, 0x09, 0x28, 0x37,
    0x80, 0x9f, 0xbe, 0xa1, 0xfc, 0xe3, 0xc2, 0xdd, 0x78, 0x67, 0x46, 0x59, 0x04, 0x1b, 0x3a, 0x25,
    0x71, 0x6e, 0x4f, 0x50, 0x0d, 0x12, 0x33, 0x2c, 0x89, 0x96, 0xb7, 0xa8, 0xf5, 0xea, 0xcb, 0xd4,
    0xa4, 0xbb, 0x9a, 0x85, 0xd8, 0xc7, 0xe6, 0xf9, 0x5c, 0x43, 0x62, 0x7d, 0x20, 0x3f, 0x1e, 0x01,
    0x55, 0x4a, 0x6b, 0x74, 0x29, 0x36, 0x17, 0x08, 0xad, 0xb2, 0x93, 0x8c, 0xd1, 0xce, 0xef, 0xf0,
    0x47, 0x58, 0x79, 0x66, 0x3b, 0x24, 0x05, 0x1a, 0xbf, 0xa0, 0x81, 0x9e, 0xc3, 0xdc, 0xfd, 0xe2,
    0xb6, 0xa9, 0x88, 0x97, 0xca, 0xd5, 0xf4, 0xeb, 0x4e, 0x51, 0x70, 0x6f, 0x32, 0x2d, 0x0c, 0x13,
    0xec, 0xf3, 0xd2, 0xcd, 0x90, 0x8f, 0xae, 0xb1, 0x14, 0x0b, 0x2a, 0x35, 0x68, 0x77, 0x56, 0x49,
    0x1d, 0x02, 0x23, 0x3c, 0x61, 0x7e, 0x5f, 0x40, 0xe5, 0xfa, 0xdb, 0xc4, 0x99, 0x86, 0xa7, 0xb8,
    0x0f, 0x10, 0x31, 0x2e, 0x73, 0x6c, 0x4d, 0x52, 0xf7, 0xe8, 0xc9, 0xd6, 0x8b, 0x94, 0xb5, 0xaa,
    0xfe, 0xe1, 0xc0, 0xdf, 0x82, 0x9d, 0xbc, 0xa3, 0x06, 0x19, 0x38, 0x27, 0x7a, 0x65, 0x44, 0x5b,
    0x2b, 0x34, 0x15, 0x0a, 0x57, 0x48, 0x69, 0x76, 0xd3, 0xcc, 0xed, 0xf2, 0xaf, 0xb0, 0x91, 0x8e,
    0xda, 0xc5, 0xe4, 0xfb, 0xa6, 0xb9, 0x98, 0x87, 0x22, 0x3d, 0x1c, 0x03, 0x5e, 0x41, 0x60, 0x7f,
    0xc8, 0xd7, 0xf6, 0xe9, 0xb4, 0xab, 0x8a, 0x95, 0x30, 0x2f, 0x0e, 0x11, 0x4c, 0x53, 0x72, 0x6d,
    0x39, 0x26, 0x07, 0x18, 0x45, 0x5a, 0x7b, 0x64, 0xc1, 0xde, 0xff, 0xe0, 0xbd, 0xa2, 0x83, 0x9c};

static const bf8_t affine_incomplete_table[256] = {
    0x00, 0x1f, 0x3e, 0x21, 0x7c, 0x63, 0x42, 0x5d, 0xf8, 0xe7, 0xc6, 0xd9, 0x84, 0x9b, 0xba, 0xa5,
    0xf1, 0xee, 0xcf, 0xd0, 0x8d, 0x92, 0xb3, 0xac, 0x09, 0x16, 0x37, 0x28, 0x75, 0x6a, 0x4b, 0x54,
    0xe3, 0xfc, 0xdd, 0xc2, 0x9f, 0x80, 0xa1, 0xbe, 0x1b, 0x04, 0x25, 0x3a, 0x67, 0x78, 0x59, 0x46,
    0x12, 0x0d, 0x2c, 0x33, 0x6e, 0x71, 0x50, 0x4f, 0xea, 0xf5, 0xd4, 0xcb, 0x96, 0x89, 0xa8, 0xb7,
    0xc7, 0xd8, 0xf9, 0xe6, 0xbb, 0xa4, 0x85, 0x9a, 0x3f, 0x20, 0x01, 0x1e, 0x43, 0x5c, 0x7d, 0x62,
    0x36, 0x29, 0x08, 0x17, 0x4a, 0x55, 0x74, 0x6b, 0xce, 0xd1, 0xf0, 0xef, 0xb2, 0xad, 0x8c, 0x93,
    0x24, 0x3b, 0x1a, 0x05, 0x58, 0x47, 0x66, 0x79, 0xdc, 0xc3, 0xe2, 0xfd, 0xa0, 0xbf, 0x9e, 0x81,
    0xd5, 0xca, 0xeb, 0xf4, 0xa9, 0xb6, 0x97, 0x88, 0x2d, 0x32, 0x13, 0x0c, 0x51, 0x4e, 0x6f, 0x70,
    0x8f, 0x90, 0xb1, 0xae, 0xf3, 0xec, 0xcd, 0xd2, 0x77, 0x68, 0x49, 0x56, 0x0b, 0x14, 0x35, 0x2a,
    0x7e, 0x61, 0x40, 0x5f, 0x02, 0x1d, 0x3c, 0x23, 0x86, 0x99, 0xb8, 0xa7, 0xfa, 0xe5, 0xc4, 0xdb,
    0x6c, 0x73, 0x52, 0x4d, 0x10, 0x0f, 0x2e, 0x31, 0x94, 0x8b, 0xaa, 0xb5, 0xe8, 0xf7, 0xd6, 0xc9,
    0x9d, 0x82, 0xa3, 0xbc, 0xe1, 0xfe, 0xdf, 0xc0, 0x65, 0x7a, 0x5b, 0x44, 0x19, 0x06, 0x27, 0x38,
    0x48, 0x57, 0x76, 0x69, 0x34, 0x2b, 0x0a, 0x15, 0xb0, 0xaf, 0x8e, 0x91, 0xcc, 0xd3, 0xf2, 0xed,
    0xb9, 0xa6, 0x87, 0x98, 0xc5, 0xda, 0xfb, 0xe4, 0x41, 0x5e, 0x7f, 0x60, 0x3d, 0x22, 0x03, 0x1c,
    0xab, 0xb4, 0x95, 0x8a, 0xd7, 0xc8, 0xe9, 0xf6, 0x53, 0x4c, 0x6d, 0x72, 0x2f, 0x30, 0x11, 0x0e,
    0x5a, 0x45, 0x64, 0x7b, 0x26, 0x39, 0x18, 0x07, 0xa2, 0xbd, 0x9c, 0x83, 0xde, 0xc1, 0xe0, 0xff};

bf8_t affine(bf8_t b) {
  return affine_table[b];
}

bf8_t affine_incomplete(bf8_t b) {
  return affine_incomplete_table[b];
}

static bf8_t compute_sbox(bf8_t in) {
  bf8_t t = bf8_inv(in);
  return affine_table[t];
}

// Moved to assmelby to ensure constant time
void aes_increment_iv(uint8_t* iv);
/*
void aes_increment_iv(uint8_t* iv) {
  uint8_t val = 1;
  // Constant time increment
  for (unsigned int i = 16; i > 0; i--) {
    iv[i - 1] += val;
    val = val * !iv[i - 1];
  }
}
*/

// ## AES ##
// Round Functions
static void __attribute__ ((noinline)) add_round_key(unsigned int round, aes_block_t state, const aes_round_keys_t* round_key,
                          unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    xor_u8_array(&state[c][0], &round_key->round_keys[round][c][0], &state[c][0], AES_NR);
  }
}

static int sub_bytes(aes_block_t state, unsigned int block_words) {
  int ret = 0;

  for (unsigned int c = 0; c < block_words; c++) {
    ret |= contains_zero(&state[c][0]);
    for (unsigned int r = 0; r < AES_NR; r++) {
      state[c][r] = compute_sbox(state[c][r]);
    }
  }

  return ret;
}

void shift_row(aes_block_t state, unsigned int block_words);
/*
void shift_row(aes_block_t state, unsigned int block_words) {
  aes_block_t new_state = {0};
  switch (block_words) {
  case 4:
  case 6:
    for (unsigned int i = 0; i < block_words; ++i) {
      new_state[i][0] = state[i][0];
      new_state[i][1] = state[(i + 1) % block_words][1];
      new_state[i][2] = state[(i + 2) % block_words][2];
      new_state[i][3] = state[(i + 3) % block_words][3];
    }
    break;
  case 8:
    for (unsigned int i = 0; i < block_words; i++) {
      new_state[i][0] = state[i][0];
      new_state[i][1] = state[(i + 1) % 8][1];
      new_state[i][2] = state[(i + 3) % 8][2];
      new_state[i][3] = state[(i + 4) % 8][3];
    }
    break;
  }

  for (unsigned int i = 0; i < block_words; ++i) {
    memcpy(&state[i][0], &new_state[i][0], AES_NR);
  }
}
*/

void mix_column(aes_block_t state, unsigned int block_words);
/*
void mix_column(aes_block_t state, unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    bf8_t tmp = bf8_mul(state[c][0], 0x02) ^ bf8_mul(state[c][1], 0x03) ^ state[c][2] ^ state[c][3];
    bf8_t tmp_1 =
        state[c][0] ^ bf8_mul(state[c][1], 0x02) ^ bf8_mul(state[c][2], 0x03) ^ state[c][3];
    bf8_t tmp_2 =
        state[c][0] ^ state[c][1] ^ bf8_mul(state[c][2], 0x02) ^ bf8_mul(state[c][3], 0x03);
    bf8_t tmp_3 =
        bf8_mul(state[c][0], 0x03) ^ state[c][1] ^ state[c][2] ^ bf8_mul(state[c][3], 0x02);

    state[c][0] = tmp;
    state[c][1] = tmp_1;
    state[c][2] = tmp_2;
    state[c][3] = tmp_3;
  }
}
*/

// Key Expansion functions
void sub_words(bf8_t* words) {
  words[0] = compute_sbox(words[0]);
  words[1] = compute_sbox(words[1]);
  words[2] = compute_sbox(words[2]);
  words[3] = compute_sbox(words[3]);
}

void rot_word(bf8_t* words);
/*
void __attribute__ ((noinline)) rot_word(bf8_t* words) {
  bf8_t tmp = words[0];
  words[0]  = words[1];
  words[1]  = words[2];
  words[2]  = words[3];
  words[3]  = tmp;
}
*/

int expand_key(aes_round_keys_t* round_keys, const uint8_t* key, unsigned int key_words,
               unsigned int block_words, unsigned int num_rounds) {
  int ret = 0;

  for (unsigned int k = 0; k < key_words; k++) {
    round_keys->round_keys[k / block_words][k % block_words][0] = bf8_load(&key[4 * k]);
    round_keys->round_keys[k / block_words][k % block_words][1] = bf8_load(&key[(4 * k) + 1]);
    round_keys->round_keys[k / block_words][k % block_words][2] = bf8_load(&key[(4 * k) + 2]);
    round_keys->round_keys[k / block_words][k % block_words][3] = bf8_load(&key[(4 * k) + 3]);
  }

  for (unsigned int k = key_words; k < block_words * (num_rounds + 1); ++k) {
    bf8_t tmp[AES_NR];
    memcpy(tmp, round_keys->round_keys[(k - 1) / block_words][(k - 1) % block_words], sizeof(tmp));

    if (k % key_words == 0) {
      rot_word(tmp);
      ret |= contains_zero(tmp);
      sub_words(tmp);
      tmp[0] ^= round_constants((k / key_words) - 1);
    }

    if (key_words > 6 && (k % key_words) == 4) {
      ret |= contains_zero(tmp);
      sub_words(tmp);
    }

    unsigned int m = k - key_words;
    round_keys->round_keys[k / block_words][k % block_words][0] =
        round_keys->round_keys[m / block_words][m % block_words][0] ^ tmp[0];
    round_keys->round_keys[k / block_words][k % block_words][1] =
        round_keys->round_keys[m / block_words][m % block_words][1] ^ tmp[1];
    round_keys->round_keys[k / block_words][k % block_words][2] =
        round_keys->round_keys[m / block_words][m % block_words][2] ^ tmp[2];
    round_keys->round_keys[k / block_words][k % block_words][3] =
        round_keys->round_keys[m / block_words][m % block_words][3] ^ tmp[3];
  }

  return ret;
}

// Calling Functions

int aes128_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_128, AES_BLOCK_WORDS, ROUNDS_128);
}

int aes192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_192, AES_BLOCK_WORDS, ROUNDS_192);
}

int aes256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_256, AES_BLOCK_WORDS, ROUNDS_256);
}

int rijndael192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_192, RIJNDAEL_BLOCK_WORDS_192, ROUNDS_192);
}

int rijndael256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_256, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
}

static void load_state(aes_block_t state, const uint8_t* src, unsigned int block_words) {
  for (unsigned int i = 0; i != block_words * 4; ++i) {
    state[i / 4][i % 4] = bf8_load(&src[i]);
  }
}

static void __attribute__ ((noinline)) store_state(uint8_t* dst, aes_block_t state, unsigned int block_words) {
  for (unsigned int i = 0; i != block_words * 4; ++i) {
    bf8_store(&dst[i], state[i / 4][i % 4]);
  }
}

static int aes_encrypt(const aes_round_keys_t* keys, aes_block_t state, unsigned int block_words,
                       unsigned int num_rounds) {
  int ret = 0;

  // first round
  add_round_key(0, state, keys, block_words);

  for (unsigned int round = 1; round < num_rounds; ++round) {
    ret |= sub_bytes(state, block_words);
    shift_row(state, block_words);
    mix_column(state, block_words);
    add_round_key(round, state, keys, block_words);
  }

  // last round
  ret |= sub_bytes(state, block_words);
  shift_row(state, block_words);
  add_round_key(num_rounds, state, keys, block_words);

  return ret;
}

int aes128_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  const int ret = aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_128);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
  return ret;
}

int aes192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  const int ret = aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_192);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
  return ret;
}

int aes256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  const int ret = aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_256);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
  return ret;
}

int rijndael192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                              uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, RIJNDAEL_BLOCK_WORDS_192);
  const int ret = aes_encrypt(key, state, RIJNDAEL_BLOCK_WORDS_192, ROUNDS_192);
  store_state(ciphertext, state, RIJNDAEL_BLOCK_WORDS_192);
  return ret;
}

int rijndael256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                              uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, RIJNDAEL_BLOCK_WORDS_256);
  const int ret = aes_encrypt(key, state, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
  store_state(ciphertext, state, RIJNDAEL_BLOCK_WORDS_256);
  return ret;
}

void prg(const uint8_t* key, const uint8_t* iv, uint8_t* out, unsigned int seclvl, size_t outlen) {
#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher;
  switch (seclvl) {
  case 256:
    cipher = EVP_aes_256_ctr();
    break;
  case 192:
    cipher = EVP_aes_192_ctr();
    break;
  default:
    cipher = EVP_aes_128_ctr();
    break;
  }

  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, iv);

  static const uint8_t plaintext[16] = {0};

  int len = 0;
  for (size_t idx = 0; idx < outlen / 16; idx += 1, out += 16) {
    EVP_EncryptUpdate(ctx, out, &len, plaintext, sizeof(plaintext));
  }
  if (outlen % 16) {
    EVP_EncryptUpdate(ctx, out, &len, plaintext, outlen % 16);
  }
  EVP_EncryptFinal_ex(ctx, out, &len);
  EVP_CIPHER_CTX_free(ctx);
#elif defined(PQCLEAN)
  uint8_t internal_iv[16];
  memcpy(internal_iv, iv, sizeof(internal_iv));

  switch (seclvl) {
  case 256:
    aes256ctx_publicinputs ctx256;
    aes256_ecb_keyexp_publicinputs(&ctx256, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes256_ecb_publicinputs(out, internal_iv, 1, &ctx256);
      aes_increment_iv(internal_iv);
    }
    if (outlen % 16) {
      uint8_t tmp[16];
      aes256_ecb_publicinputs(tmp, internal_iv, 1, &ctx256);
      memcpy(out, tmp, outlen);
    }
    aes256_ctx_release_publicinputs(&ctx256);
    break;
  case 192:
    aes192ctx_publicinputs ctx192;
    aes192_ecb_keyexp_publicinputs(&ctx192, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes192_ecb_publicinputs(out, internal_iv, 1, &ctx192);
      aes_increment_iv(internal_iv);
    }
    if (outlen % 16) {
      uint8_t tmp[16];
      aes192_ecb_publicinputs(tmp, internal_iv, 1, &ctx192);
      memcpy(out, tmp, outlen);
    }
    aes192_ctx_release_publicinputs(&ctx192);
    break;
  default:
    aes128ctx_publicinputs ctx128;
    aes128_ecb_keyexp_publicinputs(&ctx128, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes128_ecb_publicinputs(out, internal_iv, 1, &ctx128);
      aes_increment_iv(internal_iv);
    }
    if (outlen % 16) {
      uint8_t tmp[16];
      aes128_ecb_publicinputs(tmp, internal_iv, 1, &ctx128);
      memcpy(out, tmp, outlen);
    }
    aes128_ctx_release_publicinputs(&ctx128);
    break;
  }

#else
  uint8_t internal_iv[16];
  memcpy(internal_iv, iv, sizeof(internal_iv));

  aes_round_keys_t round_key;

  switch (seclvl) {
  case 256:
    aes256_init_round_keys(&round_key, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_256);
      store_state(out, state, AES_BLOCK_WORDS);
      aes_increment_iv(internal_iv);
    }
    if (outlen) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_256);
      uint8_t tmp[16];
      store_state(tmp, state, AES_BLOCK_WORDS);
      memcpy(out, tmp, outlen);
    }
    return;
  case 192:
    aes192_init_round_keys(&round_key, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_192);
      store_state(out, state, AES_BLOCK_WORDS);
      aes_increment_iv(internal_iv);
    }
    if (outlen) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_192);
      uint8_t tmp[16];
      store_state(tmp, state, AES_BLOCK_WORDS);
      memcpy(out, tmp, outlen);
    }
    return;
  default:
    aes128_init_round_keys(&round_key, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_128);
      store_state(out, state, AES_BLOCK_WORDS);
      aes_increment_iv(internal_iv);
    }
    if (outlen) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_128);
      uint8_t tmp[16];
      store_state(tmp, state, AES_BLOCK_WORDS);
      memcpy(out, tmp, outlen);
    }
    return;
  }
#endif
}

uint8_t* aes_extend_witness(const uint8_t* key, const uint8_t* in, const faest_paramset_t* params,
                            uint8_t* w) {
  const unsigned int lambda = params->faest_param.lambda;
  // const unsigned int l          = params->faest_param.l;
  const unsigned int L_ke       = params->faest_param.Lke;
  const unsigned int S_ke       = params->faest_param.Ske;
  const unsigned int num_rounds = params->faest_param.R;

  // uint8_t* w           = malloc((l + 7) / 8);
  uint8_t* const w_out = w;

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

  if (!L_ke) {
    // switch input and key for EM
    const uint8_t* tmp = key;
    key                = in;
    in                 = tmp;
  }

  // Step 3
  aes_round_keys_t round_keys;
  switch (lambda) {
  case 256:
    if (block_words == RIJNDAEL_BLOCK_WORDS_256) {
      rijndael256_init_round_keys(&round_keys, key);
    } else {
      aes256_init_round_keys(&round_keys, key);
    }
    break;
  case 192:
    if (block_words == RIJNDAEL_BLOCK_WORDS_192) {
      rijndael192_init_round_keys(&round_keys, key);
    } else {
      aes192_init_round_keys(&round_keys, key);
    }
    break;
  default:
    aes128_init_round_keys(&round_keys, key);
    break;
  }

  // Step 4
  if (L_ke > 0) {
    // Key schedule constraints only needed for normal AES, not EM variant.
    for (unsigned int i = 0; i != params->faest_param.Nwd; ++i) {
      memcpy(w, round_keys.round_keys[i / 4][i % 4], sizeof(aes_word_t));
      w += sizeof(aes_word_t);
    }

    for (unsigned int j = 0, ik = params->faest_param.Nwd; j < S_ke / 4; ++j) {
      memcpy(w, round_keys.round_keys[ik / 4][ik % 4], sizeof(aes_word_t));
      w += sizeof(aes_word_t);
      ik += lambda == 192 ? 6 : 4;
    }
  } else {
    // saving the OWF key to the extended witness
    memcpy(w, in, lambda / 8);
    w += lambda / 8;
  }

  // Step 10
  for (unsigned b = 0; b < beta; ++b, in += sizeof(aes_word_t) * block_words) {
    // Step 12
    aes_block_t state;
    load_state(state, in, block_words);

    // Step 13
    add_round_key(0, state, &round_keys, block_words);

    for (unsigned int round = 1; round < num_rounds; ++round) {
      // Step 15
      sub_bytes(state, block_words);
      // Step 16
      shift_row(state, block_words);
      // Step 17
      store_state(w, state, block_words);
      w += sizeof(aes_word_t) * block_words;
      // Step 18
      mix_column(state, block_words);
      // Step 19
      add_round_key(round, state, &round_keys, block_words);
    }
    // last round is not commited to, so not computed
  }

  return w_out;
}

// #######################
//  Masked implementation
// #######################

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
