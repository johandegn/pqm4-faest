/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_AES_H
#define FAEST_AES_H

#include "macros.h"
#include "instances.h"
#include "fields.h"

#include <stdint.h>
#include <stdlib.h>

FAEST_BEGIN_C_DECL

#define AES_MAX_ROUNDS 14

typedef uint8_t aes_word_t[4];
// round key with 4 (AES) up to 8 (Rijndael-256) units
// TODO: aes_round_key_t[8] should this be 8 ?
typedef aes_word_t aes_round_key_t[8];

// # of rows
#define AES_NR 4

// block with 4 (AES) up to 8 (Rijndael-256) units
typedef aes_word_t aes_block_t[8];

typedef struct {
  aes_round_key_t round_keys[AES_MAX_ROUNDS + 1];
} aes_round_keys_t;

int aes128_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
int aes192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
int aes256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
int rijndael192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
int rijndael256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);

int aes128_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext);
int aes192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext);
int aes256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext);
int rijndael192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                              uint8_t* ciphertext);
int rijndael256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                              uint8_t* ciphertext);

void aes_increment_iv(uint8_t* iv);

uint8_t* aes_extend_witness(const uint8_t* key, const uint8_t* in, const faest_paramset_t* params, uint8_t* w);

int expand_key(aes_round_keys_t* round_keys, const uint8_t* key, unsigned int key_words,
               unsigned int block_words, unsigned int num_rounds);

void prg(const uint8_t* key, const uint8_t* iv, uint8_t* out, unsigned int bits, size_t outlen);

void sub_words(bf8_t* words);
void rot_word(bf8_t* words);
void store_state(uint8_t* dst, aes_block_t state, unsigned int block_words);
void load_state(aes_block_t state, const uint8_t* src, unsigned int block_words);
void mix_column(aes_block_t state, unsigned int block_words);
void shift_row(aes_block_t state, unsigned int block_words);
void add_round_key(unsigned int round, aes_block_t state, const aes_round_keys_t* round_key,
                          unsigned int block_words);
bf8_t affine(bf8_t b);
bf8_t affine_incomplete(bf8_t b);
bf8_t round_constants(bf8_t b);
FAEST_END_C_DECL

#endif
