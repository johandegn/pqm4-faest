#include <stdint.h>
#include "fields.h"

void shuffle_16(uint32_t arr[16], uint32_t masks[16]);
void byte_combine_bits_shuffle_6(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* real_k, const uint8_t* real_v_w, bf128_t* k_out, bf128_t* v_w_out);

void byte_combine_bits_shuffle_32_enc(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* k, const uint8_t* v_w, bf128_t* k_out, bf128_t* v_w_out, unsigned int ix);
void byte_combine_bits_shuffle_32_key(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* k, const uint8_t* v_w, bf128_t* k_out, bf128_t* v_w_out, unsigned int iwd);