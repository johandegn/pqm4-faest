#include <stdint.h>
#include "fields.h"

void shuffle_16(uint8_t arr[16]);
void shuffle_4(uint8_t arr[4]);
void shuffle_both_11(const uint8_t* arr[11], bf128_t* arr_out[11]);
void shuffle_both_6(const uint8_t* arr[6], bf128_t* arr_out[6]);
void byte_combine_bits_shuffle_6(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* real_k, const uint8_t* real_v_w, bf128_t* k_out, bf128_t* v_w_out);
void byte_combine_bits_shuffle_11(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* real_k, const uint8_t* real_v_w, bf128_t* k_out, bf128_t* v_w_out);