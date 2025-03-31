#include "shuffle.h"
#include "fields.h"

//(Fisher-Yates shuffle)
void shuffle_16(uint8_t arr[16]) {
    // Manual loop unrolling to make compiler 
    // produce constant time code

    uint8_t r = bf8_rand();

    int i = 15;
    int j = (r & 0x0f) % (i+1);
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 14;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    r = bf8_rand();

    i = 13;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 12;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    r = bf8_rand();

    i = 11;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 10;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    r = bf8_rand();

    i = 9;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 8;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    r = bf8_rand();

    i = 7;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 6;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    r = bf8_rand();

    i = 5;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 4;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    r = bf8_rand();

    i = 3;
    j = (r & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 2;
    j = ((r >> 2) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 1;
    j = ((r >> 4) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

void shuffle_4(uint8_t arr[4]){
    uint8_t r = bf8_rand();
    int i = 3;
    int j = (r & 0x03) % (i+1);
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 2;
    j = ((r >> 2) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    i = 1;
    j = ((r >> 4) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

void shuffle_both_11(const uint8_t* arr[11], bf128_t* arr_out[11]){
    uint8_t r = bf8_rand();

    int i = 10;
    int j = r % (i+1);
    const uint8_t* temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    bf128_t* temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    r = bf8_rand();

    i = 9;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 8;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    r = bf8_rand();

    i = 7;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 6;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    r = bf8_rand();

    i = 5;
    j = (r & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 4;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    r = bf8_rand();

    i = 3;
    j = (r & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 2;
    j = ((r >> 2) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 1;
    j = ((r >> 4) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;
}

void shuffle_both_6(const uint8_t* arr[6], bf128_t* arr_out[6]){
    uint8_t r = bf8_rand();

    int i = 5;
    int j = (r & 0x0f) % (i+1);
    const uint8_t* temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    bf128_t* temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 4;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    r = bf8_rand();

    i = 3;
    j = (r & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 2;
    j = ((r >> 2) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;

    i = 1;
    j = ((r >> 4) & 0x03) % (i+1);
    temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    temp_out = arr_out[i];
    arr_out[i] = arr_out[j];
    arr_out[j] = temp_out;
}

void byte_combine_bits_shuffle_6(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* real_k, const uint8_t* real_v_w, bf128_t* k_out, bf128_t* v_w_out) {
    const uint8_t* in_places[6];
    bf128_t* out_places[6];

    // byte_combine_bits shuffling
    for(unsigned int l = 0; l < 4; l++){
      in_places[l] = decoy_in + l;
      out_places[l] = decoy_out + l;
    }
    in_places[4] = real_k;
    in_places[5] = real_v_w;

    out_places[4] = k_out;
    out_places[5] = v_w_out;

    shuffle_both_6(in_places, out_places);

    *out_places[0] = bf128_byte_combine_bits_sclf(*in_places[0]);
    *out_places[1] = bf128_byte_combine_bits_sclf(*in_places[1]);
    *out_places[2] = bf128_byte_combine_bits_sclf(*in_places[2]);
    *out_places[3] = bf128_byte_combine_bits_sclf(*in_places[3]);
    *out_places[4] = bf128_byte_combine_bits_sclf(*in_places[4]);
    *out_places[5] = bf128_byte_combine_bits_sclf(*in_places[5]);
}

void byte_combine_bits_shuffle_11(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* real_k, const uint8_t* real_v_w, bf128_t* k_out, bf128_t* v_w_out) {
    const uint8_t* in_places[11];
    bf128_t* out_places[11];

    // byte_combine_bits shuffling
    for(unsigned int l = 0; l < 9; l++){
      in_places[l] = decoy_in + l;
      out_places[l] = decoy_out + l;
    }
    in_places[9] = real_k;
    in_places[10] = real_v_w;

    out_places[9] = k_out;
    out_places[10] = v_w_out;

    shuffle_both_11(in_places, out_places);

    *out_places[0] = bf128_byte_combine_bits_sclf(*in_places[0]);
    *out_places[1] = bf128_byte_combine_bits_sclf(*in_places[1]);
    *out_places[2] = bf128_byte_combine_bits_sclf(*in_places[2]);
    *out_places[3] = bf128_byte_combine_bits_sclf(*in_places[3]);
    *out_places[4] = bf128_byte_combine_bits_sclf(*in_places[4]);
    *out_places[5] = bf128_byte_combine_bits_sclf(*in_places[5]);
    *out_places[6] = bf128_byte_combine_bits_sclf(*in_places[6]);
    *out_places[7] = bf128_byte_combine_bits_sclf(*in_places[7]);
    *out_places[8] = bf128_byte_combine_bits_sclf(*in_places[8]);
    *out_places[9] = bf128_byte_combine_bits_sclf(*in_places[9]);
    *out_places[10] = bf128_byte_combine_bits_sclf(*in_places[10]);
}