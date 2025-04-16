#include "shuffle.h"
#include "fields.h"
#include "randomness.h"

//(Fisher-Yates shuffle)
void shuffle_16(uint32_t arr[16], uint32_t masks[16]) {
    // Manual loop unrolling to make compiler 
    // produce constant time code

    uint8_t r = bf8_rand();
    int i = 15;
    int j = (r & 0x0f) % (i+1);
    int temp_i = arr[i];
    int temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 14;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 13;
    j = (r & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 12;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 11;
    j = (r & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 10;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 9;
    j = (r & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 8;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 7;
    j = (r & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 6;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 5;
    j = (r & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 4;
    j = ((r >> 4) & 0x0f) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 3;
    j = (r & 0x03) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 2;
    j = ((r >> 2) & 0x03) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    r = bf8_rand();
    i = 1;
    j = ((r >> 4) & 0x03) % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = ((uint32_t)temp_j ^ masks[i]);

    arr[0] = ((uint32_t)arr[0] ^ masks[0]);
}


void shuffle_both_17_mask(const uint8_t* arr[17], bf128_t* arr_out[17], uint32_t masks[17]){
    uint8_t r = bf8_rand();

    int i = 16;
    int j = r % (i+1);
    const uint8_t* temp_i = arr[i];
    const uint8_t* temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    bf128_t* temp_out_i = arr_out[i];
    bf128_t* temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 15;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 14;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 13;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 12;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 11;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 10;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 9;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 8;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 7;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);
    
    i = 6;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 5;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 4;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 3;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 2;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 1;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    arr[0] = (uint8_t*)((uint32_t)arr[0] ^ masks[0]);
    arr_out[0] = (bf128_t*)((uint32_t)arr_out[0] ^ masks[0]);
}

void shuffle_both_6_mask(const uint8_t* arr[6], bf128_t* arr_out[6], uint32_t masks[6]){
    uint8_t r = bf8_rand();

    int i = 5;
    int j = r % (i+1);
    const uint8_t* temp_i = arr[i];
    const uint8_t* temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    bf128_t* temp_out_i = arr_out[i];
    bf128_t* temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 4;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 3;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 2;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    i = 1;
    r = bf8_rand();
    j = r % (i+1);
    temp_i = arr[i];
    temp_j = arr[j];
    arr[j] = temp_i;
    arr[i] = (uint8_t*)((uint32_t)temp_j ^ masks[i]);
    temp_out_i = arr_out[i];
    temp_out_j = arr_out[j];
    arr_out[j] = temp_out_i;
    arr_out[i] = (bf128_t*)((uint32_t) temp_out_j ^ masks[i]);

    arr[0] = (uint8_t*)((uint32_t)arr[0] ^ masks[0]);
    arr_out[0] = (bf128_t*)((uint32_t)arr_out[0] ^ masks[0]);
}

void shuffle_both(uint8_t* arr[], bf128_t* arr_out[], int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = bf8_rand() % (i + 1);  // Pick a random index from 0 to i
        uint8_t* temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
        bf128_t* temp_out = arr_out[i];
        arr_out[i] = arr_out[j];
        arr_out[j] = temp_out;
    }
}

void byte_combine_bits_shuffle_6(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* real_k, const uint8_t* real_v_w, bf128_t* k_out, bf128_t* v_w_out) {
    const uint8_t* in_places[6];
    bf128_t* out_places[6];

    uint32_t masks[6];
    uint32_t fixed_mask;
    rand_mask((uint8_t*)masks, 6*4);
    rand_mask((uint8_t*)&fixed_mask, 4);

    for(unsigned int l = 0; l < 4; l++){
      in_places[l] = (uint8_t*)((uint32_t)(decoy_in + l)^fixed_mask);
      out_places[l] = (bf128_t*)((uint32_t)(decoy_out + l)^fixed_mask);
    }
    in_places[4] = (uint8_t*)((uint32_t)real_k ^ fixed_mask);
    in_places[5] = (uint8_t*)((uint32_t)real_v_w ^ fixed_mask);

    out_places[4] = (bf128_t*)((uint32_t)k_out ^ fixed_mask);
    out_places[5] = (bf128_t*)((uint32_t)v_w_out ^ fixed_mask);

    shuffle_both_6_mask(in_places, out_places, masks);

    *(bf128_t*)((uint32_t)out_places[0]^fixed_mask^masks[0]) = bf128_byte_combine_bits_sclf(*(uint8_t*)((uint32_t)in_places[0]^fixed_mask^masks[0]));
    *(bf128_t*)((uint32_t)out_places[1]^fixed_mask^masks[1]) = bf128_byte_combine_bits_sclf(*(uint8_t*)((uint32_t)in_places[1]^fixed_mask^masks[1]));
    *(bf128_t*)((uint32_t)out_places[2]^fixed_mask^masks[2]) = bf128_byte_combine_bits_sclf(*(uint8_t*)((uint32_t)in_places[2]^fixed_mask^masks[2]));
    *(bf128_t*)((uint32_t)out_places[3]^fixed_mask^masks[3]) = bf128_byte_combine_bits_sclf(*(uint8_t*)((uint32_t)in_places[3]^fixed_mask^masks[3]));
    *(bf128_t*)((uint32_t)out_places[4]^fixed_mask^masks[4]) = bf128_byte_combine_bits_sclf(*(uint8_t*)((uint32_t)in_places[4]^fixed_mask^masks[4]));
    *(bf128_t*)((uint32_t)out_places[5]^fixed_mask^masks[5]) = bf128_byte_combine_bits_sclf(*(uint8_t*)((uint32_t)in_places[5]^fixed_mask^masks[5]));
}

void byte_combine_bits_shuffle_32_key(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* k, const uint8_t* v_w, bf128_t* k_out, bf128_t* v_w_out, unsigned int iwd) {
    const uint8_t* in_places[9+4+4];
    bf128_t* out_places[9+4+4];

    uint32_t masks[17];
    uint32_t fixed_mask;
    rand_mask((uint8_t*)masks, 17*4);
    rand_mask((uint8_t*)&fixed_mask, 4);

    // byte_combine_bits shuffling
    for(unsigned int l = 0; l < 9; l++){
      in_places[l] = (uint8_t*)((uint32_t)(decoy_in + l)^fixed_mask);
      out_places[l] = (bf128_t*)((uint32_t)(decoy_out + l)^fixed_mask);
    }
    for(unsigned int r = 0; r < 4; r++){
      in_places[9 + r] = (uint8_t*)((uint32_t)(k + (96 + iwd*8 + 8 * r) / 8)^fixed_mask);
      out_places[9 + r] = (bf128_t*)((uint32_t)(k_out + (r + 3) % 4)^fixed_mask);
    }
    for(unsigned int r = 0; r < 4; r++){
      in_places[9 + 4 + r] = (uint8_t*)((uint32_t)(v_w + (8 * r) / 8)^fixed_mask);
      out_places[9 + 4 + r] = (bf128_t*)((uint32_t)(v_w_out + r)^fixed_mask);
    }

    shuffle_both_17_mask(in_places, out_places, masks);

    for(int i = 0; i < 9+4+4; i++){
      *(bf128_t*)((uint32_t)out_places[i]^fixed_mask^masks[i]) = bf128_byte_combine_bits(*(uint8_t*)((uint32_t)in_places[i]^fixed_mask^masks[i]));
    }
}

void byte_combine_bits_shuffle_32_enc(uint8_t* decoy_in, bf128_t* decoy_out, const uint8_t* k, const uint8_t* v_w, bf128_t* k_out, bf128_t* v_w_out, unsigned int ix) {
    const uint8_t* in_places[9+4+4];
    bf128_t* out_places[9+4+4];

    uint32_t masks[17];
    uint32_t fixed_mask;
    rand_mask((uint8_t*)masks, 17*4);
    rand_mask((uint8_t*)&fixed_mask, 4);

    for(unsigned int l = 0; l < 9; l++){
      in_places[l] = (uint8_t*)((uint32_t)(decoy_in + l)^fixed_mask);
      out_places[l] = (bf128_t*)((uint32_t)(decoy_out + l)^fixed_mask);
    }
    for(unsigned int r = 0; r < 4; r++){
      in_places[9 + r] = (uint8_t*)((uint32_t)(k + (ix + 8 * r) / 8)^fixed_mask);
      out_places[9 + r] = (bf128_t*)((uint32_t)(k_out + r)^fixed_mask);
    }
    for(unsigned int r = 0; r < 4; r++){
      in_places[9 + 4 + r] = (uint8_t*)((uint32_t)(v_w + ((ix+128) + 8 * r) / 8)^fixed_mask);
      out_places[9 + 4 + r] = (bf128_t*)((uint32_t)(v_w_out + r)^fixed_mask);
    }

    shuffle_both_17_mask(in_places, out_places, masks);

    for(int i = 0; i < 9+4+4; i++){
      *(bf128_t*)((uint32_t)out_places[i]^fixed_mask^masks[i]) = bf128_byte_combine_bits(*(uint8_t*)((uint32_t)in_places[i]^fixed_mask^masks[i]));
    }
}