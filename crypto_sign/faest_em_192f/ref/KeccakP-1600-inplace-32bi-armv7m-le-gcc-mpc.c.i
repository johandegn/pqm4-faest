/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.

Implementation by the designers, hereby denoted as "the implementer".

For more information, feedback or questions, please refer to the Keccak Team website:
https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/

---

This file implements Keccak-p[1600] in a SnP-compatible way.
Please refer to SnP-documentation.h for more details.

This implementation comes with KeccakP-1600-SnP.h in the same folder.
Please refer to LowLevel.build for the exact list of other files it must be combined with.
*/

#include "randomness.h"

//KeccakP1600_stateSizeInByts actually represents multiple states of this size
#define KeccakP1600_stateSplit         (200)

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "brg_endian.h"
#ifdef KeccakReference
#include "displayIntermediateValues.h"
#endif
#include "KeccakP-1600-SnP-mpc.h"
#include "KeccakP-1600-SnP-armv7m.h"
#include "KeccakP-1600-SnP-Relaned.h"

typedef uint64_t tKeccakLane;

#define maxNrRounds 24
#define nrLanes 25
#define index(x, y) (((x)%5)+5*((y)%5))
#define ROL32(a, offset) ((((uint32_t)a) << (offset)) ^ (((uint32_t)a) >> (32-(offset))))

/* ---------------------------------------------------------------- */

/*
State is 400 bytes, first 200 are "share 0" second 200 bytes are
"share 1"
*/
void KeccakP1600MPC_Initialize(void *state)
{
#ifdef KECCAK_MASK_PUBLIC
    rand_mask(state, 1600/8);
	memcpy(state + 1600/8, state, 1600/8);
#else
    memset(state, 0, 2*1600/8);
#endif
}

/* ---------------------------------------------------------------- */

void KeccakP1600MPC_AddByte(void *state, unsigned char *data, unsigned int offset, unsigned int dataCount)
{
    #ifdef DEBUG
    assert(offset < 200);
    #endif
    assert(dataCount == 2);

	for (unsigned int i = 0; i < dataCount; i++) {
		KeccakP1600_AddByte(state + i*KeccakP1600_stateSplit, data[i], offset);
	}
}

/* ---------------------------------------------------------------- */


void KeccakP1600MPC_AddBytes(void *state, const unsigned char **data, unsigned int offset, unsigned int dataCount, unsigned int dataByteLen)
{
    #ifdef DEBUG
    assert(offset < 200);
    assert(offset+dataByteLen <= 200);
    #endif

    assert(dataCount == 2);

	for (unsigned int i = 0; i < dataCount; i++) {
		KeccakP1600_AddBytes(state + i*KeccakP1600_stateSplit, data[i], offset, dataByteLen);
	}
}

/* ---------------------------------------------------------------- */

void KeccakP1600MPC_OverwriteBytes(void *state, const unsigned char **data, unsigned int offset, unsigned int dataCount, unsigned int dataByteLen)
{
    #ifdef DEBUG
    assert(offset < 200);
    assert(offset+dataByteLen <= 200);
    #endif
    assert( dataCount == 2);

	for (unsigned int i = 0; i < dataCount; i++) {
		KeccakP1600_OverwriteBytes(state + i*KeccakP1600_stateSplit, data[i], offset, dataByteLen);
	}
}

/* ---------------------------------------------------------------- */

//??? what is this used for ???
void KeccakP1600MPC_OverwriteWithZeroes(void *state, unsigned int byteCount)
{
    #ifdef DEBUG
    printf("KeccakP1600_OverwriteWithZeroes was called\n");
    assert(byteCount <= 200);
    #endif

	for (unsigned int i = 0; i < KeccakP1600MPC_shareCount; i++) {
		KeccakP1600_OverwriteWithZeroes(state + i*KeccakP1600_stateSplit, byteCount);
	}
}

/* ---------------------------------------------------------------- */

//static void fromBytesToWords(tKeccakLane *stateAsWords, const unsigned char *state);
//static void fromWordsToBytes(unsigned char *state, const tKeccakLane *stateAsWords);

void KeccakP1600MPC_ExtractBytes(const void *state, unsigned char **data, unsigned int offset, unsigned int dataCount, unsigned int dataByteLen)
{
    #ifdef DEBUG
    assert(offset < 200);
    assert(offset+dataByteLen <= 200);
    #endif

    assert(dataCount == 2);

    /*
        Copy from the respective state into the respective output data array.
        The caller can take care of reassembling.
    */
	for (unsigned int i = 0; i < dataCount; i++) {
		KeccakP1600_ExtractBytes(state + i*KeccakP1600_stateSplit, data[i], offset, dataByteLen);
	}
}

/* ---------------------------------------------------------------- */

void KeccakP1600MPC_ExtractAndAddBytes(const void *state, const unsigned char **inputs, unsigned char **outputs, unsigned int offset, unsigned int shareCount, unsigned int byteLen)
{
    #ifdef DEBUG
    assert(offset < 200);
    assert(offset+byteLen <= 200);
    #endif

    assert(shareCount == 2);

	for (unsigned int i = 0; i < shareCount; i++) {
		KeccakP1600_ExtractAndAddBytes(state + i*KeccakP1600_stateSplit, inputs[i], outputs[i], offset, byteLen);
	}
}

/* ---------------------------------------------------------------- */

void KeccakP1600MPCHalf_Permute_24rounds(void *state)
{
	KeccakP1600MPC_Permute_12rounds(state);
	for (unsigned int i = 0; i < KeccakP1600_stateSizeInBytes/sizeof(uint32_t); i++) {
		((uint32_t *)state)[i] ^= ((uint32_t *)state)[i + KeccakP1600_stateSizeInBytes/sizeof(uint32_t)];
		((uint32_t *)state)[i + KeccakP1600_stateSizeInBytes/sizeof(uint32_t)] = 0;
	}
	KeccakP1600_Permute_Nrounds(state, 12);
}
