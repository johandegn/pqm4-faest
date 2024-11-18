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

//KeccakP1600_stateSizeInByts actually represents multiple states of this size
#define KeccakP1600_stateSplit         200

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "brg_endian.h"
#ifdef KeccakReference
#include "displayIntermediateValues.h"
#endif
#include "KeccakP-1600-SnP.h"

typedef uint64_t tKeccakLane;

#define maxNrRounds 24
#define nrLanes 25
#define index(x, y) (((x)%5)+5*((y)%5))

#include "randomness.h"

static const tKeccakLane KeccakRoundConstants[maxNrRounds] =
{
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808a,
    0x8000000080008000,
    0x000000000000808b,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008a,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000a,
    0x000000008000808b,
    0x800000000000008b,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800a,
    0x800000008000000a,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
};

static const unsigned int KeccakRhoOffsets[nrLanes] =
{
     0,  1, 62, 28, 27, 36, 44,  6, 55, 20,  3, 10, 43, 25, 39, 41, 45, 15, 21,  8, 18,  2, 61, 56, 14
};

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
		KeccakP1600_AddByte((uint8_t *)state + i*KeccakP1600_stateSplit, data[i], offset);
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
		KeccakP1600_AddBytes((uint8_t *)state + i*KeccakP1600_stateSplit, data[i], offset, dataByteLen);
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
		KeccakP1600_OverwriteBytes((uint8_t *)state + i*KeccakP1600_stateSplit, data[i], offset, dataByteLen);
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
		KeccakP1600_OverwriteWithZeroes((uint8_t *)state + i*KeccakP1600_stateSplit, byteCount);
	}
}

/* ---------------------------------------------------------------- */

//static void fromBytesToWords(tKeccakLane *stateAsWords, const unsigned char *state);
//static void fromWordsToBytes(unsigned char *state, const tKeccakLane *stateAsWords);
void KeccakP1600MPCOnWords(tKeccakLane *state, unsigned int nrRounds);
void KeccakP1600MPCRound(tKeccakLane *state, unsigned int indexRound);
static void theta(tKeccakLane *A);
static void rho(tKeccakLane *A);
static void pi(tKeccakLane *A);
static void mpc_chi(tKeccakLane *state0, tKeccakLane * state1);
static void iota(tKeccakLane *A, unsigned int indexRound);

void KeccakP1600MPC_Permute_Nrounds(void *state, unsigned int nrounds)
{
#if (PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN)
    tKeccakLane stateAsWords[1600/64];
#endif

#ifdef KeccakReference
    displaySharedStateAsBytes(1, "Input of permutation", state, KeccakP1600_stateSplit, 2);
#endif
#if (PLATFORM_BYTE_ORDER == IS_LITTLE_ENDIAN)
    KeccakP1600MPCOnWords((tKeccakLane*)state, nrounds);
#else
    fromBytesToWords(stateAsWords, (const unsigned char *)state);
    KeccakP1600SplitOnWords(stateAsWords, nrounds);
    fromWordsToBytes((unsigned char *)state, stateAsWords);
#endif
#ifdef KeccakReference
    displaySharedStateAsBytes(1, "State after permutation", state, KeccakP1600_stateSplit, 2);
#endif
}

void KeccakP1600MPC_Permute_12rounds(void *state)
{
#if (PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN)
    tKeccakLane stateAsWords[1600/64];
#endif

#ifdef KeccakReference
    displaySharedStateAsBytes(1, "Input of permutation", state, KeccakP1600_stateSplit, 2);
#endif
#if (PLATFORM_BYTE_ORDER == IS_LITTLE_ENDIAN)
    KeccakP1600MPCOnWords((tKeccakLane*)state, 12);
#else
    fromBytesToWords(stateAsWords, (const unsigned char *)state);
    KeccakP1600MPCOnWords(stateAsWords, 12);
    fromWordsToBytes((unsigned char *)state, stateAsWords);
#endif
#ifdef KeccakReference
    displaySharedStateAsBytes(1, "State after permutation", state, KeccakP1600_stateSplit, 2);
#endif
}

void KeccakP1600MPC_Permute_24rounds(void *state)
{
#if (PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN)
    tKeccakLane stateAsWords[1600/64];
#endif

#ifdef KeccakReference
    displaySharedStateAsBytes(1, "Input of permutation", state, KeccakP1600_stateSplit, 2);;
#endif
#if (PLATFORM_BYTE_ORDER == IS_LITTLE_ENDIAN)
    KeccakP1600MPCOnWords((tKeccakLane*)state, 24);
#else
    fromBytesToWords(stateAsWords, (const unsigned char *)state);
    KeccakP1600MPCOnWords(stateAsWords, 24);
    fromWordsToBytes((unsigned char *)state, stateAsWords);
#endif
#ifdef KeccakReference
    displaySharedStateAsBytes(1, "State after permutation", state, KeccakP1600_stateSplit, 2);
#endif
}

#ifdef KECCAK_MASK_HALF

void KeccakP1600MPCHalf_Permute_24rounds(void *state)
{
#if (PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN)
    tKeccakLane stateAsWords[2*1600/64];
#endif

#ifdef KeccakReference
    displaySharedStateAsBytes(1, "Input of permutation", state, KeccakP1600_stateSplit, 2);
#endif
#if (PLATFORM_BYTE_ORDER == IS_LITTLE_ENDIAN)
	/* Run first 12 rounds in 2 shares. */
	for(int i = 0; i < 12; i++) {
		KeccakP1600MPCRound(state, i);
	}
	/* Combine shares and finish in single state. */
	for (int i = 0; i < 1600/64; i++) {
		((tKeccakLane *)state)[i] ^= ((tKeccakLane *)state)[i + 1600/64];
		((tKeccakLane *)state)[i + 1600/64] = 0;
	}
	KeccakP1600_Permute_12rounds(state);
#else
    fromBytesToWords(stateAsWords, (const unsigned char *)state);
	/* Run first 12 rounds in 2 shares. */
	for(int i = 0; i < 12; i++) {
		KeccakP1600MPCRound(stateAsWords, i);
	}
	/* Combine shares and finish in single state. */
	for (int i = 0; i < 1600/64; i++) {
		((tKeccakLane *)stateAsWords)[i] ^= ((tKeccakLane *)state)[i + 1600/64];
		((tKeccakLane *)stateAsWords)[i + 1600/64] = 0;
	}
	KeccakP1600_Permute_12rounds(stateAsWords);
    fromWordsToBytes((unsigned char *)state, stateAsWords);
#endif
#ifdef KeccakReference
    displaySharedStateAsBytes(1, "State after permutation", state, KeccakP1600_stateSplit, 2);
#endif
}

#endif

/*
static void fromBytesToWords(tKeccakLane *stateAsWords, const unsigned char *state)
{
    unsigned int i, j;

    for(i=0; i<nrLanes; i++) {
        stateAsWords[i] = 0;
        for(j=0; j<(64/8); j++)
            stateAsWords[i] |= (tKeccakLane)(state[i*(64/8)+j]) << (8*j);
    }
}
*/

/*
static void fromWordsToBytes(unsigned char *state, const tKeccakLane *stateAsWords)
{
    unsigned int i, j;

    for(i=0; i<nrLanes; i++)
        for(j=0; j<(64/8); j++)
            state[i*(64/8)+j] = (unsigned char)((stateAsWords[i] >> (8*j)) & 0xFF);
}
*/

void KeccakP1600MPCOnWords(tKeccakLane *state, unsigned int nrRounds)
{
    unsigned int i;

#ifdef KeccakReference
    displaySharedStateAsLanes(3, "Same, with lanes as 64-bit words", state, KeccakP1600_stateSplit, 2);
#endif

    for(i=(maxNrRounds-nrRounds); i<maxNrRounds; i++)
        KeccakP1600MPCRound(state, i);
}

void KeccakP1600MPCRound(tKeccakLane *state, unsigned int indexRound)
{
#ifdef KeccakReference
    displayRoundNumber(3, indexRound);
#endif

	for (unsigned int i = 0; i < KeccakP1600MPC_shareCount; i++) {
		theta(state + i*nrLanes);
	}
#ifdef KeccakReference
    displaySharedStateAsLanes(3, "State after theta", state, KeccakP1600_stateSplit, 2);
#endif

	for (unsigned int i = 0; i < KeccakP1600MPC_shareCount; i++) {
		rho(state + i*nrLanes);
	}
#ifdef KeccakReference
    displaySharedStateAsLanes(3, "State 0 after rho", state, KeccakP1600_stateSplit, 2);
#endif

	for (unsigned int i = 0; i < KeccakP1600MPC_shareCount; i++) {
		pi(state + i*nrLanes);
	}
#ifdef KeccakReference
    displaySharedStateAsLanes(3, "State 0 after pi", state, KeccakP1600_stateSplit, 2);
#endif

    mpc_chi(state,state+nrLanes);
#ifdef KeccakReference
    displaySharedStateAsLanes(3, "State 0 after chi", state,KeccakP1600_stateSplit, 2);
#endif
    //This must *not* be done on both states!
    iota(state, indexRound);
#ifdef KeccakReference
    displaySharedStateAsLanes(3, "After iota", state, KeccakP1600_stateSplit, 2);
#endif
}

#ifndef ROL64
#define ROL64(a, offset) ((offset != 0) ? ((((tKeccakLane)a) << offset) ^ (((tKeccakLane)a) >> (64-offset))) : a)
#endif

static void theta(tKeccakLane *A)
{
    unsigned int x, y;
    tKeccakLane C[5], D[5];

    for(x=0; x<5; x++) {
        C[x] = 0;
        for(y=0; y<5; y++)
            C[x] ^= A[index(x, y)];
    }
    for(x=0; x<5; x++)
        D[x] = ROL64(C[(x+1)%5], 1) ^ C[(x+4)%5];
    for(x=0; x<5; x++)
        for(y=0; y<5; y++)
            A[index(x, y)] ^= D[x];
}

static void rho(tKeccakLane *A)
{
    unsigned int x, y;

    for(x=0; x<5; x++) for(y=0; y<5; y++)
        A[index(x, y)] = ROL64(A[index(x, y)], KeccakRhoOffsets[index(x, y)]);
}

static void pi(tKeccakLane *A)
{
    unsigned int x, y;
    tKeccakLane tempA[25];

    for(x=0; x<5; x++) for(y=0; y<5; y++)
        tempA[index(x, y)] = A[index(x, y)];
    for(x=0; x<5; x++) for(y=0; y<5; y++)
        A[index(0*x+1*y, 2*x+3*y)] = tempA[index(x, y)];
}

/*
static void chi(tKeccakLane *A)
{
    unsigned int x, y;
    tKeccakLane C[5];

    for(y=0; y<5; y++) {
        for(x=0; x<5; x++)
            //perform computation on all rows of a lane at once, since we have 64bit registers
            C[x] = A[index(x, y)] ^ ((~A[index(x+1, y)]) & A[index(x+2, y)]);
        for(x=0; x<5; x++) //propgate update only after we are done with one y dimension due two the dependencies (x,x+1,x+2) in the x dimension
            A[index(x, y)] = C[x];
    }
}
*/

/*
*Assumes state0 and state1 are both 200 bytes
Uses the mpc scheme from Bertoni's Paper
*/
static void mpc_chi(tKeccakLane *state0, tKeccakLane *state1) {
    unsigned int x,y;

    //buffer for updated lane in state 0
    tKeccakLane C0[5];
    //buffer for updated lane in state 1
    tKeccakLane C1[5];

#ifdef KECCAK_SNI_SECURE
    uint64_t _mask;
    for(y=0; y<5; y++) {
        for(x=0; x<5; x++) {
            rand_mask((uint8_t*)&_mask, sizeof(uint64_t));
            state0[index(x,y)] ^= _mask;
            state1[index(x,y)] ^= _mask;
        }
    }
#endif

    //Coming from paper notation: a <-> state0 , b <-> state1
    for(y=0; y<5; y++) {
#if defined(KECCAK_DOM_SECURE) || defined(KECCAK_SNI_SECURE)
        uint64_t mask;
        for(x=0; x<5; x++) {
            rand_mask((uint8_t*)&mask, sizeof(uint64_t));
            C0[x] = (state0[index(x+1,y)] & state1[index(x+2,y)]) ^ mask;
            C1[x] = (state1[index(x+1,y)] & state0[index(x+2,y)]) ^ mask;
        }
        for(x=0; x<5; x++) {
            C0[x] ^= (~state0[index(x+1,y)]) & state0[index(x+2,y)];
            C1[x] ^= (~state1[index(x+1,y)]) & state1[index(x+2,y)];
        }
        for(x=0; x<5; x++) {
            C0[x] ^= state0[index(x,y)];
            C1[x] ^= state1[index(x,y)];
        }
#else
        for(x=0; x<5; x++) {
            C0[x] = state0[index(x,y)] ^ ( (~state0[index(x+1,y)]) & state0[index(x+2,y)]);
            C1[x] = state1[index(x,y)] ^ ( (~state1[index(x+1,y)]) & state1[index(x+2,y)]);
        }
        for(x=0; x<5; x++) {
            C0[x] ^= state0[index(x+1,y)] & state1[index(x+2,y)];
            C1[x] ^= state1[index(x+1,y)] & state0[index(x+2,y)];
        }
#endif
        for(x=0; x<5; x++) {
            state0[index(x,y)] = C0[x];
            state1[index(x,y)] = C1[x];
        }
    }
}

static void iota(tKeccakLane *A, unsigned int indexRound)
{
    A[index(0, 0)] ^= KeccakRoundConstants[indexRound];
}

/* ---------------------------------------------------------------- */


void KeccakP1600MPC_ExtractBytes(const void *state, unsigned char **data, unsigned int offset, unsigned int dataCount, unsigned int dataByteLen)
{
    #ifdef DEBUG
    assert(offset < 200);
    assert(offset+dataByteLen <= 200);
    #endif

    assert(dataCount == 2);

    /*
        Copy from the respective state into the respective output data array.
        The caller can take care of reassembling
    */
	for (unsigned int i = 0; i < dataCount; i++) {
		KeccakP1600_ExtractBytes((uint8_t *)state + i*KeccakP1600_stateSplit, data[i], offset, dataByteLen);
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
		KeccakP1600_ExtractAndAddBytes((uint8_t *)state + i*KeccakP1600_stateSplit, inputs[i], outputs[i], offset, byteLen);
	}
}

/* ---------------------------------------------------------------- */

void KeccakP1600MPC_DisplayRoundConstants(FILE *f)
{
    unsigned int i;

    for(i=0; i<maxNrRounds; i++) {
        fprintf(f, "RC[%02i][0][0] = ", i);
        fprintf(f, "%08X", (unsigned int)(KeccakRoundConstants[i] >> 32));
        fprintf(f, "%08X", (unsigned int)(KeccakRoundConstants[i] & 0xFFFFFFFFULL));
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

void KeccakP1600MPC_DisplayRhoOffsets(FILE *f)
{
    unsigned int x, y;

    for(y=0; y<5; y++) for(x=0; x<5; x++) {
        fprintf(f, "RhoOffset[%i][%i] = ", x, y);
        fprintf(f, "%2i", KeccakRhoOffsets[index(x, y)]);
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}
