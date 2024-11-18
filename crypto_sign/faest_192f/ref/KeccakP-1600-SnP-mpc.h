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

Please refer to SnP-documentation.h for more details.
*/

#ifndef _KeccakP_1600_MPC_SnP_h_
#define _KeccakP_1600_MPC_SnP_h_

#define KeccakP1600MPC_implementation      "64-bit reference mpc with masked input implementation"
//doubled the state size as whe want two internal states
#define KeccakP1600MPC_stateSizeInBytes    400
#define KeccakP1600MPC_stateAlignment      8
//communicate to sponge construction how many shares this implementation allows
#define KeccakP1600MPC_shareCount            2

#ifdef KeccakReference
void KeccakP1600MPC_StaticInitialize( void );
#else
#define KeccakP1600MPC_StaticInitialize()
#endif
void KeccakP1600MPC_Initialize(void *state);
void KeccakP1600MPC_AddByte(void *state, unsigned char *data, unsigned int offset, unsigned int dataCount);
void KeccakP1600MPC_AddBytes(void *state, const unsigned char **data, unsigned int offset, unsigned int dataCount, unsigned int dataByteLen);
void KeccakP1600MPC_OverwriteBytes(void *state, const unsigned char **data, unsigned int offset, unsigned int dataCount, unsigned int dataByteLen);
void KeccakP1600MPC_OverwriteWithZeroes(void *state, unsigned int byteCount);
void KeccakP1600MPC_Permute_Nrounds(void *state, unsigned int nrounds);
void KeccakP1600MPC_Permute_12rounds(void *state);
void KeccakP1600MPC_Permute_24rounds(void *state);
void KeccakP1600MPCHalf_Permute_24rounds(void *state);
void KeccakP1600MPC_ExtractBytes(const void *state, unsigned char **data, unsigned int offset, unsigned int dataCount, unsigned int dataByteLen);
void KeccakP1600MPC_ExtractAndAddBytes(const void *state, const unsigned char **inputs, unsigned char **outputs, unsigned int offset, unsigned int shareCount, unsigned int byteLen);

#endif
