/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

Keccak, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.

Implementation by the designers, hereby denoted as "the implementer".

For more information, feedback or questions, please refer to the Keccak Team website:
https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "KeccakMPCSponge.h"

#ifdef KeccakReference
    #include "displayIntermediateValues.h"
#endif

#ifdef XKCP_has_MaskedKeccakP1600
	#include "KeccakP-1600-SnP.h"

    #define prefix KeccakMaskedWidth1600
    #define SnPMPC KeccakP1600MPC
    #define SnPMPC_width 1600
#if defined(KECCAK_MASK_HALF) && !defined(KECCAK_REDUCED_ROUNDS)
    #define SnPMPC_Permute KeccakP1600MPCHalf_Permute_24rounds
#else
	#define SnPMPC_Permute KeccakP1600MPC_Permute_24rounds
#endif
#ifndef KECCAK_MASK_NONE
    #include "KeccakMPCSponge.inc"
    #undef prefix
    #undef SnPMPC
    #undef SnPMPC_width
    #undef SnPMPC_Permute
#endif
#endif
