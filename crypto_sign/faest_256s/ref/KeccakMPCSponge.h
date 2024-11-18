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

#ifndef _KeccakMPCSponge_h_
#define _KeccakMPCSponge_h_

/* For the documentation, please follow the link: */
/* #include "KeccakSponge-documentation.h" */

#include <string.h>
#include "align.h"
#include "config.h"

#define XKCP_DeclareMaskedSpongeStructure(prefix, size, alignment) \
    ALIGN(alignment) typedef struct prefix##_MaskedSpongeInstanceStruct { \
        unsigned char state[size]; \
        unsigned int rate; \
        unsigned int byteIOIndex; \
        int squeezing; \
    } prefix##_MaskedSpongeInstance;

#define XKCP_DeclareMaskedSpongeFunctions(prefix) \
    int prefix##_MaskedSponge(unsigned int rate, unsigned int capacity,  unsigned char **inputs, size_t inputCount, size_t inputByteLen, unsigned char *suffixes, unsigned char **outputs, size_t outputByteLen); \
    int prefix##_MaskedSpongeInitialize(prefix##_MaskedSpongeInstance *spongeInstance, unsigned int rate, unsigned int capacity); \
    int prefix##_MaskedSpongeAbsorb(prefix##_MaskedSpongeInstance *spongeInstance, unsigned char **inputs, size_t inputCount, size_t inputByteLen); \
    int prefix##_MaskedSpongeAbsorbLastFewBits(prefix##_MaskedSpongeInstance *spongeInstance, unsigned char * delimitedData, size_t delimitedDataCount); \
    int prefix##_MaskedSpongeSqueeze(prefix##_MaskedSpongeInstance *spongeInstance, unsigned char **outputs, size_t outputCount, size_t outputByteLen);


#ifdef XKCP_has_MaskedKeccakP1600
    #include "KeccakP-1600-SnP-mpc.h"
    XKCP_DeclareMaskedSpongeStructure(KeccakMaskedWidth1600, KeccakP1600MPC_stateSizeInBytes, KeccakP1600MPC_stateAlignment)
    XKCP_DeclareMaskedSpongeFunctions(KeccakMaskedWidth1600)
    #define XKCP_has_Sponge_Masked_Keccak_width1600
#else
    #error "debug error"
#endif

#endif
