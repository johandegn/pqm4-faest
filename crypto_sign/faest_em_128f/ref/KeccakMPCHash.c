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

#include <string.h>
#include <stdlib.h>
#include "KeccakMPCHash.h"

#ifndef KECCAK_MASK_NONE

/* ---------------------------------------------------------------- */

MPCHashReturn KeccakMPC_HashInitialize(KeccakMPC_HashInstance *instance, unsigned int rate, unsigned int capacity, unsigned int hashbitlen, unsigned char delimitedSuffix, unsigned int shareCount )
{
    MPCHashReturn result;
    unsigned int i;

    if (delimitedSuffix == 0)
        return KECCAKMPC_FAIL;
    result = (MPCHashReturn)KeccakMaskedWidth1600_MaskedSpongeInitialize(&instance->sponge, rate, capacity);
    if (result != KECCAKMPC_SUCCESS)
        return result;
    instance->fixedOutputLength = hashbitlen;
    instance->delimitedSuffixShares = malloc(shareCount);
    instance->shareCount = shareCount;

    //first share is inited to non null value so that the first share has a correct padding
    //while the other shares will only store the payload bytes without any padding
    instance->delimitedSuffixShares[0] = delimitedSuffix;
    for( i=1; i<shareCount; ++i) {
        instance->delimitedSuffixShares[i] = 0;
    }
    return KECCAKMPC_SUCCESS;
}

/* ---------------------------------------------------------------- */

MPCHashReturn KeccakMPC_HashUpdate(KeccakMPC_HashInstance *instance, BitSequence **data, BitLength databitlen)
{
    if ((databitlen % 8) == 0)
        return (MPCHashReturn)KeccakMaskedWidth1600_MaskedSpongeAbsorb(&instance->sponge, data, instance->shareCount, databitlen/8);
    else {
        MPCHashReturn ret = (MPCHashReturn)KeccakMaskedWidth1600_MaskedSpongeAbsorb(&instance->sponge, data, instance->shareCount, databitlen/8);
        if (ret == KECCAKMPC_SUCCESS) {
            BitLength i;
            BitSequence *lastByteShares = malloc(instance->shareCount);
            unsigned short *delimitedLastBytes = malloc(sizeof(unsigned short) * instance->shareCount);
            memset(delimitedLastBytes, 0, sizeof(unsigned short) * instance->shareCount);
            /* The last partial byte is assumed to be aligned on the least significant bits */
            for(i = 0; i < instance->shareCount; ++i) {
                lastByteShares[i] = data[i][databitlen/8];
            }
            /* Concatenate the last few bits provided here with those of the suffix */
            /*
             * (databitLen % 8) -1 : mask that gets all the trailing (over the last full byte) bits from lastByte; at most the lowest 7 bits (otherwise we would not be in this case)
             * then we shift delimitedSuffix so is does not overlap with theese bytes, so basically
             * delmitedLastBytes = (MSB)  delimitedSuffix concat "trailing bits from last byte"  (LSB)
             * delimitedLastBytes is at least two bytes large
            */
            for(i = 0; i < instance->shareCount; ++i) {
                delimitedLastBytes[i] = (unsigned short)((unsigned short)(lastByteShares[i] & ((1 << (databitlen % 8)) - 1))
                | ((unsigned short)instance->delimitedSuffixShares[i] << (databitlen % 8)));
            }

            if ((delimitedLastBytes[0] & 0xFF00) == 0x0000) { //no more than one byte, store
                for(i = 0; i < instance->shareCount; ++i) {
                    instance->delimitedSuffixShares[i] = delimitedLastBytes[i] & 0xFF;
                }
            }
            else { //got more than one byte => process the full byte
                unsigned char **sharedLastByte = malloc( sizeof(unsigned char**) * instance->shareCount);
                unsigned char *slab = malloc(instance->shareCount);
                for(i = 0; i < instance->shareCount; ++i) {
                    sharedLastByte[i] = slab+i; //ok because byte objects
                    sharedLastByte[i][0] = delimitedLastBytes[i] & 0xFF;
                    instance->delimitedSuffixShares[i] = (delimitedLastBytes[i] >> 8) & 0xFF;
                }
                ret = (MPCHashReturn)KeccakMaskedWidth1600_MaskedSpongeAbsorb(&instance->sponge, sharedLastByte, instance->shareCount, 1);
                free(slab);
                free(sharedLastByte);
            }
            free(lastByteShares);
            free(delimitedLastBytes);
        }
        return ret;
    }
}

/* ---------------------------------------------------------------- */

MPCHashReturn KeccakMPC_HashFinal(KeccakMPC_HashInstance *instance, BitSequence **hashvals)
{
    MPCHashReturn ret = (MPCHashReturn)KeccakMaskedWidth1600_MaskedSpongeAbsorbLastFewBits(&instance->sponge, instance->delimitedSuffixShares, instance->shareCount);
    if (ret == KECCAKMPC_SUCCESS) {
        free(instance->delimitedSuffixShares); //cannot return to absorbing from here. thus this array is no longer needed. as there is no explicit destory function it's a good point to deallocate
        return (MPCHashReturn)KeccakMaskedWidth1600_MaskedSpongeSqueeze(&instance->sponge, hashvals, instance->shareCount, instance->fixedOutputLength/8);
    }
    else
        return ret;
}

/* ---------------------------------------------------------------- */

MPCHashReturn KeccakMPC_HashSqueeze(KeccakMPC_HashInstance *instance, BitSequence **data, BitLength databitlen)
{
    if ((databitlen % 8) != 0)
        return KECCAKMPC_FAIL;
    return (MPCHashReturn)KeccakMaskedWidth1600_MaskedSpongeSqueeze(&instance->sponge, data, instance->shareCount, databitlen/8);
}

#endif