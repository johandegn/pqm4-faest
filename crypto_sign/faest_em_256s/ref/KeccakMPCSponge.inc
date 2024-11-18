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

#include<alloca.h>
#include<assert.h>
#include<stdlib.h>
#include<stdio.h>

#define JOIN0(a, b)                     a ## b
#define JOIN(a, b)                      JOIN0(a, b)

#define MaskedSponge                          JOIN(prefix, _MaskedSponge)
#define MaskedSpongeInstance                  JOIN(prefix, _MaskedSpongeInstance)
#define MaskedSpongeInitialize                JOIN(prefix, _MaskedSpongeInitialize)
#define MaskedSpongeAbsorb                    JOIN(prefix, _MaskedSpongeAbsorb)
#define MaskedSpongeAbsorbLastFewBits         JOIN(prefix, _MaskedSpongeAbsorbLastFewBits)
#define MaskedSpongeSqueeze                   JOIN(prefix, _MaskedSpongeSqueeze)

//the following code uses the SnP abstraction which is resolved here
//at compile time
#define SnPMPC_stateSizeInBytes            JOIN(SnPMPC, _stateSizeInBytes)
#define SnPMPC_stateAlignment              JOIN(SnPMPC, _stateAlignment)
#define SnPMPC_shareCount                  JOIN(SnPMPC, _shareCount)
#define SnPMPC_StaticInitialize            JOIN(SnPMPC, _StaticInitialize)
#define SnPMPC_Initialize                  JOIN(SnPMPC, _Initialize)
#define SnPMPC_AddByte                     JOIN(SnPMPC, _AddByte)
#define SnPMPC_AddBytes                    JOIN(SnPMPC, _AddBytes)
#define SnPMPC_ExtractBytes                JOIN(SnPMPC, _ExtractBytes)

int MaskedSponge(unsigned int rate, unsigned int capacity, unsigned char **inputs, size_t inputCount, size_t inputByteLen, unsigned char *suffixes, unsigned char **outputs, size_t outputByteLen)
{
    ALIGN(SnPMPC_stateAlignment) unsigned char state[SnPMPC_stateSizeInBytes];
    assert( inputCount == SnPMPC_shareCount );
    unsigned int partialBlock;
    size_t i;
    unsigned char *lastSuffixBits = NULL;


    unsigned int rateInBytes = rate/8;

    if (rate+capacity != SnPMPC_width)
        return 1;
    if ((rate <= 0) || (rate > SnPMPC_width) || ((rate % 8) != 0))
        return 1;

    //0 is an invalid suffix encoding
    //only suffixes[0] contains the padding bits. the other shares only contain payload bits (filled with zeroes)
    if( suffixes[0] == 0) {
        return 1;
    }

    unsigned char *curInputs[SnPMPC_shareCount];
    unsigned char *curOutputs[SnPMPC_shareCount];
    for( i = 0; i < inputCount; ++i ) {
        curInputs[i] = inputs[i];
        curOutputs[i] = outputs[i];
    }

    /* Initialize the state */
    SnPMPC_StaticInitialize();
    SnPMPC_Initialize(state);

    /* First, absorb whole blocks */
#ifdef SnP_FastLoop_Absorb
    if (((rateInBytes % (SnP_width/200)) == 0) && (inputByteLen >= rateInBytes)) {
        /* fast lane: whole lane rate */
        size_t j;
        j = SnPMPC_FastLoop_Absorb(state, rateInBytes/(SnP_width/200), curInputs, inputCount ,inputByteLen );
        //advance inputs by number of processed bytes
         for( i = 0; i < inputCount; ++i) {
            curInputs[i] += j;
        }
        inputByteLen -= j;

    }
#endif

    while(inputByteLen >= (size_t)rateInBytes) {
        #ifdef KeccakReference
        displaySharedBytes(1, "Block to be absorbed", curInputs, inputCount, rateInBytes);
        #endif
        SnPMPC_AddBytes(state, (const unsigned char**)curInputs, 0, inputCount, rateInBytes);
        SnPMPC_Permute(state);

        //advance offset in input array
        for( i = 0; i < inputCount; ++i) {
            curInputs[i] += rateInBytes;
        }
        inputByteLen -= rateInBytes;
    }

    /* Then, absorb what remains */
    partialBlock = (unsigned int)inputByteLen;
    #ifdef KeccakReference
    displaySharedBytes(1, "Block to be absorbed (part)", curInputs, inputCount, partialBlock);
    #endif

    SnPMPC_AddBytes(state, (const unsigned char **)curInputs, 0, inputCount, partialBlock);

    /* Finally, absorb the suffix */
    #ifdef KeccakReference
        displaySharedByte(1, "Block to be absorbed (last few bits + first bit of padding)", suffixes, inputCount);
    #endif

    //Note Luca:
    //We use 10^{j}1 i.e. to add 4 padding bits we would use "1001"
    //According to KeccakSponge-documentation.h the format of suffix (MSB (index 7) .. LSB (index 0))
    //is 0...1payload. Thus after adding suffix, we only need to add a single "1" to complete the padding

    /* Last few bits, whose delimiter coincides with first bit of padding */
    SnPMPC_AddByte(state, suffixes, partialBlock, inputCount);
    /* If the first bit of padding is at position rate-1, we need a whole new block for the second bit of padding */


    if ((suffixes[0] >= 0x80) && (partialBlock == (rateInBytes-1)))
        SnPMPC_Permute(state);
    /* Second bit of padding */
    lastSuffixBits = alloca(inputCount);
    memset(lastSuffixBits, 0, inputCount);
    lastSuffixBits[0] = 0x80; //not worth using a random masking
    SnPMPC_AddByte(state, lastSuffixBits, rateInBytes-1, inputCount);
    #ifdef KeccakReference
    {
        unsigned char block[SnPMPC_width/8];
        size_t l;
        for( l = 0; l< inputCount; l++) {
            memset(block, 0, SnPMPC_width/8);
            block[rateInBytes-1] = lastSuffixBits[l];
            displayBytes(1, "Second bit of padding", block, rateInBytes);
        }
    }
    #endif

    SnPMPC_Permute(state);
    #ifdef KeccakReference
        displayText(1, "--- Switching to squeezing phase ---");
    #endif
    /* First, output whole blocks */
    while(outputByteLen > (size_t)rateInBytes) {
        SnPMPC_ExtractBytes(state, curOutputs, 0, inputCount, rateInBytes);
        SnPMPC_Permute(state);
        #ifdef KeccakReference
            displaySharedBytes(1, "Squeezed block", curOutputs, inputCount, rateInBytes);
        #endif

        //advance
        for( i = 0; i < inputCount; ++i) {
            curOutputs[i] += rateInBytes;
        }
        outputByteLen -= rateInBytes;
    }

    /* Finally, output what remains */
    partialBlock = (unsigned int)outputByteLen;
    SnPMPC_ExtractBytes(state, curOutputs, 0, inputCount, partialBlock);
    #ifdef KeccakReference
        displaySharedBytes(1, "Squeezed block (part)", curOutputs, inputCount, partialBlock);
    #endif

    return 0;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

int MaskedSpongeInitialize(MaskedSpongeInstance *instance, unsigned int rate, unsigned int capacity)
{
    if (rate+capacity != SnPMPC_width)
        return 1;
    if ((rate <= 0) || (rate > SnPMPC_width) || ((rate % 8) != 0))
        return 1;
    SnPMPC_StaticInitialize();
    SnPMPC_Initialize(instance->state);
    instance->rate = rate;
    instance->byteIOIndex = 0;
    instance->squeezing = 0;

    return 0;
}

/* ---------------------------------------------------------------- */

int MaskedSpongeAbsorb(MaskedSpongeInstance *instance, unsigned char **inputs, size_t inputCount, size_t inputByteLen)
{
    size_t i, j,k;
    unsigned int partialBlock;
    unsigned int rateInBytes = instance->rate/8;

    //ireversible state switch between absorbing and squeezing
    if (instance->squeezing)
        return 1; /* Too late for additional input */

    unsigned char *curInputs[SnPMPC_shareCount];
    for(i = 0; i < inputCount; ++i) {
        curInputs[i] = inputs[i];
    }

    i = 0; //i represents the number of processed bytes
    while(i < inputByteLen) {
        //rate is the chunk of message data we can pass to our f permutation at once
        //if this condition is true we can process full blocks in an aligned manner
        if ((instance->byteIOIndex == 0) && (inputByteLen >= (i + rateInBytes))) {
#ifdef SnP_FastLoop_Absorb
            /* processing full blocks first; SnP_width/200 is the number of lanes in the state */
            if ((rateInBytes % (SnP_width/200)) == 0) {
                /* fast lane: whole lane rate
                SnPMPC_FastLoop_Absorb should be defined as
                SnPMPC_FastLoop_Absorb(void * state, unsigned in laneCount, unsigned char **data, unsigned int dataCount, unsigned int dataByteLen, )
                function returns number of processed bytes
                */
                j = SnPMPC_FastLoop_Absorb(instance->state, rateInBytes/(SnP_width/200), curInputs, inputCount ,inputByteLen - i);
                i += j;
                //advance offset in input array
                for(k = 0; k < inputCount; ++k) {
                    curInputs[k] += j;
                }
            }
            else {
#endif
                //process chunks of size rateInBytes
                // as i is the number of processes bytes, j now is the number of
                //unprocessed bytes
                for(j=inputByteLen-i; j>=rateInBytes; j-=rateInBytes) {
                    #ifdef KeccakReference
                        displaySharedBytes(1, "Block to be absorbed", curInputs, inputCount, rateInBytes);
                    #endif
                    SnPMPC_AddBytes(instance->state, (const unsigned char**)curInputs, 0, inputCount, rateInBytes);
                    SnPMPC_Permute(instance->state);

                    for(k = 0; k < inputCount; ++k) {
                        curInputs[k] += rateInBytes;
                    }
                }
                //as j was initialized to the number of *unprocessed* bytes in the
                //preceding loop, this updates i to the new number of *processed* bytes
                i = inputByteLen - j;
#ifdef SnP_FastLoop_Absorb
            }
#endif
        }
        else { //we either have < rateInBytes bytes left or byteIOIndex was not zero "unaligned"
            /* normal lane: using the message queue */

            //(inputByteLen - i) <-> unprocessed bytes
            partialBlock = (unsigned int)(inputByteLen - i);

            //byteIOIndex keeps track of a previously started "unfinished" block
            if (partialBlock+instance->byteIOIndex > rateInBytes)
                partialBlock = rateInBytes - instance->byteIOIndex; //number of bytes from the new inputs to "finish" the previously unfinished block

            #ifdef KeccakReference
            displaySharedBytes(1, "Block to be absorbed (part)", curInputs, inputCount, partialBlock);
            #endif
            i += partialBlock;

            SnPMPC_AddBytes(instance->state, (const unsigned char**) curInputs, instance->byteIOIndex, inputCount, partialBlock);
            for(k = 0; k < inputCount; ++k) {
                curInputs[k] += partialBlock;
            }

            instance->byteIOIndex += partialBlock;
            if (instance->byteIOIndex == rateInBytes) {
                SnPMPC_Permute(instance->state);
                instance->byteIOIndex = 0;
            }
        }
    }

    return 0;
}

/* ---------------------------------------------------------------- */

int MaskedSpongeAbsorbLastFewBits(MaskedSpongeInstance *instance,  unsigned char * delimitedData, size_t delimitedDataCount)
{
    unsigned int rateInBytes = instance->rate/8;
    unsigned char * lastSuffixBits = NULL;

    if (delimitedData[0] == 0) {
        return 1;
    }

    if (instance->squeezing)
        return 1; /* Too late for additional input */
    #ifdef KeccakReference
        displaySharedByte(1, "Block to be absorbed (last few bits + first bit of padding)", delimitedData, delimitedDataCount);
    #endif

    /* Last few bits, whose delimiter coincides with first bit of padding */
    //just but them whereever we left of with byteIOIndex.
    SnPMPC_AddByte(instance->state, delimitedData, instance->byteIOIndex, delimitedDataCount);

    /* If the first bit of padding is at position rate-1, we need a whole new block for the second bit of padding */
    if ((delimitedData[0] >= 0x80) && (instance->byteIOIndex == (rateInBytes-1)))
        SnPMPC_Permute(instance->state);

    /* Second bit of padding */
    lastSuffixBits = alloca(delimitedDataCount);
    memset(lastSuffixBits, 0, delimitedDataCount);
    lastSuffixBits[0] = 0x80; //not worth using a random masking
    SnPMPC_AddByte(instance->state, lastSuffixBits, rateInBytes-1, delimitedDataCount);

    #ifdef KeccakReference
    {
        unsigned char block[SnPMPC_width/8];
        size_t l;
        for( l = 0; l< delimitedDataCount; l++) {
            memset(block, 0, SnPMPC_width/8);
            block[rateInBytes-1] = lastSuffixBits[l];
            displayBytes(1, "Second bit of padding", block, rateInBytes);
        }
    }
    #endif

    SnPMPC_Permute(instance->state);
    instance->byteIOIndex = 0;
    instance->squeezing = 1;
    #ifdef KeccakReference
        displayText(1, "--- Switching to squeezing phase ---");
    #endif

    return 0;
}

/* ---------------------------------------------------------------- */

int MaskedSpongeSqueeze(MaskedSpongeInstance *instance, unsigned char **outputs, size_t outputCount, size_t outputByteLen)
{
    size_t i, j,k;
    unsigned int partialBlock;
    unsigned int rateInBytes = instance->rate/8;
    unsigned char **curOutputs;

    if(outputByteLen == 0 ) {
        return 0;
    }

    curOutputs = alloca(sizeof(unsigned char*) * outputCount);
    for( k = 0; k < outputCount; ++k ) {
        curOutputs[k] = outputs[k];
    }

    if (!instance->squeezing) {
        unsigned char * delimitedData = alloca(outputCount);
        memset(delimitedData, 0, outputCount);
        delimitedData[0] = 0x01;
        MaskedSpongeAbsorbLastFewBits(instance, delimitedData, outputCount);
    }

    i = 0;
    while(i < outputByteLen) {
        if ((instance->byteIOIndex == rateInBytes) && (outputByteLen >= (i + rateInBytes))) {
            for(j=outputByteLen-i; j>=rateInBytes; j-=rateInBytes) {
                SnPMPC_Permute(instance->state);
                SnPMPC_ExtractBytes(instance->state, curOutputs, 0, outputCount, rateInBytes);

                #ifdef KeccakReference
                displaySharedBytes(1, "Squeezed block", curOutputs, outputCount, rateInBytes);
                #endif

                for( k = 0; k < outputCount; ++k ) {
                    curOutputs[k] += rateInBytes;
                }
            }
            i = outputByteLen - j;
        }
        else {
            /* normal lane: using the message queue */
            if (instance->byteIOIndex == rateInBytes) {
                SnPMPC_Permute(instance->state);
                instance->byteIOIndex = 0;
            }
            partialBlock = (unsigned int)(outputByteLen - i);
            if (partialBlock+instance->byteIOIndex > rateInBytes)
                partialBlock = rateInBytes-instance->byteIOIndex;
            i += partialBlock;

            SnPMPC_ExtractBytes(instance->state, curOutputs, instance->byteIOIndex, outputCount, partialBlock);
            #ifdef KeccakReference
            displaySharedBytes(1, "Squeezed block (part)", curOutputs, outputCount ,partialBlock);
            #endif


            for( k = 0; k < outputCount; ++k ) {
                curOutputs[k] += partialBlock;
            }
            instance->byteIOIndex += partialBlock;
        }
    }

    return 0;
}

/* ---------------------------------------------------------------- */

#undef MaskedSponge
#undef MaskedSpongeInstance
#undef MaskedSpongeSpongeInitialize
#undef MaskedSpongeSpongeAbsorb
#undef MaskedSpongeSpongeAbsorbLastFewBits
#undef MaskedSpongeSpongeSqueeze
#undef SnPMPC_stateSizeInBytes
#undef SnPMPC_stateAlignment
#undef SnPMPC_StaticInitialize
#undef SnPMPC_Initialize
#undef SnPMPC_AddByte
#undef SnPMPC_AddBytes
#undef SnPMPC_ExtractBytes
