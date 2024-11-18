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

#ifndef _KeccakMPCHashInterface_h_
#define _KeccakMPCHashInterface_h_

#include "config.h"
#ifdef XKCP_has_MaskedKeccakP1600

#include <stdint.h>
#include <string.h>
#include "KeccakMPCSponge.h"

#ifndef _Keccak_BitTypes_
#define _Keccak_BitTypes_
typedef uint8_t BitSequence;

typedef size_t BitLength;
#endif

typedef enum { KECCAKMPC_SUCCESS = 0, KECCAKMPC_FAIL = 1, KECCAKMPC_BAD_HASHLEN = 2 } MPCHashReturn;

typedef struct {
    KeccakMaskedWidth1600_MaskedSpongeInstance sponge;
    unsigned int fixedOutputLength;
    unsigned char *delimitedSuffixShares;
    unsigned int shareCount;
} KeccakMPC_HashInstance;

/**
 * Function to initialize the Keccak[r, c] sponge function instance used in sequential hashing mode.
 * @param  hashInstance    Pointer to the hash instance to be initialized.
 * @param  rate        The value of the rate r.
 * @param  capacity    The value of the capacity c.
 * @param  hashbitlen  The desired number of output bits,
 *                     or 0 for an arbitrarily-long output.
 * @param  delimitedSuffix Bits that will be automatically appended to the end
 *                         of the input message, as in domain separation.
 *                         This is a byte containing from 0 to 7 bits
 *                         formatted like the @a delimitedData parameter of
 *                         the Keccak_SpongeAbsorbLastFewBits() function.
 * @pre    One must have r+c=1600 and the rate a multiple of 8 bits in this implementation.
 * @return KECCAK_SUCCESS if successful, KECCAK_FAIL otherwise.
 */
MPCHashReturn KeccakMPC_HashInitialize(KeccakMPC_HashInstance *hashInstance, unsigned int rate, unsigned int capacity, unsigned int hashbitlen, unsigned char delimitedSuffix, unsigned int shareCount);

/*
 * Macro to initialize a SHAKE256 instance as specified in the FIPS 202 standard.
*/
#define KeccakMPC_HashInitialize_2SHARE_SHAKE256(hashInstance)        KeccakMPC_HashInitialize(hashInstance, 1088,  512,   0, 0x1F, 2)
#define KeccakMPC_HashInitialize_2SHARE_SHAKE128(hashInstance)        KeccakMPC_HashInitialize(hashInstance, 1344,  256,   0, 0x1F, 2)

/**
  * Function to give input data to be absorbed.
  * @param  hashInstance    Pointer to the hash instance initialized by Keccak_HashInitialize().
  * @param  data        Pointer to the already masked input data shares
  *                     When @a databitLen is not a multiple of 8, the last bits of data must be
  *                     in the least significant bits of the last byte (little-endian convention).
  *                     In this case, the (8 - @a databitLen mod 8) most significant bits
  *                     of the last byte are ignored.
  * @param  databitLen  The number of input bits provided in the input data.
  * @pre    In the previous call to Keccak_HashUpdate(), databitlen was a multiple of 8.
  * @return KECCAK_SUCCESS if successful, KECCAK_FAIL otherwise.
  */
MPCHashReturn KeccakMPC_HashUpdate(KeccakMPC_HashInstance *hashInstance, BitSequence **data, BitLength databitlen);

/**
  * Function to call after all input blocks have been input and to get
  * output bits if the length was specified when calling Keccak_HashInitialize().
  * @param  hashInstance    Pointer to the hash instance initialized by Keccak_HashInitialize().
  * If @a hashbitlen was not 0 in the call to Keccak_HashInitialize(), the number of
  *     output bits is equal to @a hashbitlen.
  * If @a hashbitlen was 0 in the call to Keccak_HashInitialize(), the output bits
  *     must be extracted using the Keccak_HashSqueeze() function.
  * @param  hashvals     Pointer to the buffer where to store the output shares of the data.
  * @return KECCAK_SUCCESS if successful, KECCAK_FAIL otherwise.
  */
MPCHashReturn KeccakMPC_HashFinal(KeccakMPC_HashInstance *hashInstance, BitSequence **hashvals);

/**
 * Function to squeeze output data.
 * @param  hashInstance    Pointer to the hash instance initialized by Keccak_HashInitialize().
 * @param  data        Pointer to the buffer where to store the output data.
 * @param  databitlen  The number of output bits desired (must be a multiple of 8).
 * @pre    Keccak_HashFinal() must have been already called.
 * @pre    @a databitlen is a multiple of 8.
 * @return KECCAK_SUCCESS if successful, KECCAK_FAIL otherwise.
 */
MPCHashReturn KeccakMPC_HashSqueeze(KeccakMPC_HashInstance *hashInstance, BitSequence **data, BitLength databitlen);

#else
#error This requires a Masked mpc implementation of Keccak-p[1600]
#endif

#endif
