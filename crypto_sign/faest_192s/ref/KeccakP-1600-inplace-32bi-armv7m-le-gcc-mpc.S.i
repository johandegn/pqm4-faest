@
@ The eXtended Keccak Code Package (XKCP)
@ https://github.com/XKCP/XKCP
@
@ The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.
@
@ Implementation by Ronny Van Keer, hereby denoted as "the implementer".
@
@ For more information, feedback or questions, please refer to the Keccak Team website:
@ https://keccak.team/
@
@ To the extent possible under law, the implementer has waived all copyright
@ and related or neighboring rights to the source code in this file.
@ http://creativecommons.org/publicdomain/zero/1.0/
@
@ ---
@
@ This file implements Keccak-p[1600] in a SnP-compatible way.
@ Please refer to SnP-documentation.h for more details.
@
@ This implementation comes with KeccakP-1600-SnP.h in the same folder.
@ Please refer to LowLevel.build for the exact list of other files it must be combined with.
@

@ WARNING: This implementation assumes a little endian CPU with ARMv7M architecture (e.g., Cortex-M3) and the GCC compiler.


    .thumb
    .syntax unified
.text

.equ RNG_BASE, 0x50060800
.equ RNG_SR_H, 0x5006
.equ RNG_SR_L, 0x0804
.equ RNG_DR_H, 0x5006
.equ RNG_DR_L, 0x0808
.equ RNG_SR_ERR_MASK, 0x6 @(SECS | CECS)
.equ RNG_SR_DRDY_MASK, 0x1

.align 8
@based on  rng_get_random of libopencm3
@ return rng number in \a
@ changes \a,\tmp
.macro  rndword   a, base, tmp
__KeccakP1600_RndWord_\@:
    @read memory mapped register
    ldr     \a, [\base, RNG_SR_L]
    @check error bits #1
    ands    \tmp, \a, RNG_SR_ERR_MASK
    @loop if we get error in status reg, on rng failure this will not terminate
    bne     __KeccakP1600_RndWord_\@
    @check error bits #2
    ldr     \a, [\base, RNG_SR_L]
    ands    \tmp, \a, RNG_SR_DRDY_MASK
    beq     __KeccakP1600_RndWord_\@
    @RNG_DR contains valid randomness, construct addr and load
    ldr     \a, [\tmp, RNG_DR_L]
.endm

@8018e18:       4a08            ldr     r2, [pc, #32]   ; (8018e3c <rng_get_random+0x24>)
@8018e1a:       f8d2 3804       ldr.w   r3, [r2, #2052] ; 0x804
@8018e1e:       f013 0f06       tst.w   r3, #6
@8018e22:       d109            bne.n   8018e38 <rng_get_random+0x20>
@8018e24:       f8d2 3804       ldr.w   r3, [r2, #2052] ; 0x804
@8018e28:       f013 0301       ands.w  r3, r3, #1
@8018e2c:       d004            beq.n   8018e38 <rng_get_random+0x20>
@8018e2e:       f8d2 3808       ldr.w   r3, [r2, #2056] ; 0x808
@8018e32:       6003            str     r3, [r0, #0]
@8018e34:       2001            movs    r0, #1
@8018e36:       4770            bx      lr
@8018e38:       2000            movs    r0, #0
@8018e3a:       4770            bx      lr
@8018e3c:       50060000        .word   0x50060000

@   --- offsets in state
.equ Aba0   , 0*4
.equ Aba1   , 1*4
.equ Abe0   , 2*4
.equ Abe1   , 3*4
.equ Abi0   , 4*4
.equ Abi1   , 5*4
.equ Abo0   , 6*4
.equ Abo1   , 7*4
.equ Abu0   , 8*4
.equ Abu1   , 9*4
.equ Aga0   , 10*4
.equ Aga1   , 11*4
.equ Age0   , 12*4
.equ Age1   , 13*4
.equ Agi0   , 14*4
.equ Agi1   , 15*4
.equ Ago0   , 16*4
.equ Ago1   , 17*4
.equ Agu0   , 18*4
.equ Agu1   , 19*4
.equ Aka0   , 20*4
.equ Aka1   , 21*4
.equ Ake0   , 22*4
.equ Ake1   , 23*4
.equ Aki0   , 24*4
.equ Aki1   , 25*4
.equ Ako0   , 26*4
.equ Ako1   , 27*4
.equ Aku0   , 28*4
.equ Aku1   , 29*4
.equ Ama0   , 30*4
.equ Ama1   , 31*4
.equ Ame0   , 32*4
.equ Ame1   , 33*4
.equ Ami0   , 34*4
.equ Ami1   , 35*4
.equ Amo0   , 36*4
.equ Amo1   , 37*4
.equ Amu0   , 38*4
.equ Amu1   , 39*4
.equ Asa0   , 40*4
.equ Asa1   , 41*4
.equ Ase0   , 42*4
.equ Ase1   , 43*4
.equ Asi0   , 44*4
.equ Asi1   , 45*4
.equ Aso0   , 46*4
.equ Aso1   , 47*4
.equ Asu0   , 48*4
.equ Asu1   , 49*4

.equ _Aba0   , 0*4+50*4
.equ _Aba1   , 1*4+50*4
.equ _Abe0   , 2*4+50*4
.equ _Abe1   , 3*4+50*4
.equ _Abi0   , 4*4+50*4
.equ _Abi1   , 5*4+50*4
.equ _Abo0   , 6*4+50*4
.equ _Abo1   , 7*4+50*4
.equ _Abu0   , 8*4+50*4
.equ _Abu1   , 9*4+50*4
.equ _Aga0   , 10*4+50*4
.equ _Aga1   , 11*4+50*4
.equ _Age0   , 12*4+50*4
.equ _Age1   , 13*4+50*4
.equ _Agi0   , 14*4+50*4
.equ _Agi1   , 15*4+50*4
.equ _Ago0   , 16*4+50*4
.equ _Ago1   , 17*4+50*4
.equ _Agu0   , 18*4+50*4
.equ _Agu1   , 19*4+50*4
.equ _Aka0   , 20*4+50*4
.equ _Aka1   , 21*4+50*4
.equ _Ake0   , 22*4+50*4
.equ _Ake1   , 23*4+50*4
.equ _Aki0   , 24*4+50*4
.equ _Aki1   , 25*4+50*4
.equ _Ako0   , 26*4+50*4
.equ _Ako1   , 27*4+50*4
.equ _Aku0   , 28*4+50*4
.equ _Aku1   , 29*4+50*4
.equ _Ama0   , 30*4+50*4
.equ _Ama1   , 31*4+50*4
.equ _Ame0   , 32*4+50*4
.equ _Ame1   , 33*4+50*4
.equ _Ami0   , 34*4+50*4
.equ _Ami1   , 35*4+50*4
.equ _Amo0   , 36*4+50*4
.equ _Amo1   , 37*4+50*4
.equ _Amu0   , 38*4+50*4
.equ _Amu1   , 39*4+50*4
.equ _Asa0   , 40*4+50*4
.equ _Asa1   , 41*4+50*4
.equ _Ase0   , 42*4+50*4
.equ _Ase1   , 43*4+50*4
.equ _Asi0   , 44*4+50*4
.equ _Asi1   , 45*4+50*4
.equ _Aso0   , 46*4+50*4
.equ _Aso1   , 47*4+50*4
.equ _Asu0   , 48*4+50*4
.equ _Asu1   , 49*4+50*4

@   --- offsets on stack
.equ mDa0   , 0*4
.equ mDa1   , 1*4
.equ mDe0   , 2*4
.equ mDe1   , 3*4
.equ mDi0   , 4*4
.equ mDi1   , 5*4
.equ mDo0   , 6*4
.equ mDo1   , 7*4
.equ mDu0   , 8*4
.equ mDu1   , 9*4
.equ _mDa0  , 10*4
.equ _mDa1  , 11*4
.equ _mDe0  , 12*4
.equ _mDe1  , 13*4
.equ _mDi0  , 14*4
.equ _mDi1  , 15*4
.equ _mDo0  , 16*4
.equ _mDo1  , 17*4
.equ _mDu0  , 18*4
.equ _mDu1  , 19*4
.equ mRC    , 20*4
#ifdef KECCAK_DOM_SECURE
.equ mZ0    , 21*4
.equ mZ1    , 22*4
.equ mZ2    , 23*4
.equ mZ3    , 24*4
.equ mZ4    , 25*4
.equ mSize  , 26*4
#elif defined(KECCAK_SNI_SECURE)
.equ mZ0    , 21*4
.equ mZ1    , 22*4
.equ mZ2    , 23*4
.equ mZ3    , 24*4
.equ mZ4    , 25*4
.equ mR0    , 26*4
.equ mR1    , 27*4
.equ mR2    , 28*4
.equ mR3    , 29*4
.equ mR4    , 30*4
.equ mSize  , 31*4
#else
.equ mSize  , 21*4
#endif

.macro    xandnot2     resofs, aa, bb, cc, dd, rzofs
    and         r2, \bb, \dd
#if defined(KECCAK_DOM_SECURE) || defined(KECCAK_SNI_SECURE)
    ldr         r1, [sp, #\rzofs]
    eor         r2, r2, r1
#endif
    bic         r1, \cc, \bb
    eor         r1, r1, \aa
    eor         r1, r1, r2
    str         r1, [r0, #\resofs]
.endm

.macro    _KeccakThetaRhoPiChiIota aA1, aDax, aA2, aDex, rot2, aA3, aDix, rot3, aA4, aDox, rot4, aA5, aDux, rot5, offset, last
    ldr     r3, [r0, #\aA1]
    ldr     r4, [r0, #\aA2]
    ldr     r5, [r0, #\aA3]
    ldr     r6, [r0, #\aA4]
    ldr     r7, [r0, #\aA5]
    ldr     r2, [sp, #\aDax]
    eor     r3, r3, r2
    ldr     r2, [sp, #\aDex]
    eor     r4, r4, r2
    ldr     r2, [sp, #\aDix]
    eor     r5, r5, r2
    ldr     r2, [sp, #\aDox]
    eor     r6, r6, r2
    ldr     r2, [sp, #\aDux]
    eor     r7, r7, r2
    rors    r4, #32-\rot2
    rors    r5, #32-\rot3
    rors    r6, #32-\rot4
    rors    r7, #32-\rot5

#if defined(KECCAK_DOM_SECURE) || defined(KECCAK_SNI_SECURE)
    mov     r11, #0
    movt    r11, #:upper16:RNG_BASE
#endif
#ifdef KECCAK_SNI_SECURE
    /* Gather randomness for refreshing permutation output. */
    rndword r1, r11, r12
    rndword r2, r11, r12
    rndword r8, r11, r12
    rndword r9, r11, r12
    rndword r10, r11, r12
    /* Pipeline writes. */
    str     r1, [sp, #mR0]
    str     r2, [sp, #mR1]
    str     r8, [sp, #mR2]
    str     r9, [sp, #mR3]
    str     r10, [sp, #mR4]
#endif
#if defined(KECCAK_DOM_SECURE) || defined(KECCAK_SNI_SECURE)
    /* Gather randomness for refreshing permutation output. */
    rndword r1, r11, r12
    rndword r2, r11, r12
    rndword r8, r11, r12
    rndword r9, r11, r12
    rndword r10, r11, r12
    /* Pipeline writes. */
    str     r1, [sp, #mZ0]
    str     r2, [sp, #mZ1]
    str     r8, [sp, #mZ2]
    str     r9, [sp, #mZ3]
    str     r10, [sp, #mZ4]
#endif

    ldr     r8, [r0, #\aA1+50*4]
    ldr     r9, [r0, #\aA2+50*4]
    ldr     r10, [r0, #\aA3+50*4]
    ldr     r11, [r0, #\aA4+50*4]
    ldr     r12, [r0, #\aA5+50*4]

#ifdef KECCAK_SNI_SECURE
    ldr     r1, [sp, #mR0]
    eor     r8, r8, r1
    eor     r3, r3, r1
    ldr     r1, [sp, #mR1]
    eor     r9, r9, r1
    eor     r4, r4, r1
    ldr     r1, [sp, #mR2]
    eor     r10, r10, r1
    eor     r5, r5, r1
    ldr     r1, [sp, #mR3]
    eor     r11, r11, r1
    eor     r6, r6, r1
    ldr     r1, [sp, #mR4]
    eor     r12, r12, r1
    eor     r7, r7, r1
#endif

    ldr     r2,  [sp, #\aDax+4*10]
    eor     r8, r8, r2
    ldr     r2,  [sp, #\aDex+4*10]
    eor     r9, r9, r2
    ldr     r2, [sp, #\aDix+4*10]
    eor     r10, r10, r2
    ldr     r2, [sp, #\aDox+4*10]
    eor     r11, r11, r2
    ldr     r2, [sp, #\aDux+4*10]
    eor     r12, r12, r2
    rors    r9,  #32-\rot2
    rors    r10, #32-\rot3
    rors    r11, #32-\rot4
    rors    r12, #32-\rot5

    xandnot2 \aA1+50*4, r8, r9, r10, r5, mZ0
    xandnot2 \aA2+50*4, r9, r10, r11, r6, mZ1
    xandnot2 \aA3+50*4, r10, r11, r12, r7, mZ2
    xandnot2 \aA4+50*4, r11, r12, r8, r3, mZ3
    xandnot2 \aA5+50*4, r12, r8, r9, r4, mZ4

    xandnot2 \aA2, r4, r5, r6, r11, mZ1
    xandnot2 \aA3, r5, r6, r7, r12, mZ2
    xandnot2 \aA4, r6, r7, r3, r8, mZ3
    xandnot2 \aA5, r7, r3, r4, r9, mZ4

    and     r10, r10, r4
    ldr     r1, [sp, #mRC]
    bics    r5, r5, r4
    ldr     r4, [r1, #\offset]
    eor     r3, r3, r5
    eor     r3, r3, r4
    eor     r3, r3, r10
#if defined(KECCAK_DOM_SECURE) || defined(KECCAK_SNI_SECURE)
    ldr     r10, [sp, #mZ0]
    eor     r3, r3, r10
#endif
    .if  \last == 1
    ldr     r4, [r1, #32]!
    str     r1, [sp, #mRC]
    cmp     r4, #0xFF
    .endif
    str     r3, [r0, #\aA1]
    .endm

.macro    _KeccakThetaRhoPiChi aB1, _aB1, aA1, aDax, rot1, aB2, _aB2, aA2, aDex, rot2, aB3, _aB3, aA3, aDix, rot3, aB4, _aB4, aA4, aDox, rot4, aB5, _aB5, aA5, aDux, rot5
    ldr     \aB1, [r0, #\aA1]
    ldr     \aB2, [r0, #\aA2]
    ldr     \aB3, [r0, #\aA3]
    ldr     \aB4, [r0, #\aA4]
    ldr     \aB5, [r0, #\aA5]
    ldr     r2,  [sp, #\aDax]
    eor     \aB1, \aB1, r2
    ldr     r2,  [sp, #\aDex]
    eor     \aB2, \aB2, r2
    ldr     r2, [sp, #\aDix]
    eor     \aB3, \aB3, r2
    ldr     r2, [sp, #\aDox]
    eor     \aB4, \aB4, r2
    ldr     r2, [sp, #\aDux]
    eor     \aB5, \aB5, r2
    rors    \aB1, #32-\rot1
    .if  \rot2 > 0
    rors    \aB2, #32-\rot2
    .endif
    rors    \aB3, #32-\rot3
    rors    \aB4, #32-\rot4
    rors    \aB5, #32-\rot5

#if defined(KECCAK_DOM_SECURE) || defined(KECCAK_SNI_SECURE)
    mov     r11, #0
    movt    r11, #:upper16:RNG_BASE
#endif
#ifdef KECCAK_SNI_SECURE
    /* Gather randomness for refreshing permutation output. */
    rndword r1, r11, r12
    rndword r2, r11, r12
    rndword r8, r11, r12
    rndword r9, r11, r12
    rndword r10, r11, r12
    /* Pipeline writes. */
    str     r1, [sp, #mR0]
    str     r2, [sp, #mR1]
    str     r8, [sp, #mR2]
    str     r9, [sp, #mR3]
    str     r10, [sp, #mR4]
#endif
#if defined(KECCAK_DOM_SECURE) || defined(KECCAK_SNI_SECURE)
    /* Gather randomness for refreshing permutation output. */
    rndword r1, r11, r12
    rndword r2, r11, r12
    rndword r8, r11, r12
    rndword r9, r11, r12
    rndword r10, r11, r12
    /* Pipeline writes. */
    str     r1, [sp, #mZ0]
    str     r2, [sp, #mZ1]
    str     r8, [sp, #mZ2]
    str     r9, [sp, #mZ3]
    str     r10, [sp, #mZ4]
#endif

    ldr     \_aB1, [r0, #\aA1+50*4]
    ldr     \_aB2, [r0, #\aA2+50*4]
    ldr     \_aB3, [r0, #\aA3+50*4]
    ldr     \_aB4, [r0, #\aA4+50*4]
    ldr     \_aB5, [r0, #\aA5+50*4]

#ifdef KECCAK_SNI_SECURE
    ldr     r1, [sp, #mR0]
    eor     r8, r8, r1
    eor     r3, r3, r1
    ldr     r1, [sp, #mR1]
    eor     r9, r9, r1
    eor     r4, r4, r1
    ldr     r1, [sp, #mR2]
    eor     r10, r10, r1
    eor     r5, r5, r1
    ldr     r1, [sp, #mR3]
    eor     r11, r11, r1
    eor     r6, r6, r1
    ldr     r1, [sp, #mR4]
    eor     r12, r12, r1
    eor     r7, r7, r1
#endif

    ldr     r2,  [sp, #\aDax+4*10]
    eor     \_aB1, \_aB1, r2
    ldr     r2,  [sp, #\aDex+4*10]
    eor     \_aB2, \_aB2, r2
    ldr     r2, [sp, #\aDix+4*10]
    eor     \_aB3, \_aB3, r2
    ldr     r2, [sp, #\aDox+4*10]
    eor     \_aB4, \_aB4, r2
    ldr     r2, [sp, #\aDux+4*10]
    eor     \_aB5, \_aB5, r2
    rors    \_aB1, #32-\rot1
    .if  \rot2 > 0
    rors    \_aB2, #32-\rot2
    .endif
    rors    \_aB3, #32-\rot3
    rors    \_aB4, #32-\rot4
    rors    \_aB5, #32-\rot5
    xandnot2 \aA1, r3, r4, r5, r10, mZ0
    xandnot2 \aA2, r4, r5, r6, r11, mZ1
    xandnot2 \aA3, r5, r6, r7, r12, mZ2
    xandnot2 \aA4, r6, r7, r3, r8, mZ3
    xandnot2 \aA5, r7, r3, r4, r9, mZ4
    xandnot2 \aA1+50*4, r8, r9, r10, r5, mZ0
    xandnot2 \aA2+50*4, r9, r10, r11, r6, mZ1
    xandnot2 \aA3+50*4, r10, r11, r12, r7, mZ2
    xandnot2 \aA4+50*4, r11, r12, r8, r3, mZ3
    xandnot2 \aA5+50*4, r12, r8, r9, r4, mZ4
.endm

.macro    _KeccakRound0

    xor5        r3, Abu0, Agu0, Aku0, Amu0, Asu0
    xor5        r7, Abe1, Age1, Ake1, Ame1, Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #mDa0]
    xor5        r6,  Abu1, Agu1, Aku1, Amu1, Asu1
    xor5        lr, Abe0, Age0, Ake0, Ame0, Ase0
    eor         r8, r6, lr
    str         r8, [sp, #mDa1]

    xor5        r5,  Abi0, Agi0, Aki0, Ami0, Asi0
    xorrol      r9, r5, r6
    str         r9, [sp, #mDo0]
    xor5        r4,  Abi1, Agi1, Aki1, Ami1, Asi1
    eor         r3, r3, r4
    str         r3, [sp, #mDo1]

    xor5        r3,  Aba0, Aga0, Aka0, Ama0, Asa0
    xorrol      r10, r3, r4
    str         r10, [sp, #mDe0]
    xor5        r6,  Aba1, Aga1, Aka1, Ama1, Asa1
    eor         r11, r6, r5
    str         r11, [sp, #mDe1]

    xor5        r4,  Abo1, Ago1, Ako1, Amo1, Aso1
    xorrol      r5, lr, r4
    str         r5, [sp, #mDi0]
    xor5        r5,  Abo0, Ago0, Ako0, Amo0, Aso0
    eor         r7, r7, r5
    str         r7, [sp, #mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #mDu0]
    str         r4, [sp, #mDu1]

    xor5        r3,  _Abu0, _Agu0, _Aku0, _Amu0, _Asu0
    xor5        r7, _Abe1, _Age1, _Ake1, _Ame1, _Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #_mDa0]
    xor5        r6,  _Abu1, _Agu1, _Aku1, _Amu1, _Asu1
    xor5        lr, _Abe0, _Age0, _Ake0, _Ame0, _Ase0
    eor         r9, r6, lr
    str         r9, [sp, #_mDa1]

    xor5        r5,  _Abi0, _Agi0, _Aki0, _Ami0, _Asi0
    xorrol      r8, r5, r6
    str         r8, [sp, #_mDo0]
    xor5        r4,  _Abi1, _Agi1, _Aki1, _Ami1, _Asi1
    eor         r3, r3, r4
    str         r3, [sp, #_mDo1]

    xor5        r3,  _Aba0, _Aga0, _Aka0, _Ama0, _Asa0
    xorrol      r11, r3, r4
    str         r11, [sp, #_mDe0]
    xor5        r6,  _Aba1, _Aga1, _Aka1, _Ama1, _Asa1
    eor         r10, r6, r5
    str         r10, [sp, #_mDe1]

    xor5        r4,  _Abo1, _Ago1, _Ako1, _Amo1, _Aso1
    xorrol      r5, lr, r4
    str         r5, [sp, #_mDi0]
    xor5        r5,  _Abo0, _Ago0, _Ako0, _Amo0, _Aso0
    eor         r7, r7, r5
    str         r7, [sp, #_mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #_mDu0]
    str         r4, [sp, #_mDu1]

    _KeccakThetaRhoPiChi r5, r10, Aka1, mDa1,  2, r6, r11, Ame1, mDe1, 23, r7, r12, Asi1, mDi1, 31, r3, r8,  Abo0, mDo0, 14, r4, r9,  Agu0, mDu0, 10
    _KeccakThetaRhoPiChi r7, r12, Asa1, mDa1,  9, r3, r8,  Abe0, mDe0,  0, r4, r9,  Agi1, mDi1,  3, r5, r10, Ako0, mDo0, 12, r6, r11, Amu1, mDu1,  4
    _KeccakThetaRhoPiChi r4, r9,  Aga0, mDa0, 18, r5, r10, Ake0, mDe0,  5, r6, r11, Ami1, mDi1,  8, r7, r12, Aso0, mDo0, 28, r3, r8,  Abu1, mDu1, 14
    _KeccakThetaRhoPiChi r6, r11, Ama0, mDa0, 20, r7, r12, Ase1, mDe1,  1, r3, r8,  Abi1, mDi1, 31, r4, r9,  Ago0, mDo0, 27, r5, r10, Aku0, mDu0, 19
    _KeccakThetaRhoPiChiIota Aba0, mDa0,         Age0, mDe0, 22,      Aki1, mDi1, 22,    Amo1, mDo1, 11,     Asu0, mDu0,  7, 0, 0

    _KeccakThetaRhoPiChi r5, r10, Aka0, mDa0,  1, r6, r11, Ame0, mDe0, 22, r7, r12, Asi0, mDi0, 30, r3, r8,  Abo1, mDo1, 14, r4, r9,  Agu1, mDu1, 10
    _KeccakThetaRhoPiChi r7, r12, Asa0, mDa0,  9, r3, r8,  Abe1, mDe1,  1, r4, r9,  Agi0, mDi0,  3, r5, r10, Ako1, mDo1, 13, r6, r11, Amu0, mDu0,  4
    _KeccakThetaRhoPiChi r4, r9,  Aga1, mDa1, 18, r5, r10, Ake1, mDe1,  5, r6, r11, Ami0, mDi0,  7, r7, r12, Aso1, mDo1, 28, r3, r8,  Abu0, mDu0, 13
    _KeccakThetaRhoPiChi r6, r11, Ama1, mDa1, 21, r7, r12, Ase0, mDe0,  1, r3, r8,  Abi0, mDi0, 31, r4, r9,  Ago1, mDo1, 28, r5, r10, Aku1, mDu1, 20
    _KeccakThetaRhoPiChiIota Aba1, mDa1,         Age1, mDe1, 22,     Aki0, mDi0, 21,     Amo0, mDo0, 10,     Asu1, mDu1,  7, 4, 0
    .endm

.macro    _KeccakRound1

    xor5        r3, Asu0, Agu0, Amu0, Abu1, Aku1
    xor5        r7, Age1, Ame0, Abe0, Ake1, Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #mDa0]
    xor5        r6,  Asu1, Agu1, Amu1, Abu0, Aku0
    xor5        lr, Age0, Ame1, Abe1, Ake0, Ase0
    eor         r8, r6, lr
    str         r8, [sp, #mDa1]

    xor5        r5,  Aki1, Asi1, Agi0, Ami1, Abi0
    xorrol      r9, r5, r6
    str         r9, [sp, #mDo0]
    xor5        r4,  Aki0, Asi0, Agi1, Ami0, Abi1
    eor         r3, r3, r4
    str         r3, [sp, #mDo1]

    xor5        r3,  Aba0, Aka1, Asa0, Aga0, Ama1
    xorrol      r10, r3, r4
    str         r10, [sp, #mDe0]
    xor5        r6,  Aba1, Aka0, Asa1, Aga1, Ama0
    eor         r11, r6, r5
    str         r11, [sp, #mDe1]

    xor5        r4,  Amo0, Abo1, Ako0, Aso1, Ago0
    xorrol      r5, lr, r4
    str         r5, [sp, #mDi0]
    xor5        r5,  Amo1, Abo0, Ako1, Aso0, Ago1
    eor         r7, r7, r5
    str         r7, [sp, #mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #mDu0]
    str         r4, [sp, #mDu1]

    xor5        r3,  _Asu0, _Agu0, _Amu0, _Abu1, _Aku1
    xor5        r7, _Age1, _Ame0, _Abe0, _Ake1, _Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #_mDa0]
    xor5        r6,  _Asu1, _Agu1, _Amu1, _Abu0, _Aku0
    xor5        lr, _Age0, _Ame1, _Abe1, _Ake0, _Ase0
    eor         r9, r6, lr
    str         r9, [sp, #_mDa1]

    xor5        r5,  _Aki1, _Asi1, _Agi0, _Ami1, _Abi0
    xorrol      r8, r5, r6
    str         r8, [sp, #_mDo0]
    xor5        r4,  _Aki0, _Asi0, _Agi1, _Ami0, _Abi1
    eor         r3, r3, r4
    str         r3, [sp, #_mDo1]

    xor5        r3,  _Aba0, _Aka1, _Asa0, _Aga0, _Ama1
    xorrol      r11, r3, r4
    str         r11, [sp, #_mDe0]
    xor5        r6,  _Aba1, _Aka0, _Asa1, _Aga1, _Ama0
    eor         r10, r6, r5
    str         r10, [sp, #_mDe1]

    xor5        r4,  _Amo0, _Abo1, _Ako0, _Aso1, _Ago0
    xorrol      r5, lr, r4
    str         r5, [sp, #_mDi0]
    xor5        r5,  _Amo1, _Abo0, _Ako1, _Aso0, _Ago1
    eor         r7, r7, r5
    str         r7, [sp, #_mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #_mDu0]
    str         r4, [sp, #_mDu1]

    _KeccakThetaRhoPiChi r5, r10, Asa1, mDa1,  2, r6, r11, Ake1, mDe1, 23, r7, r12, Abi1, mDi1, 31, r3, r8,  Amo1, mDo0, 14, r4, r9,  Agu0, mDu0, 10
    _KeccakThetaRhoPiChi r7, r12, Ama0, mDa1,  9, r3, r8,  Age0, mDe0,  0, r4, r9,  Asi0, mDi1,  3, r5, r10, Ako1, mDo0, 12, r6, r11, Abu0, mDu1,  4
    _KeccakThetaRhoPiChi r4, r9,  Aka1, mDa0, 18, r5, r10, Abe1, mDe0,  5, r6, r11, Ami0, mDi1,  8, r7, r12, Ago1, mDo0, 28, r3, r8,  Asu1, mDu1, 14
    _KeccakThetaRhoPiChi r6, r11, Aga0, mDa0, 20, r7, r12, Ase1, mDe1,  1, r3, r8,  Aki0, mDi1, 31, r4, r9,  Abo0, mDo0, 27, r5, r10, Amu0, mDu0, 19
    _KeccakThetaRhoPiChiIota Aba0, mDa0,         Ame1, mDe0, 22,     Agi1, mDi1, 22,     Aso1, mDo1, 11,     Aku1, mDu0,  7, 8, 0

    _KeccakThetaRhoPiChi r5, r10, Asa0, mDa0,  1, r6, r11, Ake0, mDe0, 22, r7, r12, Abi0, mDi0, 30, r3, r8,  Amo0, mDo1, 14, r4, r9,  Agu1, mDu1, 10
    _KeccakThetaRhoPiChi r7, r12, Ama1, mDa0,  9, r3, r8,  Age1, mDe1,  1, r4, r9,  Asi1, mDi0,  3, r5, r10, Ako0, mDo1, 13, r6, r11, Abu1, mDu0,  4
    _KeccakThetaRhoPiChi r4, r9,  Aka0, mDa1, 18, r5, r10, Abe0, mDe1,  5, r6, r11, Ami1, mDi0,  7, r7, r12, Ago0, mDo1, 28, r3, r8,  Asu0, mDu0, 13
    _KeccakThetaRhoPiChi r6, r11, Aga1, mDa1, 21, r7, r12, Ase0, mDe0,  1, r3, r8,  Aki1, mDi0, 31, r4, r9,  Abo1, mDo1, 28, r5, r10, Amu1, mDu1, 20
    _KeccakThetaRhoPiChiIota Aba1, mDa1,         Ame0, mDe1, 22,     Agi0, mDi0, 21,     Aso0, mDo0, 10,     Aku0, mDu1,  7, 12, 0
    .endm

.macro    _KeccakRound2

    xor5        r3, Aku1, Agu0, Abu1, Asu1, Amu1
    xor5        r7, Ame0, Ake0, Age0, Abe0, Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #mDa0]
    xor5        r6,  Aku0, Agu1, Abu0, Asu0, Amu0
    xor5        lr, Ame1, Ake1, Age1, Abe1, Ase0
    eor         r8, r6, lr
    str         r8, [sp, #mDa1]

    xor5        r5,  Agi1, Abi1, Asi1, Ami0, Aki1
    xorrol      r9, r5, r6
    str         r9, [sp, #mDo0]
    xor5        r4,  Agi0, Abi0, Asi0, Ami1, Aki0
    eor         r3, r3, r4
    str         r3, [sp, #mDo1]

    xor5        r3,  Aba0, Asa1, Ama1, Aka1, Aga1
    xorrol      r10, r3, r4
    str         r10, [sp, #mDe0]
    xor5        r6,  Aba1, Asa0, Ama0, Aka0, Aga0
    eor         r11, r6, r5
    str         r11, [sp, #mDe1]

    xor5        r4,  Aso0, Amo0, Ako1, Ago0, Abo0
    xorrol      r5, lr, r4
    str         r5, [sp, #mDi0]
    xor5        r5,  Aso1, Amo1, Ako0, Ago1, Abo1
    eor         r7, r7, r5
    str         r7, [sp, #mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #mDu0]
    str         r4, [sp, #mDu1]

    xor5        r3, _Aku1, _Agu0, _Abu1, _Asu1, _Amu1
    xor5        r7, _Ame0, _Ake0, _Age0, _Abe0, _Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #_mDa0]
    xor5        r6,  _Aku0, _Agu1, _Abu0, _Asu0, _Amu0
    xor5        lr, _Ame1, _Ake1, _Age1, _Abe1, _Ase0
    eor         r9, r6, lr
    str         r9, [sp, #_mDa1]

    xor5        r5,  _Agi1, _Abi1, _Asi1, _Ami0, _Aki1
    xorrol      r8, r5, r6
    str         r8, [sp, #_mDo0]
    xor5        r4,  _Agi0, _Abi0, _Asi0, _Ami1, _Aki0
    eor         r3, r3, r4
    str         r3, [sp, #_mDo1]

    xor5        r3,  _Aba0, _Asa1, _Ama1, _Aka1, _Aga1
    xorrol      r11, r3, r4
    str         r11, [sp, #_mDe0]
    xor5        r6,  _Aba1, _Asa0, _Ama0, _Aka0, _Aga0
    eor         r10, r6, r5
    str         r10, [sp, #_mDe1]

    xor5        r4,  _Aso0, _Amo0, _Ako1, _Ago0, _Abo0
    xorrol      r5, lr, r4
    str         r5, [sp, #_mDi0]
    xor5        r5,  _Aso1, _Amo1, _Ako0, _Ago1, _Abo1
    eor         r7, r7, r5
    str         r7, [sp, #_mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #_mDu0]
    str         r4, [sp, #_mDu1]

    _KeccakThetaRhoPiChi r5, r10, Ama0, mDa1,  2, r6, r11, Abe0, mDe1, 23, r7, r12, Aki0, mDi1, 31, r3, r8,  Aso1, mDo0, 14, r4, r9,  Agu0, mDu0, 10
    _KeccakThetaRhoPiChi r7, r12, Aga0, mDa1,  9, r3, r8,  Ame1, mDe0,  0, r4, r9,  Abi0, mDi1,  3, r5, r10, Ako0, mDo0, 12, r6, r11, Asu0, mDu1,  4
    _KeccakThetaRhoPiChi r4, r9,  Asa1, mDa0, 18, r5, r10, Age1, mDe0,  5, r6, r11, Ami1, mDi1,  8, r7, r12, Abo1, mDo0, 28, r3, r8,  Aku0, mDu1, 14
    _KeccakThetaRhoPiChi r6, r11, Aka1, mDa0, 20, r7, r12, Ase1, mDe1,  1, r3, r8,  Agi0, mDi1, 31, r4, r9,  Amo1, mDo0, 27, r5, r10, Abu1, mDu0, 19
    _KeccakThetaRhoPiChiIota Aba0, mDa0,         Ake1, mDe0, 22,     Asi0, mDi1, 22,     Ago0, mDo1, 11,     Amu1, mDu0,  7, 16, 0

    _KeccakThetaRhoPiChi r5, r10, Ama1, mDa0,  1, r6, r11, Abe1, mDe0, 22, r7, r12, Aki1, mDi0, 30, r3, r8,  Aso0, mDo1, 14, r4, r9,  Agu1, mDu1, 10
    _KeccakThetaRhoPiChi r7, r12, Aga1, mDa0,  9, r3, r8,  Ame0, mDe1,  1, r4, r9,  Abi1, mDi0,  3, r5, r10, Ako1, mDo1, 13, r6, r11, Asu1, mDu0,  4
    _KeccakThetaRhoPiChi r4, r9,  Asa0, mDa1, 18, r5, r10, Age0, mDe1,  5, r6, r11, Ami0, mDi0,  7, r7, r12, Abo0, mDo1, 28, r3, r8,  Aku1, mDu0, 13
    _KeccakThetaRhoPiChi r6, r11, Aka0, mDa1, 21, r7, r12, Ase0, mDe0,  1, r3, r8,  Agi1, mDi0, 31, r4, r9,  Amo0, mDo1, 28, r5, r10, Abu0, mDu1, 20
    _KeccakThetaRhoPiChiIota Aba1, mDa1,         Ake0, mDe1, 22,     Asi1, mDi0, 21,     Ago1, mDo0, 10,     Amu0, mDu1,  7, 20, 0
    .endm

.macro    _KeccakRound3

    xor5        r3, Amu1, Agu0, Asu1, Aku0, Abu0
    xor5        r7, Ake0, Abe1, Ame1, Age0, Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #mDa0]
    xor5        r6,  Amu0, Agu1, Asu0, Aku1, Abu1
    xor5        lr, Ake1, Abe0, Ame0, Age1, Ase0
    eor         r8, r6, lr
    str         r8, [sp, #mDa1]

    xor5        r5,  Asi0, Aki0, Abi1, Ami1, Agi1
    xorrol      r9, r5, r6
    str         r9, [sp, #mDo0]
    xor5        r4,  Asi1, Aki1, Abi0, Ami0, Agi0
    eor         r3, r3, r4
    str         r3, [sp, #mDo1]

    xor5        r3,  Aba0, Ama0, Aga1, Asa1, Aka0
    xorrol      r10, r3, r4
    str         r10, [sp, #mDe0]
    xor5        r6,  Aba1, Ama1, Aga0, Asa0, Aka1
    eor         r11, r6, r5
    str         r11, [sp, #mDe1]

    xor5        r4,  Ago1, Aso0, Ako0, Abo0, Amo1
    xorrol      r5, lr, r4
    str         r5, [sp, #mDi0]
    xor5        r5,  Ago0, Aso1, Ako1, Abo1, Amo0
    eor         r7, r7, r5
    str         r7, [sp, #mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #mDu0]
    str         r4, [sp, #mDu1]

    xor5        r3,  _Amu1, _Agu0, _Asu1, _Aku0, _Abu0
    xor5        r7, _Ake0, _Abe1, _Ame1, _Age0, _Ase1
    xorrol      r6, r3, r7
    str         r6, [sp, #_mDa0]
    xor5        r6,  _Amu0, _Agu1, _Asu0, _Aku1, _Abu1
    xor5        lr, _Ake1, _Abe0, _Ame0, _Age1, _Ase0
    eor         r9, r6, lr
    str         r9, [sp, #_mDa1]

    xor5        r5,  _Asi0, _Aki0, _Abi1, _Ami1, _Agi1
    xorrol      r8, r5, r6
    str         r8, [sp, #_mDo0]
    xor5        r4,  _Asi1, _Aki1, _Abi0, _Ami0, _Agi0
    eor         r3, r3, r4
    str         r3, [sp, #_mDo1]

    xor5        r3,  _Aba0, _Ama0, _Aga1, _Asa1, _Aka0
    xorrol      r11, r3, r4
    str         r11, [sp, #_mDe0]
    xor5        r6,  _Aba1, _Ama1, _Aga0, _Asa0, _Aka1
    eor         r10, r6, r5
    str         r10, [sp, #_mDe1]

    xor5        r4,  _Ago1, _Aso0, _Ako0, _Abo0, _Amo1
    xorrol      r5, lr, r4
    str         r5, [sp, #_mDi0]
    xor5        r5,  _Ago0, _Aso1, _Ako1, _Abo1, _Amo0
    eor         r7, r7, r5
    str         r7, [sp, #_mDi1]

    xorrol      r5, r5, r6
    eor         r4, r4, r3
    str         r5, [sp, #_mDu0]
    str         r4, [sp, #_mDu1]

    _KeccakThetaRhoPiChi r5, r10, Aga0, mDa1,  2, r6, r11, Age0, mDe1, 23, r7, r12, Agi0, mDi1, 31, r3, r8,  Ago0, mDo0, 14, r4, r9,  Agu0, mDu0, 10
    _KeccakThetaRhoPiChi r7, r12, Aka1, mDa1,  9, r3, r8,  Ake1, mDe0,  0, r4, r9,  Aki1, mDi1,  3, r5, r10, Ako1, mDo0, 12, r6, r11, Aku1, mDu1,  4
    _KeccakThetaRhoPiChi r4, r9,  Ama0, mDa0, 18, r5, r10, Ame0, mDe0,  5, r6, r11, Ami0, mDi1,  8, r7, r12, Amo0, mDo0, 28, r3, r8,  Amu0, mDu1, 14
    _KeccakThetaRhoPiChi r6, r11, Asa1, mDa0, 20, r7, r12, Ase1, mDe1,  1, r3, r8,  Asi1, mDi1, 31, r4, r9,  Aso1, mDo0, 27, r5, r10, Asu1, mDu0, 19
    _KeccakThetaRhoPiChiIota Aba0, mDa0,         Abe0, mDe0, 22,     Abi0, mDi1, 22,     Abo0, mDo1, 11,     Abu0, mDu0,  7, 24, 0

    _KeccakThetaRhoPiChi r5, r10, Aga1, mDa0,  1, r6, r11, Age1, mDe0, 22, r7, r12, Agi1, mDi0, 30, r3, r8,  Ago1, mDo1, 14, r4, r9,  Agu1, mDu1, 10
    _KeccakThetaRhoPiChi r7, r12, Aka0, mDa0,  9, r3, r8,  Ake0, mDe1,  1, r4, r9,  Aki0, mDi0,  3, r5, r10, Ako0, mDo1, 13, r6, r11, Aku0, mDu0,  4
    _KeccakThetaRhoPiChi r4, r9,  Ama1, mDa1, 18, r5, r10, Ame1, mDe1,  5, r6, r11, Ami1, mDi0,  7, r7, r12, Amo1, mDo1, 28, r3, r8,  Amu1, mDu0, 13
    _KeccakThetaRhoPiChi r6, r11, Asa0, mDa1, 21, r7, r12, Ase0, mDe0,  1, r3, r8,  Asi0, mDi0, 31, r4, r9,  Aso0, mDo1, 28, r5, r10, Asu0, mDu1, 20
    _KeccakThetaRhoPiChiIota Aba1, mDa1,         Abe1, mDe1, 22,     Abi1, mDi0, 21,     Abo1, mDo0, 10,     Abu1, mDu1,  7, 28, 1
    .endm


@ ----------------------------------------------------------------------------
@
@  void KeccakP1600MPC_Permute_Nrounds(void *state, unsigned int nrounds)
@
.align 8
.global   KeccakP1600MPC_Permute_Nrounds
.type    KeccakP1600MPC_Permute_Nrounds, %function;
KeccakP1600MPC_Permute_Nrounds:
    lsls    r2, r1, #3
    adr     r1, KeccakP1600MPC_Permute_RoundConstants0Mod4
    subs    r1, r1, r2
    b       KeccakP1600MPC_Permute

@ ----------------------------------------------------------------------------
@
@  void KeccakP1600MPC_Permute_12rounds( void *state )
@
.align 8
.global   KeccakP1600MPC_Permute_12rounds
.type    KeccakP1600MPC_Permute_12rounds, %function;
KeccakP1600MPC_Permute_12rounds:
    adr     r1, KeccakP1600MPC_Permute_HalfRoundConstants24
    b       KeccakP1600MPC_Permute

@ ----------------------------------------------------------------------------
@
@  void KeccakP1600MPC_Permute_24rounds( void *state )
@
.align 8
.global   KeccakP1600MPC_Permute_24rounds
.type    KeccakP1600MPC_Permute_24rounds, %function;
KeccakP1600MPC_Permute_24rounds:
    adr     r1, KeccakP1600MPC_Permute_RoundConstants24
    b       KeccakP1600MPC_Permute

.align 8
KeccakP1600MPC_Permute_HalfRoundConstants24:
    @       0           1
        .long      0x00000001, 0x00000000
        .long      0x00000000, 0x00000089
        .long      0x00000000, 0x8000008b
        .long      0x00000000, 0x80008080
        .long      0x00000001, 0x0000008b
        .long      0x00000001, 0x00008000
        .long      0x00000001, 0x80008088
        .long      0x00000001, 0x80000082
        .long      0x00000000, 0x0000000b
        .long      0x00000000, 0x0000000a
        .long      0x00000001, 0x00008082
        .long      0x00000000, 0x00008003
        .long      0x000000FF  @terminator

.align 8
KeccakP1600MPC_Permute_RoundConstants24:
    @       0           1
        .long      0x00000001, 0x00000000
        .long      0x00000000, 0x00000089
        .long      0x00000000, 0x8000008b
        .long      0x00000000, 0x80008080
        .long      0x00000001, 0x0000008b
        .long      0x00000001, 0x00008000
        .long      0x00000001, 0x80008088
        .long      0x00000001, 0x80000082
        .long      0x00000000, 0x0000000b
        .long      0x00000000, 0x0000000a
        .long      0x00000001, 0x00008082
        .long      0x00000000, 0x00008003
KeccakP1600MPC_Permute_RoundConstants12:
        .long      0x00000001, 0x0000808b
        .long      0x00000001, 0x8000000b
        .long      0x00000001, 0x8000008a
        .long      0x00000001, 0x80000081
        .long      0x00000000, 0x80000081
        .long      0x00000000, 0x80000008
        .long      0x00000000, 0x00000083
        .long      0x00000000, 0x80008003
KeccakP1600MPC_Permute_RoundConstants0:
        .long      0x00000001, 0x80008088
        .long      0x00000000, 0x80000088
        .long      0x00000001, 0x00008000
        .long      0x00000000, 0x80008082
KeccakP1600MPC_Permute_RoundConstants0Mod4:
        .long      0x000000FF  @terminator

@----------------------------------------------------------------------------
@
@ void KeccakP1600MPC_Permute( void *state, void * rc )
@
.align 8
KeccakP1600MPC_Permute:
    push    { r4 - r12, lr }
    sub     sp, #mSize
    str     r1, [sp, #mRC]
KeccakP1600MPC_Permute_RoundLoop:
    _KeccakRound0
    _KeccakRound1
    _KeccakRound2
    _KeccakRound3
#ifndef KECCAK_REDUCED_ROUNDS
    bne     KeccakP1600MPC_Permute_RoundLoop
#endif
    add     sp, #mSize
    pop     { r4 - r12, pc }
