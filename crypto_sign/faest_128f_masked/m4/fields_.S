.syntax unified
.thumb


/*
    R0 = ret_ptr : bf128_t*
    R1 = x : uint8_t
 */

@ NOTES
@ Implementation using conditional move, in a high level sense
@ We dont use conditional execution, but arthimetic to handle the condition
.global bf128_byte_combine_bits_sclf
bf128_byte_combine_bits_sclf:
    push    {r4, r5, r6, r7, r8, r9, r10, r11, r12, lr}
    sub     sp, #40

    mov     r9, r1
    str     r0, [sp, #0]
    mov     r0, #0
    mov     r1, #0
    mov     r2, #0
    mov     r3, #0

    mov     r1, #16
    add     r0, sp, #20
    //bl      rand_mask


    // use r5-r8 as result registers
    mov     r5, #0
    mov     r6, #0
    mov     r7, #0
    mov     r8, #0
    ldr     r5, [sp, #20]
    ldr     r6, [sp, #24]
    ldr     r7, [sp, #28]
    ldr     r8, [sp, #32]


    // bit counter
    mov     r4, #1

    mov     r1, #4
    add     r0, sp, #36
    //bl      rand_mask
    ldr     r11, [sp, #36]
    
    //ldr     r11, =#0x5aa5a55a


    @ Start looping at 1300 sample points

bit_loopr:
    // load alpha onto stack
    add     r0, sp, #4
    sub     r1, r4, #1
    bl      bf128_alpha_wrapper

    // load all four alpha words
    ldr     r0, [sp, #4]
    ldr     r1, [sp, #8]
    ldr     r2, [sp, #12]
    ldr     r3, [sp, #16]

    // Apply alpha
    eor     r0, r0, r5
    eor     r1, r1, r6
    eor     r2, r2, r7
    eor     r3, r3, r8

    // extract bit from x
    lsr     r10, r9, r4
    and     r10, r10, #1

    // conditional move based on bit
    // X+Y
    add     lr, r0, r5

    // (b+r)*Y
    add     r12, r11, r10
    mul     r0, r0, r12   

    // ((1-b)+r) * X
    add     r12, r11, #1
    sub     r12, r12, r10
    mla     r0, r5, r12, r0

    // -r*(X+Y)
    mls     r5, r11, lr, r0


    // X+Y
    add     lr, r1, r6

    // (b+r)*Y
    add     r12, r11, r10
    mul     r1, r1, r12   

    // ((1-b)+r) * X
    add     r12, r11, #1
    sub     r12, r12, r10
    mla     r1, r6, r12, r1

    // -r*(X+Y)
    mls     r6, r11, lr, r1


    // X+Y
    add     lr, r2, r7

    // (b+r)*Y
    add     r12, r11, r10
    mul     r2, r2, r12   

    // ((1-b)+r) * X
    add     r12, r11, #1
    sub     r12, r12, r10
    mla     r2, r7, r12, r2

    // -r*(X+Y)
    mls     r7, r11, lr, r2



        // X+Y
    add     lr, r3, r8

    // (b+r)*Y
    add     r12, r11, r10
    mul     r3, r3, r12   

    // ((1-b)+r) * X
    add     r12, r11, #1
    sub     r12, r12, r10
    mla     r3, r8, r12, r3

    // -r*(X+Y)
    mls     r8, r11, lr, r3

    cmp     r4, #8
    add     r4, r4, #1
    bne     bit_loopr


    // Add lsb
    and     r9, r9, #1
    eor     r5, r5, r9

    ldr     r1, [sp, #20]
    ldr     r2, [sp, #24]
    ldr     r3, [sp, #28]
    ldr     r4, [sp, #32]
    ldr     r0, [sp, #0]

    eor     r5, r5, r1
    eor     r6, r6, r2
    eor     r7, r7, r3
    eor     r8, r8, r4

    // store result
    str     r5, [r0]
    str     r6, [r0, #4]
    str     r7, [r0, #8]
    str     r8, [r0, #12]

    mov     r1, #0
    mov     r2, r1
    mov     r3, r1
    mov     r4, r1
    mov     r5, r1
    mov     r6, r1
    mov     r7, r1
    mov     r8, r1
    mov     r9, r1
    mov     r10, r1
    mov     r11, r1
    mov     r12, r1

    //clear_stack 36
    add     sp, #40
    pop     {r4, r5, r6, r7, r8, r9, r10, r11, r12, pc}