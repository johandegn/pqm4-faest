.syntax unified
.thumb

.align 8
.macro clear_stack, size, label
    mov     r1, #\size
    mov     r0, #0
\label:
    subs    r1, #4
    str     r0, [sp, r1]
    bne     \label
    ldr     r0, [sp]
    mov     r2, #0
    mov     r3, #0
.endm

.global compute_sbox_masked
compute_sbox_masked:
    push    {r4, r5, r6, lr}
    sub     sp, #4
    mov     r6, #0
    str     r6, [sp]

    mov     r5, r2
    mov     r4, r3
    bl      bf8_inv_masked
    
    strh    r6, [sp]    //CLEAR
    mov     r3, #0      //CLEAR
    mov     r0, r6      //CLEAR
    ldrb    r1, [sp]    //CLEAR

    ldrb    r0, [r5]
    bl      affine
    strb    r0, [r5]
    
    strh    r6, [sp]    //CLEAR
    mov     r3, #0      //CLEAR
    mov     r0, r6      //CLEAR
    ldrb    r1, [sp]    //CLEAR

    ldrb    r0, [r4]
    bl      affine_incomplete
    strb    r0, [r4]

    str     r6, [sp]    //CLEAR
    mov     r3, #0      //CLEAR
    mov     r0, r6      //CLEAR
    ldrb    r1, [sp]    //CLEAR

    add     sp, #4
    pop     {r4, r5, r6, pc}


.global aes_extend_witness_masked_output
aes_extend_witness_masked_output:
    adds    r2, #7
    lsrs    r2, r2, #3
    beq     aes_extend_witness_masked_output_early_return
    push    {r4, r5, lr}
    adds    r5, r0, r2
    add     ip, r0, #4294967295     @ 0xffffffff
    add     lr, r0, r2, lsl #1
    mov     r3, r5
aes_extend_witness_masked_output_loop: 
    ldr     r4, [r1, #0]
    subs    r2, r3, r5
    ldrb    r4, [r4, r2]
    strb    r4, [ip, #1]!
    strb    r1, [r3]                // CLEAR
    ldr     r4, [r1, #4]
    ldrb    r2, [r4, r2]
    strb    r2, [r3], #1
    cmp     r3, lr
    bne     aes_extend_witness_masked_output_loop
    mov     r4, #7                  // CLEAR
    pop     {r4, r5, pc}
aes_extend_witness_masked_output_early_return:
    bx      lr


.global shift_row
shift_row:
    push    {r4, r5, r6, lr}
    sub     sp, #32
    mov     r4, r1
    movs    r2, #32
    movs    r1, #0
    mov     r5, r0
    mov     r0, sp
    bl      memset
    cmp     r4, #6
    beq     sr_l2
    cmp     r4, #8
    beq     sr_l1
    cmp     r4, #4
    beq     sr_l2
    cbz     r4, sr_exit
sr_rexit:
    lsls    r2, r4, #2
    mov     r1, sp
    mov     r0, r5
    bl      memcpy
sr_exit:
    // CLEAR Stack
    mov     r1, #0
    str     r1, [sp]
    str     r1, [sp, #4]
    str     r1, [sp, #8]
    str     r1, [sp, #12]
    str     r1, [sp, #16]
    str     r1, [sp, #20]
    str     r1, [sp, #24]
    str     r1, [sp, #28]
    add     sp, #32
    pop     {r4, r5, r6, pc}
sr_l1:
    sub     ip, sp, #4
    movs    r2, #0
sr_loop:
    ldrb    lr, [r5, r2, lsl #2]
    mov     r3, r2
    adds    r2, #1
    adds    r0, r3, #3
    and     r1, r2, #7
    add     r1, r5, r1, lsl #2
    and     r0, r0, #7
    adds    r3, #4
    add     r0, r5, r0, lsl #2
    and     r3, r3, #7
    ldrb    r6, [r1, #1]
    ldrb    r0, [r0, #2]
    movs    r1, #0
    add     r3, r5, r3, lsl #2
    bfi     r1, lr, #0, #8
    ldrb    r3, [r3, #3]
    bfi     r1, r6, #8, #8
    bfi     r1, r0, #16, #8
    bfi     r1, r3, #24, #8
    cmp     r2, #8
    str     r1, [ip, #4]!
    bne     sr_loop
    b       sr_rexit
sr_l2:
    sub     ip, sp, #4
    movs    r2, #0
sr_loop2:
    ldrb    lr, [r5, r2, lsl #2]
    mov     r3, r2
    adds    r2, #1
    adds    r1, r3, #2
    udiv    r0, r2, r4
    mls     r0, r4, r0, r2
    adds    r3, #3
    add     r0, r5, r0, lsl #2
    udiv    r6, r1, r4
    mls     r1, r4, r6, r1
    add     r1, r5, r1, lsl #2
    udiv    r6, r3, r4
    mls     r3, r4, r6, r3
    add     r3, r5, r3, lsl #2
    ldrb    r6, [r0, #1]
    ldrb    r1, [r1, #2]
    ldrb    r3, [r3, #3]
    movs    r0, #0
    bfi     r0, lr, #0, #8
    bfi     r0, r6, #8, #8
    bfi     r0, r1, #16, #8
    bfi     r0, r3, #24, #8
    cmp     r2, r4
    str     r0, [ip, #4]!
    bcc     sr_loop2
    b       sr_rexit

.global mix_column
mix_column:
    cmp     r1, #0
    beq     mc_fexit
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #20
    add     r3, r0, r1, lsl #2
    mov     r4, r0
    str     r3, [sp, #12]
mc_loop:
    ldrb    r7, [r4, #0]
    movs    r1, #2
    mov     r0, r7
    bl      bf8_mul
    ldrb    fp, [r4, #1]
    str     r0, [sp, #8]
    movs    r1, #3
    mov     r0, fp
    bl      bf8_mul
    movs    r1, #2
    mov     sl, r0
    mov     r0, fp
    bl      bf8_mul
    ldrb    r9, [r4, #2]
    ldrb    r8, [r4, #3]
    mov     r3, r0
    movs    r1, #3
    mov     r0, r9
    str     r3, [sp, #4]
    bl      bf8_mul
    movs    r6, #0
    mov     ip, r0
    movs    r1, #2
    mov     r0, r9
    str     ip, [sp]
    mov     r5, r6
    bl      bf8_mul
    bfi     r6, r9, #0, #8
    bfi     r6, r7, #8, #8
    bfi     r5, r8, #0, #8
    bfi     r6, r7, #16, #8
    bfi     r5, r8, #8, #8
    bfi     r5, fp, #16, #8
    bfi     r6, fp, #24, #8
    movs    r1, #3
    mov     fp, r0
    mov     r0, r8
    bl      bf8_mul
    bfi     r5, r9, #24, #8
    movs    r1, #3
    mov     r9, r0
    mov     r0, r7
    bl      bf8_mul
    ldr     r2, [sp, #8]
    ldr     r3, [sp, #4]
    movs    r7, #0
    bfi     r7, r2, #0, #8
    bfi     r7, r3, #8, #8
    mov     ip, r0
    bfi     r7, fp, #16, #8
    movs    r1, #2
    mov     r0, r8
    bfi     r7, ip, #24, #8
    bl      bf8_mul
    ldr     r2, [sp, #0]
    movs    r3, #0
    bfi     r3, sl, #0, #8
    bfi     r3, r2, #8, #8
    bfi     r3, r9, #16, #8
    bfi     r3, r0, #24, #8
    eors    r3, r7
    eors    r3, r6
    eors    r3, r5
    str     r3, [r4], #4
    ldr     r3, [sp, #12]
    cmp     r4, r3
    bne     mc_loop
    // CLEAR STACK
    mov     r1, #0
    str     r1, [sp]
    str     r1, [sp, #4]
    str     r1, [sp, #8]
    str     r1, [sp, #12]
    str     r1, [sp, #16]
    add     sp, #20
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
mc_fexit:
    bx      lr


.global aes_encrypt_round_masked
aes_encrypt_round_masked:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #20
    ldrd    r5, sl, [sp, #56]       @ 0x38
    add     r9, r0, #32
    ldr     fp, [sp, #64]   @ 0x40
    mov     r4, r1
    mov     r6, r2
    mov     r2, r1
    mov     r1, r9
    mov     r7, r3
    mov     r8, r0
    bl      sub_bytes_masked
    mov     r0, #7                              // CLEAR
    mov     r2, #7                              // CLEAR
    strd    r0, r2, [sp, #4]                    // CLEAR
    strd    sl, fp, [sp, #4]
    str     r5, [sp, #0]
    ldr     r3, [r7, #0]
    mov     r0, r8
    mov     r2, r6
    mov     r1, r4
    bl      aes_encrypt_round_masked_inner
    mov     r0, #7                              // CLEAR
    mov     r2, #7                              // CLEAR
    strd    r0, r2, [sp, #4]                    // CLEAR
    strd    sl, fp, [sp, #4]
    add     r2, r6, #480    @ 0x1e0
    ldr     r3, [r7, #4]
    str     r5, [sp, #0]
    mov     r0, r9
    mov     r1, r4
    bl      aes_encrypt_round_masked_inner
    mov     r0, #7                              // CLEAR
    mov     r2, #7                              // CLEAR
    strd    r0, r2, [sp, #4]                    // CLEAR
    ldr     r3, [r5, #0]
    add     r3, r3, r4, lsl #2
    str     r3, [r5, #0]
    // CLEAR STACK
    mov     r1, #0
    str     r1, [sp]
    str     r1, [sp, #4]
    str     r1, [sp, #8]
    str     r1, [sp, #12]
    str     r1, [sp, #16]
    add     sp, #20
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global bf8_mul
bf8_mul:
    push    {lr}
    lsrs    r3, r1, #7
    add     r3, r3, r3, lsl #1
    add     r3, r3, r3, lsl #3
    sub     sp, #12
    eor     r3, r3, r1, lsl #1
    mov     r2, r1
    uxtb    r3, r3
    strb    r2, [sp, #5]
    strb    r3, [sp, #6]
    eors    r2, r3
    add     r3, sp, #8
    movs    r1, #0
    add     r3, r3, r0, lsr #6
    strb    r1, [sp, #4]
    strb    r2, [sp, #7]
    ldrb    r2, [r3, #-4]
    mov     ip, #4
 bf8_mul_loop:
    asr     r3, r0, ip
    and     r3, r3, #3
    add     r3, sp
    mov     lr, r2, asr #6
    ldrb    r3, [r3, #4]
    sub     ip, ip, #2
    orr     r1, lr, r1, lsl #2
    eor     r2, r3, r2, lsl #2
    cmn     ip, #2
    uxtb    r1, r1
    uxtb    r2, r2
    bne     bf8_mul_loop
    mov     r0, sp
    strb    r2, [sp]
    strb    r1, [sp, #1]
    bl      bf8_reduce
    // CLEAR STACK
    mov     r1, #0
    str     r1, [sp]
    str     r1, [sp, #4]
    str     r1, [sp, #8]
    add     sp, #12
    ldr     pc, [sp], #4


.global expand_128key_masked
expand_128key_masked:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, lr}
    sub     sp, #8
    mov     r6, r0
    ldr     r5, [sp, #40]   @ 0x28
    mov     r4, r3
    mov     r7, r2
    cbz     r2, expand_128key_masked_l1
    mov     r2, r1
    movs    r0, #0
expand_128key_masked_l2:
    udiv    ip, r0, r4
    mls     r3, r4, ip, r0
    add     r3, r3, ip, lsl #3
    ldrb    ip, [r2]
    strb    ip, [r6, r3, lsl #2]
    add     r3, r6, r3, lsl #2
    ldrb    ip, [r2, #1]
    strb    ip, [r3, #1]
    ldrb    ip, [r2, #2]
    strb    ip, [r3, #2]
    mov     ip, r0
    adds    r0, #1
    ldrb    lr, [r2, #3]
    strb    lr, [r3, #3]
    cmp     r7, r0
    add     r2, r2, #4
    bne     expand_128key_masked_l2
    movs    r2, #0
expand_128key_masked_l3:
    udiv    r0, r2, r4
    mls     r3, r4, r0, r2
    add     r3, r3, r0, lsl #3
    add     r3, r6, r3, lsl #2
    ldrb    r0, [r1, #32]
    strb    r0, [r3, #480]  @ 0x1e0
    ldrb    r0, [r1, #33]   @ 0x21
    strb    r0, [r3, #481]  @ 0x1e1
    ldrb    r0, [r1, #34]   @ 0x22
    strb    r0, [r3, #482]  @ 0x1e2
    ldrb    lr, [r1, #35]   @ 0x23
    strb    lr, [r3, #483]  @ 0x1e3
    cmp     ip, r2
    add     r1, r1, #4
    add     r2, r2, #1
    bne     expand_128key_masked_l3
expand_128key_masked_l1:
    mla     r9, r5, r4, r4
    cmp     r7, r9
    bcs     expand_128key_masked_l4_step
    mov     r5, r7
    add     r8, r6, #480    @ 0x1e0
    b       expand_128key_masked_l5
expand_128key_masked_l6:
    uxtb    r1, r3
    mov     r3, #0                          // CLEAR
    subs    r3, r5, r7
    udiv    r2, r3, r4
    mls     r3, r4, r2, r3
    add     r3, r3, r2, lsl #3
    ldrb    lr, [sp, #1]
    ldrb    ip, [sp, #2]
    udiv    r0, r5, r4
    movs    r2, #0
    bfi     r2, r1, #0, #8
    mov     r1, #0                          // CLEAR
    bfi     r2, lr, #8, #8
    mov     lr, #0                          // CLEAR
    ldrb    lr, [sp, #3]
    mls     r1, r4, r0, r5
    bfi     r2, ip, #16, #8
    mov     ip, #0                          // CLEAR
    ldr     ip, [r6, r3, lsl #2]
    add     r1, r1, r0, lsl #3
    bfi     r2, lr, #24, #8
    mov     lr, #0                          // CLEAR
    eor     r2, r2, ip
    str     r2, [r6, r1, lsl #2]
    nop                                     // CLEAR
    ldr     ip, [sp, #40]                   // CLEAR
    nop                                     // CLEAR
    mov     r2, #0                          // CLEAR
    nop                                     // CLEAR
    str     r6, [sp, #0]                    // CLEAR
    nop                                     // CLEAR
    ldr     r2, [sp, #4]
    ldr     ip, [sp, #40]                   // CLEAR
    ldr     r3, [r8, r3, lsl #2]
    ldr     ip, [sp, #40]                   // CLEAR
    adds    r5, #1
    eors    r3, r2
    cmp     r5, r9
    str     r3, [r8, r1, lsl #2]
    str     r6, [sp, #4]                    // CLEAR
    bne     expand_128key_masked_l5
expand_128key_masked_l4_step:
    b       expand_128key_masked_l4
expand_128key_masked_l5:
    mov     r3, #0                          // CLEAR
    nop                                     // CLEAR
    mov     r1, #0                          // CLEAR
    nop                                     // CLEAR
    mov     r2, #0                          // CLEAR
    nop                                     // CLEAR
    mov     lr, #0                          // CLEAR
    nop                                     // CLEAR
    mov     ip, #0                          // CLEAR
    nop                                     // CLEAR
    mov     r0, #0                          // CLEAR
    str     r0, [sp]                        // CLEAR
    str     r0, [sp, #4]                    // CLEAR
    udiv    sl, r5, r7
    subs    r3, r5, #1
    udiv    r2, r3, r4
    mls     r3, r4, r2, r3
    add     r3, r3, r2, lsl #3
    add     r2, r3, #120    @ 0x78
    ldr     r3, [r6, r3, lsl #2]
    ldr     r1, [sp, #40]                   // CLEAR
    str     r3, [sp, #0]
    str     r6, [sp, #4]                    // CLEAR
    nop                                     // CLEAR
    ldr     r1, [r6, r2, lsl #2]
    ldr     r2, [sp, #40]                   // CLEAR
    str     r1, [sp, #4]
    str     r2, [sp, #40]                   // CLEAR
    mls     r2, r7, sl, r5
    cmp     r2, #0
    bne     expand_128key_masked_l6
    mov     r0, sp
    bl      rot_word
    add     r0, sp, #4
    bl      rot_word
    mov     r0, sp
    bl      sub_words_masked
    add     r0, sl, #4294967295     @ 0xffffffff
    uxtb    r0, r0
    bl      round_constants
    ldrb    r3, [sp]
    eors    r3, r0
    b       expand_128key_masked_l6
expand_128key_masked_l4:
    // CLEAR STACK
    mov     r3, #0
    nop
    mov     r2, #0
    nop
    mov     r0, #0
    nop
    mov     r1, #0
    str     r1, [sp]
    str     r1, [sp, #4]
    add     sp, #8
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, pc}


.global rot_word
rot_word:
    sub     sp, #4

    ldrb    ip, [r0, #1]
    ldrb    r2, [r0, #2]
    ldrb    r1, [r0, #3]
    movs    r3, #0
    bfi     r3, ip, #0, #8
    bfi     r3, r2, #8, #8
    ldrb    r2, [r0, #0]
    bfi     r3, r1, #16, #8
    bfi     r3, r2, #24, #8
    str     r3, [r0, #0]

    mov     r0, #0
    str     r0, [sp]                        // CLEAR
    mov     ip, #0                          // CLEAR
    mov     r3, #0                          // CLEAR
    mov     r2, #0                          // CLEAR
    mov     r1, #0                          // CLEAR
    add     sp, #4
    bx      lr


.global aes_key_schedule_backward_128_vbb_vk_round
aes_key_schedule_backward_128_vbb_vk_round:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #236        @ 0xec
// CLEAR STACK
    mov     r4, #236
    mov     r5, #0
aes_key_init_clear_loop:
    subs    r4, #4
    str     r5, [sp, r4]
    bne     aes_key_init_clear_loop
// END OF CLEAR
    mov     r4, r2
    movs    r2, #0
    mov     r8, r3
    ldr     r3, [sp]                        // CLEAR
    movs    r3, #0
    str     r2, [sp]                        // CLEAR
    strd    r2, r3, [sp, #48]       @ 0x30
    strd    r2, r3, [sp, #40]       @ 0x28
    lsls    r3, r4, #7
    lsrs    r3, r3, #2
    str     r3, [sp, #24]
    add     r3, sp, #48     @ 0x30
    mov     r6, r0
    mov     r7, r1
    ldmia   r3, {r0, r1}
    stmia   sp, {r0, r1}
    movs    r5, #0
    str     r5, [sp, #8]
    add     r3, sp, #40     @ 0x28
    ldmia   r3, {r2, r3}
    add     r0, sp, #56     @ 0x38
    lsls    r4, r4, #3
    bl      bf128_mul_bit
    add     r3, r4, #128    @ 0x80
    str     r3, [sp, #28]
    ldr     r3, [sp, #56]   @ 0x38
    str     r5, [sp, #20]
    eor     lr, r3, #1
    add     r3, sp, #56     @ 0x38
    ldmia   r3, {r0, r1, r2, r3}
    add     r5, sp, #104    @ 0x68
    stmia   r5, {r0, r1, r2, r3}
    str     lr, [sp, #104]  @ 0x68
    ldmia   r5, {r0, r1, r2, r3}
    add     r9, sp, #56     @ 0x38
    ldr     ip, [sp, #60]   @ 0x3c
    stmia   r9, {r0, r1, r2, r3}
    strd    ip, lr, [sp, #32]
    add     r4, sp, #72     @ 0x48
    mov     sl, r6
aes_key_loop:
    ldr     r3, [sp, #28]
    ldr     r2, [sp, #20]
    mov     fp, r5
    add     r9, r3, r2
    mov     r6, r5
aes_key_l1:
    mov     r1, r9
    mov     r2, r8
    mov     r0, sl
    bl      get_vole_aes_128_share
    ldmia   r0, {r0, r1, r2, r3}
    stmia   r4, {r0, r1, r2, r3}
    mov     r1, r9
    mov     r3, r8
    mov     r2, r4
    mov     r0, sl
    bl      add_vole_to_vk_cache_share
    ldmia   r4, {r0, r1, r2, r3}
    stmia   r6, {r0, r1, r2, r3}
    adds    r6, #16
    add     r3, sp, #232    @ 0xe8
    cmp     r6, r3
    add     r9, r9, #1
    bne     aes_key_l1
    ldrd    r2, r3, [sp, #20]
    adds    r6, r3, r2
aes_key_l2:
    mov     r1, r6
    mov     r2, r8
    mov     r0, sl
    bl      get_vk_128_share
    ldrd    r1, r3, [r0, #8]
    ldrd    ip, r2, [r0]
    ldr     r0, [fp, #8]
    eors    r1, r0
    str     r1, [fp, #8]
    ldr     r1, [fp, #12]
    eors    r3, r1
    str     r3, [fp, #12]
    ldr     r3, [fp]
    eor     r3, r3, ip
    str     r3, [fp]
    ldr     r3, [fp, #4]
    eors    r3, r2
    str     r3, [fp, #4]
    add     fp, fp, #16
    add     r3, sp, #232    @ 0xe8
    cmp     fp, r3
    add     r6, r6, #1
    bne     aes_key_l2
    mov     r2, r7
    movs    r3, #7
    mov     fp, r7
aes_key_l3:
    subs    r6, r3, #2
    and     r0, r3, #7
    and     r6, r6, #7
    add     ip, r5, r6, lsl #4
    add     r9, r5, r0, lsl #4
    lsls    r6, r6, #4
    lsls    r0, r0, #4
    ldr     r6, [r5, r6]
    ldr     r0, [r5, r0]
    ldr     r7, [r9, #8]
    subs    r1, r3, #5
    and     r1, r1, #7
    eors    r0, r6
    add     r6, r5, r1, lsl #4
    lsls    r1, r1, #4
    adds    r3, #1
    ldr     r1, [r5, r1]
    eors    r0, r1
    ldr     r1, [ip, #8]
    eor     lr, r7, r1
    ldr     r7, [ip, #12]
    ldr     r1, [r9, #12]
    ldr     r9, [r9, #4]
    eors    r1, r7
    ldr     r7, [ip, #4]
    eor     ip, r9, r7
    ldr     r7, [r6, #8]
    eor     lr, lr, r7
    ldr     r7, [r6, #12]
    eors    r7, r1
    ldr     r1, [r6, #4]
    str     r0, [r2, #0]
    eor     r1, ip, r1
    cmp     r3, #15
    strd    r1, lr, [r2, #4]
    str     r7, [r2, #12]
    add     r2, r2, #16
    bne     aes_key_l3
    ldr     r3, [sp, #36]   @ 0x24
    str     r3, [sp, #56]   @ 0x38
    movs    r6, #0
    ldr     r3, [sp, #32]
    str     r3, [sp, #60]   @ 0x3c
    str     r6, [sp, #8]
    add     r3, sp, #64     @ 0x40
    ldmia   r3, {r0, r1}
    stmia   sp, {r0, r1}
    add     r3, sp, #56     @ 0x38
    ldmia   r3, {r2, r3}
    add     r0, sp, #88     @ 0x58
    bl      bf128_mul_bit
    str     r6, [sp, #8]
    ldrd    r1, r2, [fp, #8]
    ldr     r0, [sp, #96]   @ 0x60
    eors    r1, r0
    str     r1, [fp, #8]
    ldr     r1, [sp, #100]  @ 0x64
    ldrd    r3, r6, [fp]
    eors    r2, r1
    str     r2, [fp, #12]
    ldr     r2, [sp, #88]   @ 0x58
    eors    r3, r2
    str     r3, [fp]
    ldr     r3, [sp, #92]   @ 0x5c
    eors    r6, r3
    add     r3, sp, #64     @ 0x40
    ldmia   r3, {r0, r1}
    str     r6, [fp, #4]
    stmia   sp, {r0, r1}
    add     r3, sp, #56     @ 0x38
    ldmia   r3, {r2, r3}
    add     r0, sp, #88     @ 0x58
    bl      bf128_mul_bit
    ldr     r3, [sp, #20]
    ldrd    r0, r1, [fp, #40]       @ 0x28
    adds    r3, #8
    str     r3, [sp, #20]
    ldr     r3, [sp, #96]   @ 0x60
    eors    r3, r0
    str     r3, [fp, #40]   @ 0x28
    ldr     r3, [sp, #100]  @ 0x64
    ldrd    r2, r6, [fp, #32]
    eors    r3, r1
    str     r3, [fp, #44]   @ 0x2c
    ldr     r3, [sp, #88]   @ 0x58
    eors    r3, r2
    str     r3, [fp, #32]
    ldr     r3, [sp, #92]   @ 0x5c
    eors    r3, r6
    str     r3, [fp, #36]   @ 0x24
    ldr     r3, [sp, #20]
    mov     r7, fp
    cmp     r3, #32
    add     r7, r7, #128    @ 0x80
    bne     aes_key_loop
// CLEAR STACK
    mov     r1, #236
    mov     r0, #0
aes_key_clear_loop:
    subs    r1, #4
    str     r0, [sp, r1]
    bne     aes_key_clear_loop
    ldr     r0, [sp]
    mov     r2, #0
    mov     r3, #0
// END OF CLEAR
    add     sp, #236        @ 0xec
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global bf128_mul
bf128_mul:
    sub     sp, #8
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #300        @ 0x12c
    add     r1, sp, #336    @ 0x150
    stmia   r1, {r2, r3}
    add     r3, sp, #352    @ 0x160
    mov     r8, r0
    ldmia   r3, {r0, r1, r2, r3}
    add     r6, sp, #56     @ 0x38
    stmia   r6, {r0, r1, r2, r3}
    movs    r4, #0
    movs    r5, #0
    movs    r7, #0
    strd    r4, r5, [sp, #40]       @ 0x28
    strd    r4, r5, [sp, #48]       @ 0x30
    mov     ip, r7
    add     r5, sp, #40     @ 0x28
    add     r4, sp, #8
    str     r8, [sp, #4]
bf128_mul_l3:
    movs    r3, #1
    lsl     r3, r3, ip
    add     r2, r5, r3, lsl #4
    ldrd    r1, r7, [r2, #8]
    lsls    r3, r3, #4
    adds    r1, r1, r1
    ldr     r3, [r5, r3]
    ldr     r0, [r2, #4]
    mov     r6, r7, asr #31
    and     r6, r6, #135    @ 0x87
    adcs    r7, r7
    adds    r3, r3, r3
    eor     r6, r6, r3
    ldr     r3, [r5, #24]
    orr     r1, r1, r0, lsr #31
    eor     r8, r1, r3
    ldr     r3, [r5, #28]
    eor     r9, r7, r3
    ldr     r3, [r5, #16]
    mov     r2, #2
    eor     sl, r3, r6
    ldr     r3, [r5, #20]
    lsl     r2, r2, ip
    adcs    r0, r0
    add     lr, r2, #1
    eor     fp, r3, r0
    add     r3, r5, r2, lsl #4
    strd    r1, r7, [r3, #8]
    lsls    r2, r2, #4
    add     r1, r5, lr, lsl #4
    mov     lr, lr, lsl #4
    str     r6, [r5, r2]
    strd    r8, r9, [r1, #8]
    str     r0, [r3, #4]
    str     sl, [r5, lr]
    str     fp, [r1, #4]
    cmp     ip, #0
    beq     bf128_mul_l1
    add     lr, r4, r2
    adds    r2, #8
    mov     r1, r5
    add     r2, r5
bf128_mul_l2:
    ldrd    r6, r7, [r1, #40]       @ 0x28
    ldr     r0, [r2, #0]
    ldr     r9, [r2, #-4]
    eors    r6, r0
    ldr     r0, [r2, #4]
    eor     r8, r0, r7
    ldr     r0, [r1, #32]
    ldr     r7, [r2, #-8]
    eors    r0, r7
    ldr     r7, [r1, #36]   @ 0x24
    str     r0, [r3, #32]
    adds    r1, #16
    eor     r7, r7, r9
    cmp     lr, r1
    strd    r7, r6, [r3, #36]       @ 0x24
    str     r8, [r3, #44]   @ 0x2c
    add     r3, r3, #16
    bne     bf128_mul_l2
bf128_mul_l1:
    add     ip, ip, #1
    cmp     ip, #3
    bne     bf128_mul_l3
    ldr     r8, [sp, #4]
    movs    r2, #32
    movs    r1, #0
    mov     r0, r4
    movs    r6, #0
    bl      memset
    mov     lr, #60 @ 0x3c
    mov     ip, r6
    str     r8, [sp, #4]
bf128_mul_l6:
    ldr     r9, [r4]
    ldr     r2, [sp, #12]
    add     r0, sp, #336    @ 0x150
    mov     r1, r4
    rsb     r8, lr, #32
    sub     r7, lr, #32
bf128_mul_l4:
    ldrd    r3, sl, [r0]
    lsl     fp, sl, r8
    lsr     r3, r3, lr
    lsr     sl, sl, r7
    orr     r3, r3, fp
    orr     r3, r3, sl
    and     r3, r3, #15
    add     sl, r5, r3, lsl #4
    lsls    r3, r3, #4
    adds    r0, #8
    ldr     r3, [r5, r3]
    eor     r3, r3, r9
    str     r3, [r1, #0]
    ldr     r3, [sl, #4]
    eors    r3, r2
    str     r3, [r1, #4]
    ldr     r2, [r1, #8]!
    ldr     r3, [sl, #8]
    eor     r9, r2, r3
    ldr     r3, [r1, #4]
    ldr     r2, [sl, #12]
    str     r9, [r1]
    eors    r2, r3
    add     r3, sp, #352    @ 0x160
    cmp     r3, r0
    str     r2, [r1, #4]
    bne     bf128_mul_l4
    sub     lr, lr, #4
    cmn     lr, #4
    beq     bf128_mul_l5
    ldr     r1, [r4, #20]
    ldr     r2, [r4, #12]
    ldr     r8, [r4, #16]
    ldr     r3, [r4, #4]
    ldr     r7, [r4, #8]
    lsrs    r0, r1, #28
    mov     ip, ip, lsl #4
    lsls    r1, r1, #4
    orr     r1, r1, r8, lsr #28
    orr     ip, ip, r6, lsr #28
    orr     r6, r0, r6, lsl #4
    lsrs    r0, r2, #28
    orr     r0, r0, r8, lsl #4
    str     r1, [r4, #20]
    ldr     r1, [r4, #0]
    str     r0, [r4, #16]
    lsls    r2, r2, #4
    lsrs    r0, r3, #28
    lsls    r3, r3, #4
    orr     r3, r3, r1, lsr #28
    orr     r2, r2, r7, lsr #28
    orr     r0, r0, r7, lsl #4
    lsls    r1, r1, #4
    strd    r0, r2, [r4, #8]
    strd    r1, r3, [r4]
    b       bf128_mul_l6
bf128_mul_l5:
    ldr     r5, [sp, #24]
    ldr     r2, [sp, #16]
    ldr     r7, [sp, #20]
    ldr     r0, [sp, #28]
    ldr     r1, [sp, #12]
    ldr     r8, [sp, #4]
    mov     r3, ip, lsr #30
    eor     r3, r3, ip, lsr #31
    adds    lr, r6, r6
    eor     r3, r3, r5
    eor     r2, r2, r6
    eor     r3, r3, ip, lsr #25
    ldr     r5, [r4, #0]
    eor     r2, r2, lr
    eor     r4, ip, r7
    adc     lr, ip, ip
    eor     r2, r2, r6, lsl #2
    adds    r9, r3, r3
    mov     r7, ip, lsl #2
    mov     ip, ip, lsl #7
    orr     r7, r7, r6, lsr #30
    orr     ip, ip, r6, lsr #25
    eor     r2, r2, r6, lsl #7
    eor     r5, r5, r3
    adc     r6, r0, r0
    eors    r1, r0
    eor     r4, r4, lr
    mov     lr, r0, lsl #2
    eor     r5, r5, r9
    eor     r2, r2, r0, lsr #31
    eors    r1, r6
    orr     lr, lr, r3, lsr #30
    lsls    r6, r0, #7
    eors    r4, r7
    eor     r2, r2, r0, lsr #30
    eor     r5, r5, r3, lsl #2
    eor     r1, r1, lr
    orr     r6, r6, r3, lsr #25
    eor     r2, r2, r0, lsr #25
    eor     r4, r4, ip
    eor     r5, r5, r3, lsl #7
    eors    r1, r6
    mov     r0, r8
    strd    r5, r1, [r8]
    strd    r2, r4, [r8, #8]
    clear_stack 300 bf128_mul_clear_loop
    add     sp, #300        @ 0x12c
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    str     r0, [sp]                        // CLEAR
    str     r2, [sp, #4]                    // CLEAR
    add     sp, #8
    bx      lr


.global bf128_add_wrapper
bf128_add_wrapper:
    push    {r4, lr}
    mov     ip, r0
    mov     lr, r2
    ldmia   r1, {r0, r1, r2, r3}
    sub     sp, #32
    stmia   sp, {r0, r1, r2, r3}
    ldmia   lr, {r0, r1, r2, r3}
    add     r4, sp, #32
    stmdb   r4, {r0, r1, r2, r3}
    ldrd    r2, r0, [sp, #8]
    ldrd    r3, r1, [sp, #24]
    eors    r3, r2
    eors    r1, r0
    str     ip, [ip, #8]                  // CLEAR
    str     r3, [ip, #8]
    str     ip, [ip, #12]                 // CLEAR
    str     r1, [ip, #12]
    ldrd    r2, r3, [sp, #16]
    ldrd    r0, r1, [sp]
    eors    r2, r0
    eors    r3, r1
    str     ip, [ip]                      // CLEAR
    str     ip, [ip, #4]                  // CLEAR
    strd    r2, r3, [ip]
    mov     r0, #0                        // CLEAR
    mov     r1, r0                        // CLEAR
    mov     r2, r0                        // CLEAR
    mov     r3, r0                        // CLEAR
    stmia   sp, {r0, r1, r2, r3}          // CLEAR
    add     sp, #32
    stmdb   sp, {r0, r1, r2, r3}          // CLEAR
    pop     {r4, pc}


.global bf128_mul_wrapper
bf128_mul_wrapper:
    push    {r4, r5, r6, r7, lr}
    sub     sp, #44 @ 0x2c
    mov     r4, r0
    mov     ip, r1
    add     r6, r1, #8
    ldmia   r2, {r0, r1, r2, r3}
    add     lr, sp, #8
    stmia   lr, {r0, r1, r2, r3}
    ldmia   r6, {r0, r1}
    stmia   sp, {r0, r1}
    add     r5, sp, #24
    ldmia   ip, {r2, r3}
    mov     r0, r5

    mov     r6, #0
    mov     r7, #0
    strd    r6, r7, [r4, #8]
    strd    r6, r7, [r4]

    bl      bf128_mul
    ldmia   r5, {r0, r1, r2, r3}
    stmia   r4, {r0, r1, r2, r3}

    mov     r0, #0                        // CLEAR
    mov     r1, r0                        // CLEAR
    mov     r2, r0                        // CLEAR
    mov     r3, r0                        // CLEAR
    strd    r6, r7, [sp]                  // CLEAR
    strd    r6, r7, [sp, #8]              // CLEAR
    strd    r6, r7, [sp, #16]             // CLEAR
    strd    r6, r7, [sp, #24]             // CLEAR
    strd    r6, r7, [sp, #32]             // CLEAR
    str     r6, [sp, #40]                 // CLEAR
    add     sp, #44 @ 0x2c
    pop     {r4, r5, r6, r7, pc}


.global bf128_byte_combine
bf128_byte_combine:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    ldrd    fp, sl, [r1]
    sub     sp, #68 @ 0x44
    ldr     r9, [r1, #8]
    ldr     r8, [r1, #12]
    str     r0, [sp, #28]
    add     r4, r1, #16
    movs    r7, #0
    add     r5, sp, #32
    add     r6, sp, #8
bf128_byte_combine_loop:
    mov     r1, r7
    mov     r0, r5
    bl      bf128_alpha_wrapper
    ldmia   r5, {r0, r1, r2, r3}
    stmia   r6, {r0, r1, r2, r3}
    add     ip, r4, #8
    ldmia   ip, {r0, r1}
    stmia   sp, {r0, r1}
    add     r0, sp, #48     @ 0x30
    ldmia   r4, {r2, r3}
    bl      bf128_mul
    ldrd    r2, r3, [sp, #48]       @ 0x30
    adds    r7, #1
    eor     fp, r2, fp
    eor     sl, r3, sl
    ldr     r2, [sp, #56]   @ 0x38
    ldr     r3, [sp, #60]   @ 0x3c
    cmp     r7, #7
    eor     r9, r2, r9
    add     r4, r4, #16
    eor     r8, r3, r8
    bne     bf128_byte_combine_loop
    ldr     r3, [sp, #28]
    mov     r0, r3
    strd    fp, sl, [r3]
    strd    r9, r8, [r3, #8]
    clear_stack 68 bf128_byte_combine_clear_loop
    add     sp, #68 @ 0x44
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global bf128_byte_combine_vk_share
bf128_byte_combine_vk_share:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #76 @ 0x4c
    mov     r5, r1
    mov     fp, r2
    str     r1, [sp, #28]
    mov     r2, r3
    mov     r4, r0
    mov     r1, fp
    mov     r0, r5
    str     r3, [sp, #24]
    bl      get_vk_128_share
    add     r3, sp, #8
    ldrd    sl, r9, [r0]
    ldrd    r8, r7, [r0, #8]
    strd    r3, r4, [sp, #32]
    add     fp, fp, #1
    movs    r6, #0
    add     r5, sp, #40     @ 0x28
bf128_byte_combine_vk_share_loop:
    ldr     r2, [sp, #24]
    ldr     r0, [sp, #28]
    add     r1, fp, r6
    bl      get_vk_128_share
    mov     r1, r6
    mov     r4, r0
    mov     r0, r5
    bl      bf128_alpha_wrapper
    ldmia   r5, {r0, r1, r2, r3}
    ldr     lr, [sp, #32]
    stmia   lr, {r0, r1, r2, r3}
    add     ip, r4, #8
    ldmia   ip, {r0, r1}
    stmia   sp, {r0, r1}
    add     r0, sp, #56     @ 0x38
    ldmia   r4, {r2, r3}
    bl      bf128_mul
    ldr     r3, [sp, #56]   @ 0x38
    eor     sl, sl, r3
    ldr     r3, [sp, #60]   @ 0x3c
    eor     r9, r9, r3
    ldr     r3, [sp, #64]   @ 0x40
    adds    r6, #1
    eor     r8, r8, r3
    ldr     r3, [sp, #68]   @ 0x44
    cmp     r6, #7
    eor     r7, r7, r3
    bne     bf128_byte_combine_vk_share_loop
    ldr     r4, [sp, #36]   @ 0x24
    mov     r0, r4
    strd    sl, r9, [r4]
    strd    r8, r7, [r4, #8]
    clear_stack 76 bf128_byte_combine_vk_share_clear_loop
    add     sp, #76 @ 0x4c
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global aes_enc_forward_backward_128
aes_enc_forward_backward_128:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, sp, #772    @ 0x304
    movs    r2, #0
    str     r3, [sp, #48]   @ 0x30
    movs    r3, #0
    strd    r2, r3, [sp, #120]      @ 0x78
    mov     r8, #0
    str     r8, [sp, #8]
    add     r9, sp, #120    @ 0x78
    strd    r1, r0, [sp, #28]
    ldmia   r9, {r0, r1}
    strd    r2, r3, [sp, #112]      @ 0x70
    stmia   sp, {r0, r1}
    add     r4, sp, #192    @ 0xc0
    add     r6, sp, #112    @ 0x70
    ldmia   r6, {r2, r3}
    mov     r0, r4
    ldr     r7, [sp, #820]  @ 0x334
    bl      bf128_mul_bit
    str     r8, [sp, #8]
    ldmia   r4, {r0, r1, r2, r3}
    add     fp, sp, #640    @ 0x280
    stmia   fp, {r0, r1, r2, r3}
    add     r5, sp, #128    @ 0x80
    eor     r0, r0, #1
    stmia   r5, {r0, r1, r2, r3}
    str     r0, [sp, #640]  @ 0x280
    ldmia   r9, {r0, r1}
    stmia   sp, {r0, r1}
    mov     r0, r4
    ldmia   r6, {r2, r3}
    bl      bf128_mul_bit
    str     r8, [sp, #8]
    ldmia   r4, {r0, r1, r2, r3}
    add     r6, sp, #208    @ 0xd0
    stmia   fp, {r0, r1, r2, r3}
    eor     r0, r0, #1
    stmia   r6, {r0, r1, r2, r3}
    add     ip, sp, #216    @ 0xd8
    str     r0, [sp, #640]  @ 0x280
    ldmia   ip, {r0, r1}
    stmia   sp, {r0, r1}
    add     sl, sp, #144    @ 0x90
    ldmia   r6, {r2, r3}
    mov     r0, sl
    bl      bf128_mul_bit
    cmp     r7, #9
    beq     aes_enc_forward_backward_128_l1
    ldr     r3, [sl, #4]
    str     r3, [sp, #60]   @ 0x3c
    ldr     r3, [sl]
    str     r3, [sp, #64]   @ 0x40
    ldr     r3, [sl, #12]
    str     r3, [sp, #56]   @ 0x38
    ldr     r3, [sl, #8]
    str     r3, [sp, #68]   @ 0x44
    cmp     r7, #0
    beq     aes_enc_forward_backward_128_l2
aes_enc_forward_backward_128_l9:
    movs    r1, #2
    add     r0, sp, #160    @ 0xa0
    bl      bf128_byte_combine_bits
    adds    r7, #1
    movs    r1, #3
    add     r0, sp, #176    @ 0xb0
    bl      bf128_byte_combine_bits
    ldr     r2, [sp, #28]
    ldr     r6, [sp, #812]  @ 0x32c
    ldr     r8, [sp, #32]
    lsls    r3, r7, #7
    str     r3, [sp, #76]   @ 0x4c
    subs    r3, #128        @ 0x80
    add     r2, r3
    movs    r3, #0
    str     r3, [sp, #72]   @ 0x48
    add     r3, sp, #96     @ 0x60
    str     r3, [sp, #44]   @ 0x2c
    add     r3, sp, #8
    str     r2, [sp, #92]   @ 0x5c
    add     r9, sp, #512    @ 0x200
    add     sl, sp, #448    @ 0x1c0
    str     r3, [sp, #28]
aes_enc_forward_backward_128_l6:
    ldr     r1, [sp, #72]   @ 0x48
    ldr     r3, [sp, #92]   @ 0x5c
    ldr     r2, [sp, #76]   @ 0x4c
    str     r6, [sp, #80]   @ 0x50
    add     r3, r1
    subs    r3, r3, r2
    str     r3, [sp, #52]   @ 0x34
    adds    r3, r1, #7
    subs    r3, r3, r2
    movs    r4, #0
    ldr     r6, [sp, #824]  @ 0x338
    str     r3, [sp, #48]   @ 0x30
    str     r2, [sp, #24]
    mov     r7, r4
aes_enc_forward_backward_128_l5:
    ldr     r3, [sp, #52]   @ 0x34
    ldr     r2, [sp, #24]
    mov     r4, r9
    adds    r5, r3, r2
aes_enc_forward_backward_128_l3:
    mov     r1, r5
    mov     r2, r6
    mov     r0, r8
    bl      get_vole_aes_128_share
    ldmia   r0, {r0, r1, r2, r3}
    stmia   r4, {r0, r1, r2, r3}
    adds    r4, #16
    cmp     r4, fp
    stmia   fp, {r0, r1, r2, r3}
    add     r5, r5, #1
    bne     aes_enc_forward_backward_128_l3
    add     r3, sp, #384    @ 0x180
    add     r0, r3, r7, lsl #4
    ldr     r4, [sp, #48]   @ 0x30
    ldr     r3, [sp, #24]
    mov     r1, r9
    adds    r4, r3, r4
    lsls    r5, r7, #4
    bl      bf128_byte_combine
    ldr     r2, [sp, #24]
    mov     r3, r6
    add     r0, sl, r5
    mov     r1, r8
    ubfx    r5, r4, #3, #2
    bl      bf128_byte_combine_vk_share
    add     r4, r5, r4, lsr #5
    strd    r5, r7, [sp, #32]
    and     r4, r4, #3
    movs    r2, #7
    mov     r3, fp
    str     r4, [sp, #40]   @ 0x28
aes_enc_forward_backward_128_l4:
    subs    r0, r2, #2
    and     r4, r2, #7
    and     r0, r0, #7
    add     ip, r9, r4, lsl #4
    add     lr, r9, r0, lsl #4
    lsls    r4, r4, #4
    lsls    r0, r0, #4
    ldr     r5, [r9, r4]
    ldr     r4, [r9, r0]
    ldr     r0, [lr, #8]
    ldr     r7, [lr, #12]
    ldr     lr, [lr, #4]
    subs    r1, r2, #5
    and     r1, r1, #7
    eors    r4, r5
    add     r5, r9, r1, lsl #4
    lsls    r1, r1, #4
    adds    r2, #1
    ldr     r1, [r9, r1]
    eors    r4, r1
    ldr     r1, [ip, #8]
    eors    r1, r0
    ldr     r0, [ip, #12]
    eors    r7, r0
    ldr     r0, [ip, #4]
    eor     ip, lr, r0
    ldr     r0, [r5, #8]
    eors    r1, r0
    ldr     r0, [r5, #12]
    ldr     r5, [r5, #4]
    str     r4, [r3, #0]
    eors    r0, r7
    eor     r5, ip, r5
    cmp     r2, #15
    strd    r5, r1, [r3, #4]
    str     r0, [r3, #12]
    add     r3, r3, #16
    bne     aes_enc_forward_backward_128_l4
    ldrd    r5, r7, [sp, #32]
    ldr     r4, [sp, #40]   @ 0x28
    ldr     r3, [sp, #816]  @ 0x330
    ldr     r2, [sp, #68]   @ 0x44
    ldr     r1, [sp, #56]   @ 0x38
    ldr     r0, [sp, #64]   @ 0x40
    add     r5, r5, r4, lsl #2
    add     r5, r3, r5, lsl #4
    ldr     r3, [sp, #648]  @ 0x288
    ldr     r4, [sp, #60]   @ 0x3c
    eors    r3, r2
    str     r3, [sp, #648]  @ 0x288
    ldr     r3, [sp, #652]  @ 0x28c
    eors    r3, r1
    str     r3, [sp, #652]  @ 0x28c
    ldr     r3, [sp, #640]  @ 0x280
    eors    r3, r0
    str     r3, [sp, #640]  @ 0x280
    ldr     r3, [sp, #644]  @ 0x284
    eors    r3, r4
    str     r3, [sp, #644]  @ 0x284
    ldr     r3, [sp, #680]  @ 0x2a8
    eors    r3, r2
    str     r3, [sp, #680]  @ 0x2a8
    ldr     r3, [sp, #684]  @ 0x2ac
    eors    r3, r1
    str     r3, [sp, #684]  @ 0x2ac
    ldr     r3, [sp, #672]  @ 0x2a0
    eors    r3, r0
    str     r3, [sp, #672]  @ 0x2a0
    ldr     r3, [sp, #676]  @ 0x2a4
    eors    r3, r4
    ldr     r4, [sp, #44]   @ 0x2c
    str     r3, [sp, #676]  @ 0x2a4
    mov     r1, fp
    mov     r0, r4
    bl      bf128_byte_combine
    ldr     r3, [sp, #24]
    adds    r3, #8
    str     r3, [sp, #24]
    ldmia   r4, {r0, r1, r2, r3}
    adds    r7, #1
    cmp     r7, #4
    stmia   r5, {r0, r1, r2, r3}
    bne     aes_enc_forward_backward_128_l5
    ldr     r3, [sp, #72]   @ 0x48
    ldr     r5, [sp, #28]
    ldr     r6, [sp, #80]   @ 0x50
    adds    r3, #32
    add     r2, sp, #160    @ 0xa0
    str     r3, [sp, #72]   @ 0x48
    ldmia   r2, {r0, r1, r2, r3}
    stmia   r5, {r0, r1, r2, r3}
    add     r2, sp, #392    @ 0x188
    ldmia   r2, {r0, r1}
    stmia   sp, {r0, r1}
    add     r2, sp, #384    @ 0x180
    ldmia   r2, {r2, r3}
    add     r0, sp, #256    @ 0x100
    bl      bf128_mul
    ldr     r3, [sl, #8]
    ldr     r2, [sp, #264]  @ 0x108
    str     r5, [sp, #28]
    eors    r3, r2
    str     r3, [r6, #8]
    ldr     r2, [sl, #12]
    ldr     r3, [sp, #268]  @ 0x10c
    eors    r3, r2
    str     r3, [r6, #12]
    ldr     r2, [sp, #256]  @ 0x100
    ldr     r3, [sl]
    eors    r3, r2
    str     r3, [r6, #0]
    ldr     r2, [sp, #260]  @ 0x104
    ldr     r3, [sl, #4]
    eors    r3, r2
    add     r2, sp, #176    @ 0xb0
    str     r3, [r6, #4]
    ldmia   r2, {r0, r1, r2, r3}
    stmia   r5, {r0, r1, r2, r3}
    add     r2, sp, #408    @ 0x198
    ldmia   r2, {r0, r1}
    stmia   sp, {r0, r1}
    add     r2, sp, #400    @ 0x190
    ldmia   r2, {r2, r3}
    add     r0, sp, #272    @ 0x110
    bl      bf128_mul
    ldrd    r7, lr, [sp, #416]      @ 0x1a0
    ldrd    r2, r4, [sp, #424]      @ 0x1a8
    ldr     r3, [sp, #272]  @ 0x110
    ldr     r0, [sp, #280]  @ 0x118
    ldr     ip, [sp, #436]  @ 0x1b4
    str     r2, [sp, #52]   @ 0x34
    eors    r0, r2
    eor     r2, r7, r3
    ldr     r3, [sp, #276]  @ 0x114
    str     lr, [sp, #88]   @ 0x58
    mov     r5, lr
    ldr     lr, [sp, #432]  @ 0x1b0
    str     r7, [sp, #84]   @ 0x54
    eors    r3, r5
    mov     r7, ip
    ldr     r1, [sp, #284]  @ 0x11c
    str     r4, [sp, #80]   @ 0x50
    eors    r3, r7
    mov     r7, lr
    ldr     lr, [sp, #444]  @ 0x1bc
    str     lr, [sp, #40]   @ 0x28
    eors    r1, r4
    mov     r4, lr
    ldr     lr, [sp, #440]  @ 0x1b8
    str     lr, [sp, #36]   @ 0x24
    mov     r5, lr
    eors    r1, r4
    ldr     r4, [r6, #8]
    str     ip, [sp, #48]   @ 0x30
    eors    r0, r5
    eors    r0, r4
    str     r0, [r6, #8]
    ldr     r0, [r6, #12]
    eors    r1, r0
    str     r1, [r6, #12]
    ldr     r1, [r6, #0]
    eors    r2, r7
    eors    r2, r1
    str     r2, [r6, #0]
    ldr     r2, [r6, #4]
    ldrd    r0, r4, [sp, #392]      @ 0x188
    eors    r3, r2
    str     r3, [r6, #4]
    ldr     r3, [sl, #24]
    eors    r3, r0
    ldrd    r5, r1, [sp, #384]      @ 0x180
    str     r0, [sp, #24]
    str     r3, [r6, #24]
    ldr     r3, [sl, #28]
    str     r1, [sp, #32]
    eors    r3, r4
    str     r3, [r6, #28]
    ldr     r3, [sl, #16]
    ldr     ip, [sp, #28]
    eors    r3, r5
    str     r3, [r6, #16]
    ldr     r3, [sl, #20]
    eors    r3, r1
    add     r1, sp, #160    @ 0xa0
    str     r3, [r6, #20]
    ldmia   r1, {r0, r1, r2, r3}
    stmia   ip, {r0, r1, r2, r3}
    add     r1, sp, #408    @ 0x198
    ldmia   r1, {r0, r1}
    stmia   sp, {r0, r1}
    add     r1, sp, #288    @ 0x120
    mov     r0, r1
    add     r1, sp, #400    @ 0x190
    ldmia   r1, {r2, r3}
    bl      bf128_mul
    ldrd    r0, r1, [r6, #24]
    ldr     r2, [sp, #296]  @ 0x128
    ldr     r3, [r6, #16]
    ldr     ip, [sp, #28]
    eors    r2, r0
    str     r2, [r6, #24]
    ldr     r2, [sp, #300]  @ 0x12c
    eors    r2, r1
    str     r2, [r6, #28]
    ldr     r2, [sp, #288]  @ 0x120
    ldr     r1, [r6, #20]
    eors    r3, r2
    str     r3, [r6, #16]
    ldr     r3, [sp, #292]  @ 0x124
    eors    r3, r1
    str     r3, [r6, #20]
    add     r3, sp, #176    @ 0xb0
    ldmia   r3, {r0, r1, r2, r3}
    stmia   ip, {r0, r1, r2, r3}
    add     r3, sp, #424    @ 0x1a8
    ldmia   r3, {r0, r1}
    stmia   sp, {r0, r1}
    add     r3, sp, #416    @ 0x1a0
    ldmia   r3, {r2, r3}
    add     r0, sp, #304    @ 0x130
    bl      bf128_mul
    ldr     r1, [sp, #312]  @ 0x138
    ldr     r3, [sp, #36]   @ 0x24
    ldr     r2, [r6, #28]
    ldr     r0, [r6, #20]
    eor     lr, r1, r3
    ldr     r1, [sp, #40]   @ 0x28
    ldr     r3, [sp, #316]  @ 0x13c
    eors    r3, r1
    eors    r3, r2
    ldr     r2, [sp, #304]  @ 0x130
    str     r3, [r6, #28]
    eors    r2, r7
    ldr     r3, [sp, #308]  @ 0x134
    ldr     r7, [sp, #48]   @ 0x30
    eors    r3, r7
    eors    r3, r0
    str     r3, [r6, #20]
    ldr     r0, [sp, #24]
    ldr     r3, [sl, #40]   @ 0x28
    ldr     r7, [sp, #32]
    eors    r3, r0
    ldr     r0, [sl, #44]   @ 0x2c
    eors    r4, r0
    ldr     r0, [sl, #32]
    eors    r5, r0
    ldr     r0, [sl, #36]   @ 0x24
    eors    r7, r0
    ldr     r0, [r6, #24]
    eor     r1, lr, r0
    str     r1, [r6, #24]
    ldr     r1, [r6, #16]
    ldr     lr, [sp, #88]   @ 0x58
    eors    r2, r1
    str     r2, [r6, #16]
    ldr     r2, [sp, #404]  @ 0x194
    eors    r7, r2
    eor     r1, lr, r2
    str     r7, [r6, #36]   @ 0x24
    ldr     r2, [sp, #400]  @ 0x190
    ldr     r7, [sp, #84]   @ 0x54
    str     r1, [sp, #32]
    eor     r0, r7, r2
    eors    r5, r2
    ldr     r7, [sp, #80]   @ 0x50
    ldr     r2, [sp, #412]  @ 0x19c
    str     r5, [r6, #32]
    eors    r4, r2
    eors    r7, r2
    str     r4, [r6, #44]   @ 0x2c
    str     r0, [sp, #24]
    str     r7, [sp, #36]   @ 0x24
    ldr     r7, [sp, #408]  @ 0x198
    ldr     r2, [sp, #52]   @ 0x34
    ldr     r5, [sp, #28]
    eors    r3, r7
    add     lr, sp, #160    @ 0xa0
    str     r3, [r6, #40]   @ 0x28
    eors    r7, r2
    ldmia   lr, {r0, r1, r2, r3}
    stmia   r5, {r0, r1, r2, r3}
    add     lr, sp, #424    @ 0x1a8
    ldmia   lr, {r0, r1}
    stmia   sp, {r0, r1}
    add     lr, sp, #416    @ 0x1a0
    ldmia   lr, {r2, r3}
    add     r0, sp, #320    @ 0x140
    bl      bf128_mul
    ldrd    r4, r1, [r6, #36]       @ 0x24
    ldr     r0, [sp, #328]  @ 0x148
    ldr     r2, [r6, #44]   @ 0x2c
    ldr     r3, [r6, #32]
    ldr     ip, [sp, #28]
    eors    r1, r0
    str     r1, [r6, #40]   @ 0x28
    ldr     r1, [sp, #332]  @ 0x14c
    eors    r2, r1
    str     r2, [r6, #44]   @ 0x2c
    ldr     r2, [sp, #320]  @ 0x140
    eors    r3, r2
    str     r3, [r6, #32]
    ldr     r3, [sp, #324]  @ 0x144
    add     lr, sp, #176    @ 0xb0
    eors    r4, r3
    ldmia   lr, {r0, r1, r2, r3}
    str     r4, [r6, #36]   @ 0x24
    stmia   ip, {r0, r1, r2, r3}
    add     lr, sp, #440    @ 0x1b8
    ldmia   lr, {r0, r1}
    stmia   sp, {r0, r1}
    add     lr, sp, #432    @ 0x1b0
    ldmia   lr, {r2, r3}
    add     r0, sp, #336    @ 0x150
    bl      bf128_mul
    ldrd    r3, r0, [r6, #36]       @ 0x24
    ldr     r4, [sp, #344]  @ 0x158
    ldr     r1, [r6, #44]   @ 0x2c
    ldr     r2, [r6, #32]
    ldr     r5, [sp, #28]
    eors    r0, r4
    str     r0, [r6, #40]   @ 0x28
    ldr     r0, [sp, #348]  @ 0x15c
    eors    r1, r0
    str     r1, [r6, #44]   @ 0x2c
    ldr     r1, [sp, #336]  @ 0x150
    eors    r2, r1
    str     r2, [r6, #32]
    ldr     r2, [sp, #340]  @ 0x154
    add     lr, sp, #176    @ 0xb0
    eors    r3, r2
    str     r3, [r6, #36]   @ 0x24
    ldmia   lr, {r0, r1, r2, r3}
    stmia   r5, {r0, r1, r2, r3}
    add     ip, sp, #392    @ 0x188
    ldmia   ip, {r0, r1}
    stmia   sp, {r0, r1}
    add     ip, sp, #384    @ 0x180
    ldmia   ip, {r2, r3}
    add     r0, sp, #352    @ 0x160
    bl      bf128_mul
    ldr     r1, [sl, #60]   @ 0x3c
    ldr     r2, [sp, #36]   @ 0x24
    ldr     r4, [sp, #24]
    ldr     r3, [sl, #52]   @ 0x34
    ldr     r0, [sl, #56]   @ 0x38
    eors    r1, r2
    ldr     r2, [sl, #48]   @ 0x30
    eors    r2, r4
    ldr     r4, [sp, #32]
    eors    r3, r4
    ldr     r4, [sp, #360]  @ 0x168
    eors    r0, r7
    eors    r0, r4
    str     r0, [r6, #56]   @ 0x38
    ldr     r0, [sp, #364]  @ 0x16c
    eors    r1, r0
    str     r1, [r6, #60]   @ 0x3c
    ldr     r1, [sp, #352]  @ 0x160
    eors    r2, r1
    str     r2, [r6, #48]   @ 0x30
    ldr     r2, [sp, #356]  @ 0x164
    eors    r3, r2
    add     r2, sp, #160    @ 0xa0
    str     r3, [r6, #52]   @ 0x34
    ldmia   r2, {r0, r1, r2, r3}
    stmia   r5, {r0, r1, r2, r3}
    add     r2, sp, #440    @ 0x1b8
    ldmia   r2, {r0, r1}
    stmia   sp, {r0, r1}
    add     r2, sp, #432    @ 0x1b0
    ldmia   r2, {r2, r3}
    add     r0, sp, #368    @ 0x170
    bl      bf128_mul
    ldrd    r0, r1, [r6, #56]       @ 0x38
    ldr     r3, [sp, #376]  @ 0x178
    ldr     r2, [sp, #76]   @ 0x4c
    eors    r3, r0
    str     r3, [r6, #56]   @ 0x38
    ldr     r3, [sp, #380]  @ 0x17c
    adds    r2, #32
    eors    r3, r1
    str     r2, [sp, #76]   @ 0x4c
    ldrd    r2, r4, [r6, #48]       @ 0x30
    str     r3, [r6, #60]   @ 0x3c
    ldr     r3, [sp, #368]  @ 0x170
    eors    r3, r2
    str     r3, [r6, #48]   @ 0x30
    ldr     r3, [sp, #372]  @ 0x174
    eors    r3, r4
    str     r3, [r6, #52]   @ 0x34
    ldr     r3, [sp, #72]   @ 0x48
    cmp     r3, #128        @ 0x80
    add     r6, r6, #64     @ 0x40
    bne     aes_enc_forward_backward_128_l6
    add     sp, sp, #772    @ 0x304
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
aes_enc_forward_backward_128_l2:
    ldr     r8, [sp, #808]  @ 0x328
    str     r7, [sp, #24]
    add     r6, sp, #136    @ 0x88
    add     sl, sp, #224    @ 0xe0
    add     r9, sp, #240    @ 0xf0
aes_enc_forward_backward_128_l8:
    movs    r4, #0
aes_enc_forward_backward_128_l7:
    str     r7, [sp, #8]
    ldmia   r6, {r0, r1}
    stmia   sp, {r0, r1}
    add     r3, sp, #768    @ 0x300
    add     r3, r3, r4, lsl #4
    sub     r0, r3, #128    @ 0x80
    adds    r4, #1
    ldmia   r5, {r2, r3}
    bl      bf128_mul_bit
    cmp     r4, #8
    bne     aes_enc_forward_backward_128_l7
    mov     r1, fp
    mov     r0, sl
    bl      bf128_byte_combine
    ldr     r4, [sp, #24]
    ldr     r3, [sp, #824]  @ 0x338
    ldr     r1, [sp, #32]
    mov     r2, r4
    mov     r0, r9
    bl      bf128_byte_combine_vk_share
    ldr     r2, [sl, #8]
    ldr     r3, [r9, #8]
    eors    r3, r2
    str     r3, [r8, #8]
    ldr     r2, [sl, #12]
    ldr     r3, [r9, #12]
    eors    r3, r2
    str     r3, [r8, #12]
    ldr     r2, [sl]
    ldr     r3, [r9]
    eors    r3, r2
    str     r3, [r8]
    ldr     r2, [sl, #4]
    ldr     r3, [r9, #4]
    mov     r1, r4
    adds    r1, #8
    eors    r3, r2
    cmp     r1, #128        @ 0x80
    str     r3, [r8, #4]
    str     r1, [sp, #24]
    add     r8, r8, #16
    bne     aes_enc_forward_backward_128_l8
    b       aes_enc_forward_backward_128_l9
aes_enc_forward_backward_128_l1:
    ldr     r3, [sl, #4]
    str     r3, [sp, #56]   @ 0x38
    ldr     r3, [sl]
    str     r3, [sp, #52]   @ 0x34
    ldr     r3, [sl, #12]
    str     r3, [sp, #64]   @ 0x40
    ldr     r3, [sl, #8]
    str     r3, [sp, #60]   @ 0x3c
    ldr     r4, [sp, #816]  @ 0x330
    add     r3, sp, #96     @ 0x60
    mov     lr, r8
    add     r9, sp, #512    @ 0x200
    add     r5, sp, #152    @ 0x98
    str     r3, [sp, #44]   @ 0x2c
aes_enc_forward_backward_128_l13:
    movs    r3, #0
    str     r4, [sp, #28]
    str     r3, [sp, #24]
    str     lr, [sp, #40]   @ 0x28
    str     r4, [sp, #68]   @ 0x44
aes_enc_forward_backward_128_l12:
    ldr     r2, [sp, #24]
    ldr     r3, [sp, #40]   @ 0x28
    str     r9, [sp, #36]   @ 0x24
    subs    r3, r3, r2
    and     r3, r3, #3
    adds    r3, #36 @ 0x24
    add     r3, r2, r3, lsl #2
    bic     r8, r3, #3758096384     @ 0xe0000000
    ldr     r2, [sp, #48]   @ 0x30
    lsls    r3, r3, #3
    sub     r8, r8, #144    @ 0x90
    mov     r4, r9
    ldr     r9, [sp, #32]
    add     r8, r2
    add     r7, r3, #128    @ 0x80
    movs    r6, #0
aes_enc_forward_backward_128_l10:
    ldrb    r3, [r8]
    asrs    r3, r6
    and     r3, r3, #1
    str     r3, [sp, #8]
    ldmia   r5, {r0, r1}
    stmia   sp, {r0, r1}
    mov     r0, fp
    ldmia   sl, {r2, r3}
    bl      bf128_mul_bit
    adds    r1, r7, r6
    ldr     r2, [sp, #824]  @ 0x338
    mov     r0, r9
    bl      get_vk_128_share
    ldr     r3, [sp, #648]  @ 0x288
    ldr     r2, [r0, #8]
    ldr     r1, [sp, #652]  @ 0x28c
    eors    r2, r3
    ldr     r3, [r0, #12]
    eors    r3, r1
    ldrd    r1, ip, [r0]
    ldr     r0, [sp, #640]  @ 0x280
    eors    r0, r1
    ldr     r1, [sp, #644]  @ 0x284
    str     r0, [r4, #0]
    adds    r6, #1
    eor     r1, r1, ip
    cmp     r6, #8
    strd    r1, r2, [r4, #4]
    str     r3, [r4, #12]
    add     r4, r4, #16
    bne     aes_enc_forward_backward_128_l10
    ldr     r9, [sp, #36]   @ 0x24
    mov     r2, fp
    movs    r3, #7
aes_enc_forward_backward_128_l11:
    and     r4, r3, #7
    add     r6, r9, r4, lsl #4
    subs    r0, r3, #2
    lsls    r4, r4, #4
    and     r0, r0, #7
    ldr     r7, [r9, r4]
    ldr     r8, [r6, #4]
    add     r4, r9, r0, lsl #4
    lsls    r0, r0, #4
    subs    r1, r3, #5
    ldr     r0, [r9, r0]
    and     r1, r1, #7
    eor     ip, r0, r7
    add     r7, r9, r1, lsl #4
    lsls    r1, r1, #4
    ldr     r0, [r6, #8]
    ldr     r1, [r9, r1]
    eor     ip, ip, r1
    ldr     r1, [r4, #8]
    eor     lr, r0, r1
    ldr     r0, [r6, #12]
    ldr     r6, [r4, #12]
    eors    r0, r6
    ldr     r6, [r4, #4]
    ldr     r4, [r7, #8]
    eor     r1, lr, r4
    ldr     r4, [r7, #12]
    eors    r4, r0
    ldr     r0, [r7, #4]
    str     ip, [r2]
    adds    r3, #1
    eor     r6, r8, r6
    eors    r0, r6
    cmp     r3, #15
    strd    r0, r1, [r2, #4]
    str     r4, [r2, #12]
    add     r2, r2, #16
    bne     aes_enc_forward_backward_128_l11
    ldr     r2, [sp, #60]   @ 0x3c
    ldr     r3, [sp, #648]  @ 0x288
    ldr     r1, [sp, #64]   @ 0x40
    ldr     r0, [sp, #52]   @ 0x34
    ldr     r4, [sp, #56]   @ 0x38
    ldr     r6, [sp, #44]   @ 0x2c
    eors    r3, r2
    str     r3, [sp, #648]  @ 0x288
    ldr     r3, [sp, #652]  @ 0x28c
    eors    r3, r1
    str     r3, [sp, #652]  @ 0x28c
    ldr     r3, [sp, #640]  @ 0x280
    eors    r3, r0
    str     r3, [sp, #640]  @ 0x280
    ldr     r3, [sp, #644]  @ 0x284
    eors    r3, r4
    str     r3, [sp, #644]  @ 0x284
    ldr     r3, [sp, #680]  @ 0x2a8
    eors    r3, r2
    str     r3, [sp, #680]  @ 0x2a8
    ldr     r3, [sp, #684]  @ 0x2ac
    eors    r3, r1
    str     r3, [sp, #684]  @ 0x2ac
    ldr     r3, [sp, #672]  @ 0x2a0
    eors    r3, r0
    str     r3, [sp, #672]  @ 0x2a0
    ldr     r3, [sp, #676]  @ 0x2a4
    mov     r0, r6
    eors    r3, r4
    mov     r1, fp
    str     r3, [sp, #676]  @ 0x2a4
    bl      bf128_byte_combine
    ldmia   r6, {r0, r1, r2, r3}
    ldr     r4, [sp, #24]
    ldr     r6, [sp, #28]
    adds    r4, #1
    stmia   r6, {r0, r1, r2, r3}
    cmp     r4, #4
    add     r3, r6, #16
    str     r4, [sp, #24]
    str     r3, [sp, #28]
    bne     aes_enc_forward_backward_128_l12
    ldr     lr, [sp, #40]   @ 0x28
    ldr     r4, [sp, #68]   @ 0x44
    add     lr, lr, #1
    cmp     lr, #4
    add     r4, r4, #64     @ 0x40
    bne     aes_enc_forward_backward_128_l13
    add     sp, sp, #772    @ 0x304
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global bf8_inv_masked
bf8_inv_masked:
    stmdb   sp!, {r4, r5, r6, r7, r8, lr}
    sub     sp, #8
    mov     r6, r2
    mov     r5, r3
    mov     r8, r0
    mov     r7, r1
    bl      bf64_rand
    ldr     r3, =0x80808081
    adds    r4, r0, r1
    adc     r4, r4, #0
    umull   r3, ip, r3, r4
    mov     ip, ip, lsr #7
    rsb     ip, ip, ip, lsl #8
    sub     r4, r4, ip
    adds    r4, #1
    uxtb    r4, r4
    ldrb    r0, [r8]
    mov     r1, r4
    bl      bf8_mul
    ldrb    r1, [r7, #0]
    mov     r8, r0
    mov     r0, r4
    bl      bf8_mul
    eor     r0, r8, r0
    uxtb    r0, r0
    bl      bf8_inv
    movs    r3, #0
    mov     r7, r0
    movs    r1, #1
    add     r0, sp, #7
    strb    r3, [sp, #7]
    bl      rand_mask
    ldrb    r8, [sp, #7]
    eor     r0, r8, r7
    mov     r1, r4
    uxtb    r0, r0
    bl      bf8_mul
    mov     r1, r8
    strb    r6, [r6, #0]                // CLEAR
    strb    r0, [r6, #0]
    mov     r0, r4
    bl      bf8_mul
    strb    r5, [r5, #0]                // CLEAR
    strb    r0, [r5, #0]
    // CLEAR STACK
    mov     r1, #0
    str     r1, [sp]
    str     r1, [sp, #4]
    add     sp, #8
    ldmia   sp!, {r4, r5, r6, r7, r8, pc}
