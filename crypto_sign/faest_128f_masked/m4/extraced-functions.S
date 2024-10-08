.syntax unified
.thumb
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
    beq.n   aes_extend_witness_masked_output_early_return
    push    {r4, r5, lr}
    adds    r5, r0, r2
    add.w   ip, r0, #4294967295     @ 0xffffffff
    add.w   lr, r0, r2, lsl #1
    mov     r3, r5
aes_extend_witness_masked_output_loop: 
    ldr     r4, [r1, #0]
    subs    r2, r3, r5
    ldrb    r4, [r4, r2]
    strb.w  r4, [ip, #1]!
    strb    r1, [r3]                // CLEAR
    ldr     r4, [r1, #4]
    ldrb    r2, [r4, r2]
    strb.w  r2, [r3], #1
    cmp     r3, lr
    bne.n   aes_extend_witness_masked_output_loop
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
    nop




.global aes_encrypt_round_masked
aes_encrypt_round_masked:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #20
    ldrd    r5, sl, [sp, #56]       @ 0x38
    add.w   r9, r0, #32
    ldr.w   fp, [sp, #64]   @ 0x40
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
    add.w   r2, r6, #480    @ 0x1e0
    ldr     r3, [r7, #4]
    str     r5, [sp, #0]
    mov     r0, r9
    mov     r1, r4
    bl      aes_encrypt_round_masked_inner
    mov     r0, #7                              // CLEAR
    mov     r2, #7                              // CLEAR
    strd    r0, r2, [sp, #4]                    // CLEAR
    ldr     r3, [r5, #0]
    add.w   r3, r3, r4, lsl #2
    str     r3, [r5, #0]
    // CLEAR STACK
    mov     r1, #0
    str     r1, [sp]
    str     r1, [sp, #4]
    str     r1, [sp, #8]
    str     r1, [sp, #12]
    str     r1, [sp, #16]
    add     sp, #20
    ldmia.w sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}

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
    adc.w   r4, r4, #0
    umull   r3, ip, r3, r4
    mov.w   ip, ip, lsr #7
    rsb     ip, ip, ip, lsl #8
    sub.w   r4, r4, ip
    adds    r4, #1
    uxtb    r4, r4
    ldrb.w  r0, [r8]
    mov     r1, r4
    bl      bf8_mul
    ldrb    r1, [r7, #0]
    mov     r8, r0
    mov     r0, r4
    bl      bf8_mul
    eor.w   r0, r8, r0
    uxtb    r0, r0
    bl      bf8_inv
    movs    r3, #0
    mov     r7, r0
    movs    r1, #1
    add.w   r0, sp, #7
    strb.w  r3, [sp, #7]
    bl      rand_mask
    ldrb.w  r8, [sp, #7]
    eor.w   r0, r8, r7
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
    ldmia.w sp!, {r4, r5, r6, r7, r8, pc}


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
    add.w   r3, r3, ip, lsl #3
    ldrb.w  ip, [r2]
    strb.w  ip, [r6, r3, lsl #2]
    add.w   r3, r6, r3, lsl #2
    ldrb.w  ip, [r2, #1]
    strb.w  ip, [r3, #1]
    ldrb.w  ip, [r2, #2]
    strb.w  ip, [r3, #2]
    mov     ip, r0
    adds    r0, #1
    ldrb.w  lr, [r2, #3]
    strb.w  lr, [r3, #3]
    cmp     r7, r0
    add.w   r2, r2, #4
    bne.n   expand_128key_masked_l2
    movs    r2, #0
expand_128key_masked_l3:
    udiv    r0, r2, r4
    mls     r3, r4, r0, r2
    add.w   r3, r3, r0, lsl #3
    add.w   r3, r6, r3, lsl #2
    ldrb.w  r0, [r1, #32]
    strb.w  r0, [r3, #480]  @ 0x1e0
    ldrb.w  r0, [r1, #33]   @ 0x21
    strb.w  r0, [r3, #481]  @ 0x1e1
    ldrb.w  r0, [r1, #34]   @ 0x22
    strb.w  r0, [r3, #482]  @ 0x1e2
    ldrb.w  lr, [r1, #35]   @ 0x23
    strb.w  lr, [r3, #483]  @ 0x1e3
    cmp     ip, r2
    add.w   r1, r1, #4
    add.w   r2, r2, #1
    bne.n   expand_128key_masked_l3
expand_128key_masked_l1:
    mla     r9, r5, r4, r4
    cmp     r7, r9
    bcs.n   expand_128key_masked_l4_step
    mov     r5, r7
    add.w   r8, r6, #480    @ 0x1e0
    b.n     expand_128key_masked_l5
expand_128key_masked_l6:
    uxtb    r1, r3
    mov     r3, #0                          // CLEAR
    subs    r3, r5, r7
    udiv    r2, r3, r4
    mls     r3, r4, r2, r3
    add.w   r3, r3, r2, lsl #3
    ldrb.w  lr, [sp, #1]
    ldrb.w  ip, [sp, #2]
    udiv    r0, r5, r4
    movs    r2, #0
    bfi     r2, r1, #0, #8
    mov     r1, #0                          // CLEAR
    bfi     r2, lr, #8, #8
    mov     lr, #0                          // CLEAR
    ldrb.w  lr, [sp, #3]
    mls     r1, r4, r0, r5
    bfi     r2, ip, #16, #8
    mov     ip, #0                          // CLEAR
    ldr.w   ip, [r6, r3, lsl #2]
    add.w   r1, r1, r0, lsl #3
    bfi     r2, lr, #24, #8
    mov     lr, #0                          // CLEAR
    eor.w   r2, r2, ip
    str.w   r2, [r6, r1, lsl #2]
    nop                                     // CLEAR
    ldr     ip, [sp, #40]                   // CLEAR
    nop                                     // CLEAR
    mov     r2, #0                          // CLEAR
    nop                                     // CLEAR
    str     r6, [sp, #0]                    // CLEAR
    nop                                     // CLEAR
    ldr     r2, [sp, #4]
    ldr     ip, [sp, #40]                   // CLEAR
    ldr.w   r3, [r8, r3, lsl #2]
    ldr     ip, [sp, #40]                   // CLEAR
    adds    r5, #1
    eors    r3, r2
    cmp     r5, r9
    str.w   r3, [r8, r1, lsl #2]
    str     r6, [sp, #4]                    // CLEAR
    bne.n   expand_128key_masked_l5
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
    add.w   r3, r3, r2, lsl #3
    add.w   r2, r3, #120    @ 0x78
    ldr.w   r3, [r6, r3, lsl #2]
    ldr     r1, [sp, #40]                   // CLEAR
    str     r3, [sp, #0]
    str     r6, [sp, #4]                    // CLEAR
    nop                                     // CLEAR
    ldr.w   r1, [r6, r2, lsl #2]
    ldr     r2, [sp, #40]                   // CLEAR
    str     r1, [sp, #4]
    str     r2, [sp, #40]                   // CLEAR
    mls     r2, r7, sl, r5
    cmp     r2, #0
    bne.n   expand_128key_masked_l6
    mov     r0, sp
    bl      rot_word
    add     r0, sp, #4
    bl      rot_word
    mov     r0, sp
    bl      sub_words_masked
    add.w   r0, sl, #4294967295     @ 0xffffffff
    uxtb    r0, r0
    bl      round_constants
    ldrb.w  r3, [sp]
    eors    r3, r0
    b.n     expand_128key_masked_l6
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
    ldmia.w sp!, {r4, r5, r6, r7, r8, r9, sl, pc}


.global rot_word
rot_word:
    sub     sp, #4

    ldrb.w  ip, [r0, #1]
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
    ldmia.w r3, {r0, r1}
    stmia.w sp, {r0, r1}
    movs    r5, #0
    str     r5, [sp, #8]
    add     r3, sp, #40     @ 0x28
    ldmia   r3, {r2, r3}
    add     r0, sp, #56     @ 0x38
    lsls    r4, r4, #3
    bl      bf128_mul_bit
    add.w   r3, r4, #128    @ 0x80
    str     r3, [sp, #28]
    ldr     r3, [sp, #56]   @ 0x38
    str     r5, [sp, #20]
    eor.w   lr, r3, #1
    add     r3, sp, #56     @ 0x38
    ldmia   r3, {r0, r1, r2, r3}
    add     r5, sp, #104    @ 0x68
    stmia.w r5, {r0, r1, r2, r3}
    str.w   lr, [sp, #104]  @ 0x68
    ldmia.w r5, {r0, r1, r2, r3}
    add.w   r9, sp, #56     @ 0x38
    ldr.w   ip, [sp, #60]   @ 0x3c
    stmia.w r9, {r0, r1, r2, r3}
    strd    ip, lr, [sp, #32]
    add     r4, sp, #72     @ 0x48
    mov     sl, r6
aes_key_loop:
    ldr     r3, [sp, #28]
    ldr     r2, [sp, #20]
    mov     fp, r5
    add.w   r9, r3, r2
    mov     r6, r5
aes_key_l1:
    mov     r1, r9
    mov     r2, r8
    mov     r0, sl
    bl      get_vole_aes_128_share
    ldmia   r0, {r0, r1, r2, r3}
    stmia.w r4, {r0, r1, r2, r3}
    mov     r1, r9
    mov     r3, r8
    mov     r2, r4
    mov     r0, sl
    bl      add_vole_to_vk_cache_share
    ldmia.w r4, {r0, r1, r2, r3}
    stmia.w r6, {r0, r1, r2, r3}
    adds    r6, #16
    add     r3, sp, #232    @ 0xe8
    cmp     r6, r3
    add.w   r9, r9, #1
    bne.n   aes_key_l1
    ldrd    r2, r3, [sp, #20]
    adds    r6, r3, r2
aes_key_l2:
    mov     r1, r6
    mov     r2, r8
    mov     r0, sl
    bl      get_vk_128_share
    ldrd    r1, r3, [r0, #8]
    ldrd    ip, r2, [r0]
    ldr.w   r0, [fp, #8]
    eors    r1, r0
    str.w   r1, [fp, #8]
    ldr.w   r1, [fp, #12]
    eors    r3, r1
    str.w   r3, [fp, #12]
    ldr.w   r3, [fp]
    eor.w   r3, r3, ip
    str.w   r3, [fp]
    ldr.w   r3, [fp, #4]
    eors    r3, r2
    str.w   r3, [fp, #4]
    add.w   fp, fp, #16
    add     r3, sp, #232    @ 0xe8
    cmp     fp, r3
    add.w   r6, r6, #1
    bne.n   aes_key_l2
    mov     r2, r7
    movs    r3, #7
    mov     fp, r7
aes_key_l3:
    subs    r6, r3, #2
    and.w   r0, r3, #7
    and.w   r6, r6, #7
    add.w   ip, r5, r6, lsl #4
    add.w   r9, r5, r0, lsl #4
    lsls    r6, r6, #4
    lsls    r0, r0, #4
    ldr     r6, [r5, r6]
    ldr     r0, [r5, r0]
    ldr.w   r7, [r9, #8]
    subs    r1, r3, #5
    and.w   r1, r1, #7
    eors    r0, r6
    add.w   r6, r5, r1, lsl #4
    lsls    r1, r1, #4
    adds    r3, #1
    ldr     r1, [r5, r1]
    eors    r0, r1
    ldr.w   r1, [ip, #8]
    eor.w   lr, r7, r1
    ldr.w   r7, [ip, #12]
    ldr.w   r1, [r9, #12]
    ldr.w   r9, [r9, #4]
    eors    r1, r7
    ldr.w   r7, [ip, #4]
    eor.w   ip, r9, r7
    ldr     r7, [r6, #8]
    eor.w   lr, lr, r7
    ldr     r7, [r6, #12]
    eors    r7, r1
    ldr     r1, [r6, #4]
    str     r0, [r2, #0]
    eor.w   r1, ip, r1
    cmp     r3, #15
    strd    r1, lr, [r2, #4]
    str     r7, [r2, #12]
    add.w   r2, r2, #16
    bne.n   aes_key_l3
    ldr     r3, [sp, #36]   @ 0x24
    str     r3, [sp, #56]   @ 0x38
    movs    r6, #0
    ldr     r3, [sp, #32]
    str     r3, [sp, #60]   @ 0x3c
    str     r6, [sp, #8]
    add     r3, sp, #64     @ 0x40
    ldmia.w r3, {r0, r1}
    stmia.w sp, {r0, r1}
    add     r3, sp, #56     @ 0x38
    ldmia   r3, {r2, r3}
    add     r0, sp, #88     @ 0x58
    bl      bf128_mul_bit
    str     r6, [sp, #8]
    ldrd    r1, r2, [fp, #8]
    ldr     r0, [sp, #96]   @ 0x60
    eors    r1, r0
    str.w   r1, [fp, #8]
    ldr     r1, [sp, #100]  @ 0x64
    ldrd    r3, r6, [fp]
    eors    r2, r1
    str.w   r2, [fp, #12]
    ldr     r2, [sp, #88]   @ 0x58
    eors    r3, r2
    str.w   r3, [fp]
    ldr     r3, [sp, #92]   @ 0x5c
    eors    r6, r3
    add     r3, sp, #64     @ 0x40
    ldmia.w r3, {r0, r1}
    str.w   r6, [fp, #4]
    stmia.w sp, {r0, r1}
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
    str.w   r3, [fp, #40]   @ 0x28
    ldr     r3, [sp, #100]  @ 0x64
    ldrd    r2, r6, [fp, #32]
    eors    r3, r1
    str.w   r3, [fp, #44]   @ 0x2c
    ldr     r3, [sp, #88]   @ 0x58
    eors    r3, r2
    str.w   r3, [fp, #32]
    ldr     r3, [sp, #92]   @ 0x5c
    eors    r3, r6
    str.w   r3, [fp, #36]   @ 0x24
    ldr     r3, [sp, #20]
    mov     r7, fp
    cmp     r3, #32
    add.w   r7, r7, #128    @ 0x80
    bne.w   aes_key_loop
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
    ldmia.w sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global bf128_mul
bf128_mul:
    sub     sp, #8
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #300        @ 0x12c
    add     r1, sp, #336    @ 0x150
    stmia.w r1, {r2, r3}
    add     r3, sp, #352    @ 0x160
    mov     r8, r0
    ldmia   r3, {r0, r1, r2, r3}
    add     r6, sp, #56     @ 0x38
    stmia.w r6, {r0, r1, r2, r3}
    movs    r4, #0
    movs    r5, #0
    movs    r7, #0
    strd    r4, r5, [sp, #40]       @ 0x28
    strd    r4, r5, [sp, #48]       @ 0x30
    mov     ip, r7
    add     r5, sp, #40     @ 0x28
    add     r4, sp, #8
    str.w   r8, [sp, #4]
bf128_mul_l3:
    movs    r3, #1
    lsl.w   r3, r3, ip
    add.w   r2, r5, r3, lsl #4
    ldrd    r1, r7, [r2, #8]
    lsls    r3, r3, #4
    adds    r1, r1, r1
    ldr     r3, [r5, r3]
    ldr     r0, [r2, #4]
    mov.w   r6, r7, asr #31
    and.w   r6, r6, #135    @ 0x87
    adcs    r7, r7
    adds    r3, r3, r3
    eor.w   r6, r6, r3
    ldr     r3, [r5, #24]
    orr.w   r1, r1, r0, lsr #31
    eor.w   r8, r1, r3
    ldr     r3, [r5, #28]
    eor.w   r9, r7, r3
    ldr     r3, [r5, #16]
    mov.w   r2, #2
    eor.w   sl, r3, r6
    ldr     r3, [r5, #20]
    lsl.w   r2, r2, ip
    adcs    r0, r0
    add.w   lr, r2, #1
    eor.w   fp, r3, r0
    add.w   r3, r5, r2, lsl #4
    strd    r1, r7, [r3, #8]
    lsls    r2, r2, #4
    add.w   r1, r5, lr, lsl #4
    mov.w   lr, lr, lsl #4
    str     r6, [r5, r2]
    strd    r8, r9, [r1, #8]
    str     r0, [r3, #4]
    str.w   sl, [r5, lr]
    str.w   fp, [r1, #4]
    cmp.w   ip, #0
    beq.n   bf128_mul_l1
    add.w   lr, r4, r2
    adds    r2, #8
    mov     r1, r5
    add     r2, r5
bf128_mul_l2:
    ldrd    r6, r7, [r1, #40]       @ 0x28
    ldr     r0, [r2, #0]
    ldr.w   r9, [r2, #-4]
    eors    r6, r0
    ldr     r0, [r2, #4]
    eor.w   r8, r0, r7
    ldr     r0, [r1, #32]
    ldr.w   r7, [r2, #-8]
    eors    r0, r7
    ldr     r7, [r1, #36]   @ 0x24
    str     r0, [r3, #32]
    adds    r1, #16
    eor.w   r7, r7, r9
    cmp     lr, r1
    strd    r7, r6, [r3, #36]       @ 0x24
    str.w   r8, [r3, #44]   @ 0x2c
    add.w   r3, r3, #16
    bne.n   bf128_mul_l2
bf128_mul_l1:
    add.w   ip, ip, #1
    cmp.w   ip, #3
    bne.n   bf128_mul_l3
    ldr.w   r8, [sp, #4]
    movs    r2, #32
    movs    r1, #0
    mov     r0, r4
    movs    r6, #0
    bl      memset
    mov.w   lr, #60 @ 0x3c
    mov     ip, r6
    str.w   r8, [sp, #4]
bf128_mul_l6:
    ldr.w   r9, [r4]
    ldr     r2, [sp, #12]
    add     r0, sp, #336    @ 0x150
    mov     r1, r4
    rsb     r8, lr, #32
    sub.w   r7, lr, #32
bf128_mul_l4:
    ldrd    r3, sl, [r0]
    lsl.w   fp, sl, r8
    lsr.w   r3, r3, lr
    lsr.w   sl, sl, r7
    orr.w   r3, r3, fp
    orr.w   r3, r3, sl
    and.w   r3, r3, #15
    add.w   sl, r5, r3, lsl #4
    lsls    r3, r3, #4
    adds    r0, #8
    ldr     r3, [r5, r3]
    eor.w   r3, r3, r9
    str     r3, [r1, #0]
    ldr.w   r3, [sl, #4]
    eors    r3, r2
    str     r3, [r1, #4]
    ldr.w   r2, [r1, #8]!
    ldr.w   r3, [sl, #8]
    eor.w   r9, r2, r3
    ldr     r3, [r1, #4]
    ldr.w   r2, [sl, #12]
    str.w   r9, [r1]
    eors    r2, r3
    add     r3, sp, #352    @ 0x160
    cmp     r3, r0
    str     r2, [r1, #4]
    bne.n   bf128_mul_l4
    sub.w   lr, lr, #4
    cmn.w   lr, #4
    beq.n   bf128_mul_l5
    ldr     r1, [r4, #20]
    ldr     r2, [r4, #12]
    ldr.w   r8, [r4, #16]
    ldr     r3, [r4, #4]
    ldr     r7, [r4, #8]
    lsrs    r0, r1, #28
    mov.w   ip, ip, lsl #4
    lsls    r1, r1, #4
    orr.w   r1, r1, r8, lsr #28
    orr.w   ip, ip, r6, lsr #28
    orr.w   r6, r0, r6, lsl #4
    lsrs    r0, r2, #28
    orr.w   r0, r0, r8, lsl #4
    str     r1, [r4, #20]
    ldr     r1, [r4, #0]
    str     r0, [r4, #16]
    lsls    r2, r2, #4
    lsrs    r0, r3, #28
    lsls    r3, r3, #4
    orr.w   r3, r3, r1, lsr #28
    orr.w   r2, r2, r7, lsr #28
    orr.w   r0, r0, r7, lsl #4
    lsls    r1, r1, #4
    strd    r0, r2, [r4, #8]
    strd    r1, r3, [r4]
    b.n     bf128_mul_l6
bf128_mul_l5:
    ldr     r5, [sp, #24]
    ldr     r2, [sp, #16]
    ldr     r7, [sp, #20]
    ldr     r0, [sp, #28]
    ldr     r1, [sp, #12]
    ldr.w   r8, [sp, #4]
    mov.w   r3, ip, lsr #30
    eor.w   r3, r3, ip, lsr #31
    adds.w  lr, r6, r6
    eor.w   r3, r3, r5
    eor.w   r2, r2, r6
    eor.w   r3, r3, ip, lsr #25
    ldr     r5, [r4, #0]
    eor.w   r2, r2, lr
    eor.w   r4, ip, r7
    adc.w   lr, ip, ip
    eor.w   r2, r2, r6, lsl #2
    adds.w  r9, r3, r3
    mov.w   r7, ip, lsl #2
    mov.w   ip, ip, lsl #7
    orr.w   r7, r7, r6, lsr #30
    orr.w   ip, ip, r6, lsr #25
    eor.w   r2, r2, r6, lsl #7
    eor.w   r5, r5, r3
    adc.w   r6, r0, r0
    eors    r1, r0
    eor.w   r4, r4, lr
    mov.w   lr, r0, lsl #2
    eor.w   r5, r5, r9
    eor.w   r2, r2, r0, lsr #31
    eors    r1, r6
    orr.w   lr, lr, r3, lsr #30
    lsls    r6, r0, #7
    eors    r4, r7
    eor.w   r2, r2, r0, lsr #30
    eor.w   r5, r5, r3, lsl #2
    eor.w   r1, r1, lr
    orr.w   r6, r6, r3, lsr #25
    eor.w   r2, r2, r0, lsr #25
    eor.w   r4, r4, ip
    eor.w   r5, r5, r3, lsl #7
    eors    r1, r6
    mov     r0, r8
    strd    r5, r1, [r8]
    strd    r2, r4, [r8, #8]
// CLEAR STACK
    mov     r1, #300
    mov     r0, #0
bf128_mul_clear_loop:
    subs    r1, #4
    str     r0, [sp, r1]
    bne     bf128_mul_clear_loop
    ldr     r0, [sp]
    mov     r2, #0
    mov     r3, #0
// END OF CLEAR
    add     sp, #300        @ 0x12c
    ldmia.w sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    str     r0, [sp]                        // CLEAR
    str     r2, [sp, #4]                    // CLEAR
    add     sp, #8
    bx      lr


.global bf128_add_noinline
bf128_add_noinline:
    sub     sp, #8
    push    {r4}
    add     r1, sp, #4
    stmia   r1, {r2, r3}
    ldrd    r3, r1, [sp, #28]
    ldr     r4, [sp, #12]
    eors    r4, r3
    str     r4, [r0, #8]
    ldr     r4, [sp, #20]
    eors    r2, r4
    str     r2, [r0, #0]
    ldr     r2, [sp, #16]
    eors    r1, r2
    str     r1, [r0, #12]
    ldr     r2, [sp, #24]
    ldr     r1, [sp, #8]
    pop     {r4}
    eors    r2, r1
    str     r2, [r0, #4]
    add     r1, sp, #4                    // CLEAR
    add     sp, #8
    mov     r2, #0                        // CLEAR
    mov     r3, r2                        // CLEAR
    stmia   r1, {r2, r3}                  // CLEAR
    bx      lr