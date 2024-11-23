.syntax unified
.thumb

.global aes_increment_iv
aes_increment_iv:
    push    {r4, lr}
    add     r4, r0, #16
    ldr     r0, [r4, #-4]!
    ldr     r1, [r4, #-4]!
    ldr     r2, [r4, #-4]!
    ldr     r3, [r4, #-4]!
    
    rev     r0, r0
    adds    r0, r0, #1
    rev     r0, r0
    
    rev     r1, r1    
    adcs    r1, r1, #0
    rev     r1, r1    
    
    rev     r2, r2
    adcs    r2, r2, #0
    rev     r2, r2

    rev     r3, r3
    adcs    r3, r3, #0
    rev     r3, r3
    
    str     r3, [r4], #4
    str     r2, [r4], #4
    str     r1, [r4], #4
    str     r0, [r4]
    pop     {r4, pc}

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

.global copy_first_round_key
copy_first_round_key:
    cbz     r2, cfrk_exit
    push    {r4, lr}
    mov     lr, #0
cfrk_loop:
    udiv    r4, lr, r3
    mls     ip, r3, r4, lr
    add     ip, ip, r4, lsl #3
    ldrb    r4, [r1, #0]
    strb    r4, [r0, ip, lsl #2]
    add     ip, r0, ip, lsl #2
    ldrb    r4, [r1, #1]
    strb    r4, [ip, #1]
    ldrb    r4, [r1, #2]
    strb    r4, [ip, #2]
    add     lr, lr, #1
    ldrb    r4, [r1, #3]
    strb    r4, [ip, #3]
    cmp     r2, lr
    add     r1, r1, #4
    bne     cfrk_loop
    mov     r0, #0          // CLEAR
    mov     r1, #0          // CLEAR
    mov     r2, #0          // CLEAR
    mov     r3, #0          // CLEAR
    mov     r4, #0          // CLEAR
    pop     {r4, pc}
cfrk_exit:
    bx      lr

.global setup_round_key_tmp
setup_round_key_tmp:
    push    {r4, r5, lr}
    ldr     r5, [sp, #12]
    add     ip, r2, #4294967295     @ 0xffffffff
    udiv    r4, ip, r5
    mls     ip, r5, r4, ip
    add     ip, ip, r4, lsl #3
    add     r1, r1, ip, lsl #2
    adds    r4, r0, #3
    add     ip, r0, #4294967295     @ 0xffffffff
srkt_loop:
    ldrb    lr, [r1], #1
    strb    lr, [ip, #1]!
    cmp     ip, r4
    bne     srkt_loop
    udiv    r1, r2, r3
    mls     r2, r3, r1, r2
    cbz     r2, srkt_rotword
    mov     lr, #0          // CLEAR
    mov     r0, #0          // CLEAR
    mov     r1, #0          // CLEAR
    mov     r2, #0          // CLEAR
    mov     r3, #0          // CLEAR
    mov     r4, #0          // CLEAR
    mov     r5, #0          // CLEAR
    pop     {r4, r5, pc}
srkt_rotword:
    mov     lr, #0          // CLEAR
    mov     r1, #0          // CLEAR
    mov     r2, #0          // CLEAR
    mov     r3, #0          // CLEAR
    mov     r4, #0          // CLEAR
    mov     r5, #0          // CLEAR
    ldmia   sp!, {r4, r5, lr}
    b       rot_word

.global finalize_round_key
finalize_round_key:
    push    {r4, lr}
    ldr     r4, [sp, #8]
    ldrb    lr, [r3]
    subs    r2, r1, r2
    udiv    ip, r2, r4
    mls     r2, r4, ip, r2
    add     r2, r2, ip, lsl #3
    ldrb    ip, [r0, r2, lsl #2]
    eor     ip, ip, lr
    udiv    lr, r1, r4
    mls     r1, r4, lr, r1
    add     r1, r1, lr, lsl #3
    add     r2, r0, r2, lsl #2
    strb    ip, [r0, r1, lsl #2]
    add     r0, r0, r1, lsl #2
    ldrb    r4, [r3, #1]
    ldrb    r1, [r2, #1]
    eors    r1, r4
    strb    r1, [r0, #1]
    ldrb    r1, [r2, #2]
    ldrb    r4, [r3, #2]
    eors    r1, r4
    strb    r1, [r0, #2]
    ldrb    r2, [r2, #3]
    ldrb    r3, [r3, #3]
    eors    r3, r2
    strb    r3, [r0, #3]
    mov     lr, #0          // CLEAR
    mov     r0, #0          // CLEAR
    mov     r1, #0          // CLEAR
    mov     r2, #0          // CLEAR
    mov     r3, #0          // CLEAR
    mov     r4, #0          // CLEAR
    pop     {r4, pc}

.global expand_128key_masked
expand_128key_masked:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, lr}
    sub     sp, #16
    mov     r8, r0
    ldr     r4, [sp, #48]   @ 0x30
    mov     r5, r2
    mov     r6, r3
    mov     r7, r1
    add.w   r9, r8, #480    @ 0x1e0
    bl      copy_first_round_key
    mla     sl, r4, r6, r6
    add.w   r1, r7, #32
    mov     r3, r6
    mov     r2, r5
    mov     r0, r9
    bl      copy_first_round_key
    cmp     r5, sl
    bcs.n   ekm_exit
    mov     r4, r5
    b.n     ekm_fwd
ekm_loop: 
    mov     r1, r4
    add     r3, sp, #8
    mov     r2, r5
    str     r6, [sp, #0]
    mov     r0, r8
    bl      finalize_round_key
    mov     r1, r4
    str     r6, [sp, #0]
    add     r3, sp, #12
    mov     r2, r5
    mov     r0, r9
    adds    r4, #1
    bl      finalize_round_key
    cmp     r4, sl
    beq.n   ekm_exit
ekm_fwd: 
    mov     r0, #0              //CLEAR
    str     r0, [sp, #8]        //CLEAR
    str     r0, [sp, #12]       //CLEAR
    mov     r3, r5
    mov     r2, r4
    mov     r1, r8
    str     r6, [sp, #0]
    add     r0, sp, #8
    bl      setup_round_key_tmp
    mov     r3, r5
    str     r6, [sp, #0]
    mov     r2, r4
    mov     r1, r9
    add     r0, sp, #12
    bl      setup_round_key_tmp
    udiv    r7, r4, r5
    mls     r3, r5, r7, r4
    cmp     r3, #0
    bne.n   ekm_loop
    add     r0, sp, #8
    subs    r7, #1
    bl      sub_words_masked
    uxtb    r0, r7
    bl      round_constants
    ldrb.w  r3, [sp, #8]
    eors    r0, r3
    strb.w  r0, [sp, #8]
    mov     r0, #0              //CLEAR
    b.n     ekm_loop
ekm_exit:
    // CLEAR STACK
    mov     r0, #0
    str     r0, [sp]
    str     r0, [sp, #4]
    str     r0, [sp, #8]
    str     r0, [sp, #12]
    mov     r1, #0
    mov     r2, #0
    mov     r3, #0
    mov     r4, #0
    mov     r5, #0
    mov     r6, #0
    mov     r7, #0
    mov     r8, #0
    mov     r9, #0
    mov     sl, #0
    add     sp, #16
    ldmia.w sp!, {r4, r5, r6, r7, r8, r9, sl, pc}


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