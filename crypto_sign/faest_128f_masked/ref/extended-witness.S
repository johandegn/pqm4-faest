.syntax unified
.thumb
/*
.global bf8_mul
// r0 = lhs (this is a share value)
// r1 = rhs (this is a share value)
bf8_mul:
    mov.w   ip, #0
    orr     ip, ip, ip                  // Clear (pipeline)
    push    {lr}                        // Clear
    push    {lr}                        // Clear
    pop     {lr}                        // ^^^^^
    pop     {lr}                        // ^^^^^

    push    {r4, lr}
    eor     r3, ip, ip                  // Clear
    eor     r4, ip, ip                  // Clear
    mov     lr, ip
mul_label:
    eor     r2, ip, ip                  // Clear
    asr.w   r2, r1, ip
    orr     ip, ip, ip                  // Clear (pipeline)
    sbfx    r3, r0, #7, #1
    sxtb    r4, r0
    orr     ip, ip, ip                  // Clear (pipeline)
    sbfx    r2, r2, #0, #1
    // pipeline leak???
    and.w   r3, r3, #27
    add.w   ip, ip, #1
    ands    r2, r0
    eor.w   r3, r3, r4, lsl #1
    cmp.w   ip, #7
    eor.w   lr, lr, r2
    uxtb    r0, r3
    bne.n   mul_label
    sbfx    r1, r1, #7, #1
    orr     ip, ip, ip                  // Clear (pipeline)
    ands    r1, r0
    orr     ip, ip, ip                  // Clear (pipeline)
    eor.w   r0, lr, r1
    uxtb    r0, r0
    eor     r1, ip, ip                  // Clear
    eor     r2, ip, ip                  // Clear
    eor     r3, ip, ip                  // Clear
    eor     r4, ip, ip                  // Clear
    pop     {r4, pc}
    nop
*/
.global bf8_square_masked
bf8_square_masked:
    push    {r3, r4, r5, lr}
    mov     r4, r1
    ldrb    r1, [r0, #0]
    mov     r5, r0
    mov     r0, r1
    bl      bf8_mul
    strb    r0, [r4, #0]
    eor     r1, r4, r4                  // Clear
    ldrb    r1, [r5, #1]
    eor     r0, r4, r4                  // Clear
    mov     r0, r1
    bl      bf8_mul
    strb    r0, [r4, #1]
    pop     {r3, r4, r5, pc}

.global bf8_mul_masked 
bf8_mul_masked:
    push    {r4, r5, r6, r7, r8, lr}
    mov     r8, r0
    sub     sp, #8
    mov     r5, r1
    add.w   r0, sp, #7
    movs    r1, #1
    mov     r6, r2
    bl      rand_mask
    ldrb.w  r7, [r8]
    ldrb    r1, [r5, #0]
    mov     r0, r7
    bl      bf8_mul
    eor     r1, r5, r5                  // Clear
    ldrb    r1, [r5, #1]
    mov     r4, r0
    mov     r0, r7
    ldrb.w  r7, [sp, #7]
    eors    r4, r7
    bl      bf8_mul
    eors    r4, r0
    strb    r4, [r6, #0]
    eor     r4, r5, r5                  // Clear
    ldrb    r1, [r5, #0]
    ldrb.w  r8, [r8, #1]
    mov     r0, r8
    bl      bf8_mul
    eor     r1, r5, r5                  // Clear
    ldrb    r1, [r5, #1]
    eor     r7, r0
    mov     r0, r8
    bl      bf8_mul
    eors    r7, r0
    strb    r7, [r6, #1]
    add     sp, #8
    pop     {r4, r5, r6, r7, r8, pc}

.global bf8_inv_masked
bf8_inv_masked:
    push    {r1}
    push    {r1}
    push    {r1}
    pop     {r1}
    pop     {r1}
    pop     {r1}
    push    {r4, r5, lr}
    eor     r2, r1, r1                  // Clear
    eor     r3, r1, r1                  // Clear
    eor     r4, r1, r1                  // Clear
    eor     r5, r1, r1                  // Clear
    sub     sp, #44 @ 0x2c
    movs    r3, #0
    mov     r5, r0
    mov     r4, r1
    mov     r1, sp
    strh.w  r3, [sp]
    strh.w  r3, [sp, #4]
    strh.w  r3, [sp, #8]
    strh.w  r3, [sp, #12]
    strh.w  r3, [sp, #16]
    strh.w  r3, [sp, #20]
    strh.w  r3, [sp, #24]
    strh.w  r3, [sp, #28]
    strh.w  r3, [sp, #32]
    strh.w  r3, [sp, #36]   @ 0x24
    bl      bf8_square_masked
    add     r2, sp, #4
    mov     r0, r5
    mov     r1, sp
    bl      bf8_mul_masked
    add     r2, sp, #8
    add     r0, sp, #4
    mov     r1, sp
    bl      bf8_mul_masked
    add     r2, sp, #12
    add     r0, sp, #8
    mov     r1, sp
    bl      bf8_mul_masked
    add     r1, sp, #16
    add     r0, sp, #12
    bl      bf8_square_masked
    add     r0, sp, #16
    add     r1, sp, #20
    bl      bf8_square_masked
    add     r0, sp, #20
    add     r1, sp, #24
    bl      bf8_square_masked
    add     r2, sp, #28
    add     r1, sp, #12
    add     r0, sp, #24
    bl      bf8_mul_masked
    add     r1, sp, #32
    add     r0, sp, #28
    bl      bf8_square_masked
    add     r1, sp, #36     @ 0x24
    add     r0, sp, #32
    bl      bf8_square_masked
    mov     r2, r4
    mov     r1, sp
    add     r0, sp, #36     @ 0x24
    bl      bf8_mul_masked
    add     sp, #44 @ 0x2c
    pop     {r4, r5, pc}
    nop

.global compute_sbox_masked
// r0 = in_share_ptr
// r1 = out_share_ptr
compute_sbox_masked:
    push    {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #20
    mov     r7, r1
    movs    r3, #0
    add     r1, sp, #8
    str     r7, [sp, #4]
    strh.w  r3, [sp, #8]
    bl      bf8_inv_masked
    add     r6, sp, #8

sbox_label:
    eors    r5, r6, r6                  // Clear
    ldrb    r5, [r6, #0]
    
    eors    r0, r6, r6                  // Clear
    eors    r2, r6, r6                  // Clear
    eors    r3, r6, r6                  // Clear
    eors    r4, r6, r6                  // Clear
    eors    r7, r6, r6                  // Clear
    eors    r8, r6, r6                  // Clear
    eors    r9, r6, r6                  // Clear
    eors    fp, r6, r6                  // Clear
    eors    sl, r6, r6                  // Clear
    str     r0, [sp]                    // Clear
    
    and.w   r0, r5, #241    @ 0xf1
    bl      __paritysi2
    mov     r7, r0
    and.w   r0, r5, #227    @ 0xe3
    bl      __paritysi2
    str     r0, [sp, #0]
    and.w   r0, r5, #199    @ 0xc7
    bl      __paritysi2
    mov     r4, r0
    and.w   r0, r5, #143    @ 0x8f
    bl      __paritysi2
    mov     r8, r0
    and.w   r0, r5, #31
    bl      __paritysi2
    mov     fp, r0
    and.w   r0, r5, #62     @ 0x3e
    bl      __paritysi2
    mov     sl, r0
    and.w   r0, r5, #124    @ 0x7c
    bl      __paritysi2
    mov     r9, r0
    and.w   r0, r5, #248    @ 0xf8
    bl      __paritysi2
    ldr     r2, [sp, #0]
    lsls    r3, r4, #2
    eor.w   r3, r3, r2, lsl #1
    eors    r3, r7
    eor.w   r3, r3, r8, lsl #3
    eor.w   r3, r3, fp, lsl #4
    eor.w   r3, r3, sl, lsl #5
    eor.w   r3, r3, r9, lsl #6
    eor.w   r3, r3, r0, lsl #7
    strb.w  r3, [r6], #1
    add.w   r3, sp, #10
    cmp     r3, r6
    bne.n   sbox_label

    eors    r3, r6, r6                  // Clear
    ldrb.w  r3, [sp, #9]
    eors    r1, r6, r6                  // Clear
    ldr     r1, [sp, #4]
    eors    r2, r6, r6                  // Clear
    ldrb.w  r2, [sp, #8]
    strb    r2, [r1, #0]
    orr     r6, r6, r6                  // Clear (pipeline)
    strb    r1, [sp, #4]
    eor.w   r3, r3, #99     @ 0x63
    strb    r3, [r1, #1]
    
    eors    r0, r6, r6                  // Clear
    eors    r1, r6, r6                  // Clear
    eors    r2, r6, r6                  // Clear
    eors    r3, r6, r6                  // Clear
    eors    r4, r6, r6                  // Clear
    eors    r5, r6, r6                  // Clear
    eors    r7, r6, r6                  // Clear
    eors    r8, r6, r6                  // Clear
    eors    r9, r6, r6                  // Clear
    eors    fp, r6, r6                  // Clear
    eors    sl, r6, r6                  // Clear
    add     sp, #20
    pop     {r4, r5, r6, r7, r8, r9, sl, fp, pc}
    nop

.global sub_bytes_masked
sub_bytes_masked:
    stmdb   sp!, {r4, r5, r6, r7, r8, lr}
    eors    r2, lr, lr
    eors    r3, lr, lr
    eors    r4, lr, lr
    eors    r5, lr, lr
    eors    r6, lr, lr
    eors    r7, lr, lr
    eors    r8, lr, lr

    movs    r7, #0
    sub     sp, #8
    mov     r8, r1
    adds    r5, r0, #4
    mov     r6, r7
loop_2:
    subs    r4, r5, #4
loop_1:
    eors    r2, lr, lr                  // Clear
    strb.w  r2, [sp]
    eors    r3, lr, lr                  // Clear
    strb.w  r3, [sp, #1]
    
    ldrb.w  r2, [r4], #1
    strb.w  r2, [sp]
    ldrb    r3, [r4, #31]
    strb.w  r3, [sp, #1]
    add     r1, sp, #4
    mov     r0, sp
    strh.w  r6, [sp, #4]
    eors    r2, lr, lr                  // Clear
    eors    r3, lr, lr                  // Clear
    bl      compute_sbox_masked
    eors    r2, lr, lr                  // Clear
    eors    r3, lr, lr                  // Clear
    ldrb.w  r2, [sp, #4]
    ldrb.w  r3, [sp, #5]
    strb.w  r2, [r4, #-1]
    cmp     r4, r5
    strb    r3, [r4, #31]
    bne.n   loop_1
    adds    r7, #1
    cmp     r8, r7
    add.w   r5, r5, #4
    bne.n   loop_2
    add     sp, #8

    eors    r0, lr, lr
    eors    r1, lr, lr
    eors    r2, lr, lr
    eors    r3, lr, lr
    eors    r4, lr, lr
    eors    r5, lr, lr
    eors    r6, lr, lr
    eors    r7, lr, lr
    eors    r8, lr, lr
    ldmia.w sp!, {r4, r5, r6, r7, r8, pc}

.global sub_words_masked
sub_words_masked:
    push    {r4, r5, r6, r7, lr}
    adds    r5, r0, #3
    sub     sp, #12
    subs    r4, r0, #1
    mov     r7, r5
    movs    r6, #0
sub_words_label:
    ldrb.w  r2, [r4, #1]!
    strb.w  r2, [sp]
    push    {lr}                        // Clear LSU pipeline
    pop     {lr}                        // ^^^^^
    ldrb.w  r3, [r5, #1]!
    add     r1, sp, #4
    mov     r0, sp
    strb.w  r3, [sp, #1]
    strh.w  r6, [sp, #4]
    bl      compute_sbox_masked
    ldrb.w  r2, [sp, #4]
    strb    r2, [r4, #0]
    push    {lr}                        // Clear LSU pipeline
    pop     {lr}                        // ^^^^^
    ldrb.w  r3, [sp, #5]
    cmp     r4, r7
    strb    r3, [r5, #0]
    bne.n   sub_words_label
    add     sp, #12
    pop     {r4, r5, r6, r7, pc}
