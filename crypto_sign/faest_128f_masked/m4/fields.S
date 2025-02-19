.syntax unified
.thumb

.macro clear_stack, n
    mov     r1, #0
    .set    offset, 0
    .rept   (\n/4)
    str     r1, [sp, #offset]
    .set    offset, offset+4
    .endr
.endm

.global bf128_mul_sclf
bf128_mul_sclf:
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
    clear_stack 300
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

    bl      bf128_mul_sclf
    ldmia   r5, {r0, r1, r2, r3}
    stmia   r4, {r0, r1, r2, r3}

    mov     r0, #0                        // CLEAR
    mov     r1, r0                        // CLEAR
    mov     r2, r0                        // CLEAR
    mov     r3, r0                        // CLEAR
    clear_stack 44
    add     sp, #44 @ 0x2c
    pop     {r4, r5, r6, r7, pc}

.global bf128_byte_combine_sclf
bf128_byte_combine_sclf:
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
    bl      bf128_mul_sclf
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
    clear_stack 68
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
    bl      bf128_mul_sclf
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
    clear_stack 76
    add     sp, #76 @ 0x4c
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}

.global bf128_mul_bit_sclf
bf128_mul_bit_sclf:
    sub     sp, #8
    push    {r4}
    add     r1, sp, #4
    stmia   r1, {r2, r3}
    ldrb    r1, [sp, #20]
    ldr     r4, [sp, #8]
    sbfx    ip, r1, #0, #1
    and     r2, r2, ip
    str     r2, [r0, #0]
    ldr     r2, [sp, #12]
    mov     r1, ip
    and     r2, ip, r2
    ands    r1, r4
    ldr     r4, [sp, #16]
    str     r2, [r0, #8]
    mov     r2, ip
    ands    r2, r4
    pop     {r4}
    str     r1, [r0, #4]
    str     r2, [r0, #12]
    // CLEAR STACK AND REGS
    mov     r0, #0
    mov     r1, r0
    mov     r2, r0
    mov     r3, r0
    strd    r0, r1, [sp]
    add     sp, #8
    bx      lr

.global bf128_byte_combine_bits_sclf
bf128_byte_combine_bits_sclf:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #60 @ 0x3c
    movs    r6, #0
    mov     fp, r1
    str     r0, [sp, #20]
    and     r9, r1, #1
    mov     r8, r6
    mov     r7, r6
    movs    r4, #1
    add     r5, sp, #24
    add     sl, sp, #32
 bf128_byte_combine_bits_l:
    subs    r1, r4, #1
    mov     r0, r5
    bl      bf128_alpha_wrapper
    asr     r3, fp, r4
    and     r3, r3, #1
    str     r3, [sp, #8]
    ldmia   sl, {r0, r1}
    stmia   sp, {r0, r1}
    add     r0, sp, #40     @ 0x28
    ldmia   r5, {r2, r3}
    bl      bf128_mul_bit
    ldrd    r2, r3, [sp, #40]       @ 0x28
    adds    r4, #1
    eor     r9, r9, r2
    eors    r6, r3
    ldr     r2, [sp, #48]   @ 0x30
    ldr     r3, [sp, #52]   @ 0x34
    cmp     r4, #8
    eor     r8, r8, r2
    eor     r7, r7, r3
    bne     bf128_byte_combine_bits_l
    ldr     r3, [sp, #20]
    mov     r0, r3
    strd    r9, r6, [r3]
    strd    r8, r7, [r3, #8]
    clear_stack 60
    add     sp, #60 @ 0x3c
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global bf128_sum_poly_vbb_share
bf128_sum_poly_vbb_share:
    stmdb   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
    sub     sp, #12
    mov     r4, r1
    mov     r8, r2
    str     r1, [sp, #0]
    mov     r2, r3
    mov     fp, r0
    add     r1, r8, #127    @ 0x7f
    mov     r0, r4
    mov     sl, r3
    bl      get_vole_aes_128_share
    ldrd    r4, r6, [r0]
    ldrd    r3, r5, [r0, #8]
    add     r7, r8, #126    @ 0x7e
    str     fp, [sp, #4]
bf128_sum_poly_vbb_share_loop:
    adds    r3, r3, r3
    mov     r1, r7
    ldr     r0, [sp, #0]
    mov     r2, sl
    adc     fp, r5, r5
    orr     r9, r3, r6, lsr #31
    bl      get_vole_aes_128_share
    ldr     r3, [r0, #0]
    ldr     r1, [r0, #4]
    adds    r4, r4, r4
    eor     r4, r4, r3
    mov     r5, r5, asr #31
    ldr     r3, [r0, #8]
    ldr     r0, [r0, #12]
    adcs    r6, r6
    and     r5, r5, #135    @ 0x87
    cmp     r8, r7
    eor     r4, r4, r5
    eor     r3, r9, r3
    eor     r6, r6, r1
    eor     r5, fp, r0
    add     r7, r7, #4294967295     @ 0xffffffff
    bne     bf128_sum_poly_vbb_share_loop
    ldr     fp, [sp, #4]
    mov     r0, fp
    strd    r4, r6, [fp]
    strd    r3, r5, [fp, #8]
    mov     r0, #0                      // CLEAR
    mov     r1, #0                      // CLEAR
    mov     r2, r1                      // CLEAR
    mov     r3, r1                      // CLEAR
    str     r1, [sp, #0]                // CLEAR
    str     r1, [sp, #4]                // CLEAR
    str     r1, [sp, #8]                // CLEAR
    add     sp, #12
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}

.global xor_u8_array_wrapper
xor_u8_array_wrapper:
    cbz     r3, xor_u8_array_wrapper_exit
    push    {lr}
    subs    r3, #1
    add     ip, r0, #4294967295     @ 0xffffffff
    subs    r1, #1
    subs    r2, #1
    add     r0, r3
xor_u8_array_wrapper_loop:
    ldrb    r3, [ip, #1]!
    ldrb    lr, [r1, #1]!
    cmp     ip, r0
    eor     r3, r3, lr
    strb    r3, [r2, #1]!
    bne     xor_u8_array_wrapper_loop
    mov     r0, #0                      // CLEAR
    mov     r1, #0                      // CLEAR
    mov     r2, #0                      // CLEAR
    mov     r3, #0                      // CLEAR
    mov     lr, #0                      // CLEAR
    ldr     pc, [sp], #4
xor_u8_array_wrapper_exit:
    bx      lr