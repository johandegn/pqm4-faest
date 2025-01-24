.syntax unified
.thumb

.global compute_sbox_masked
compute_sbox_masked:
    push    {r4, r5, r6, lr}
    mov     r5, r0
    mov     r4, r1
    bl      bf8_inv_masked

    // CLEAR temp registers and ldr-str pipeline
    ldr     r6, [sp, #8]
    mov     r0, #0
    mov     r1, #0 
    mov     r2, #0
    mov     r3, #0
    str     r6, [sp, #8]

    ldrb    r0, [r5, #0]
    bl      affine_incomplete
    strb    r0, [r5, #0]

    // CLEAR temp registers and ldr-str pipeline
    ldr     r6, [sp, #8]
    mov     r0, #0
    mov     r1, #0 
    mov     r2, #0
    mov     r3, #0
    str     r6, [sp, #8]

    ldrb    r0, [r4, #0]
    bl      affine
    strb    r0, [r4, #0]

    // CLEAR temp registers and ldr-str pipeline
    ldr     r6, [sp, #8]
    mov     r0, #0
    mov     r1, #0 
    mov     r2, #0
    mov     r3, #0
    str     r6, [sp, #8]
    mov     r4, #0
    mov     r5, #0
    mov     r6, #0

    pop     {r4, r5, r6, pc}

bf8_mul_asm:
    // Input arguments:
    // r0 = a, r1 = b
    
    // Setup r2 for the product
    // r3 for bitmask
    // ip for garbage result
    mov     r2, #0
    mov     r3, #1
    mov     ip, #0

    // Iterate over bits of b
    // bit 0
    ands    ip, r3, r1, lsr #0
    ite     eq
    eoreq   ip, r2, r0, lsl #0
    eorne   r2, r2, r0, lsl #0

    // bit 1
    ands    ip, r3, r1, lsr #1
    ite     eq
    eoreq   ip, r2, r0, lsl #1
    eorne   r2, r2, r0, lsl #1

    // bit 2
    ands    ip, r3, r1, lsr #2
    ite     eq
    eoreq   ip, r2, r0, lsl #2
    eorne   r2, r2, r0, lsl #2

    // bit 3
    ands    ip, r3, r1, lsr #3
    ite     eq
    eoreq   ip, r2, r0, lsl #3
    eorne   r2, r2, r0, lsl #3

    // bit 4
    ands    ip, r3, r1, lsr #4
    ite     eq
    eoreq   ip, r2, r0, lsl #4
    eorne   r2, r2, r0, lsl #4

    // bit 5
    ands    ip, r3, r1, lsr #5
    ite     eq
    eoreq   ip, r2, r0, lsl #5
    eorne   r2, r2, r0, lsl #5

    // bit 6
    ands    ip, r3, r1, lsr #6
    ite     eq
    eoreq   ip, r2, r0, lsl #6
    eorne   r2, r2, r0, lsl #6

    // bit 7
    ands    ip, r3, r1, lsr #7
    ite     eq
    eoreq   ip, r2, r0, lsl #7
    eorne   r2, r2, r0, lsl #7

    mov     r0, #0                  // CLEAR
    mov     r1, #0                  // CLEAR
    mov     ip, #0                  // CLEAR
    mov     r3, #0                  // CLEAR

    // Field reduction
    mov     r0, r2
    mov     r1, r2, lsr #8
    eor     r0, r0, r1
    eor     r0, r0, r1, lsl #1
    eor     r0, r0, r1, lsl #3
    eor     r0, r0, r1, lsl #4

    eor     r3, r3, r1, lsr #7
    eor     r3, r3, r1, lsr #5
    eor     r3, r3, r1, lsr #4
    
    eor     r0, r0, r3
    eor     r0, r0, r3, lsl #1
    eor     r0, r0, r3, lsl #3
    eor     r0, r0, r3, lsl #4

    mov     r1, #0                  // CLEAR
    mov     r2, #0                  // CLEAR
    mov     r3, #0                  // CLEAR
    bx      lr

bf8_shift_reduce:
    // Input arguments:
    // r0 = x
    // r1 = shift

    lsl     r2, r0, r1

    // Field reduction
    mov     r0, r2
    mov     r1, r2, lsr #8
    eor     r0, r0, r1
    eor     r0, r0, r1, lsl #1
    eor     r0, r0, r1, lsl #3
    eor     r0, r0, r1, lsl #4

    eor     r3, r3, r1, lsr #7
    eor     r3, r3, r1, lsr #5
    eor     r3, r3, r1, lsr #4
    
    eor     r0, r0, r3
    eor     r0, r0, r3, lsl #1
    eor     r0, r0, r3, lsl #3
    eor     r0, r0, r3, lsl #4

    uxtb    r0, r0
    bx      lr

bf8_inv_asm:
    // Computing x^(-1) = x^(255)
    // Doing it by following order
    // 1 - 7 - 15 - 30 - 60 - 120 - 127 - 254 - 255
    // Input arguments:
    // r0 = x
    push    {lr}

    mov     r1, #7
    bl      bf8_shift_reduce
    mov     r1, #8
    bl      bf8_shift_reduce
    mov     r1, r0
    bl      bf8_mul_asm
    mov     r1, r0
    bl      bf8_mul_asm
    mov     r1, r0
    bl      bf8_mul_asm
    mov     r1, #7
    bl      bf8_shift_reduce
    mov     r1, r0
    bl      bf8_mul_asm
    mov     r1, #1
    bl      bf8_shift_reduce

    pop     {pc}

.global bf8_inv_masked
bf8_inv_masked:
    push    {r4, r5, r6, r7, r8, r10, r11, lr}

    // ####
    // Initial cleaning
    mov     r2, #0
    mov     r3, #0
    mov     r4, #0
    mov     r5, #0
    mov     r6, #0
    mov     r7, #0
    mov     r8, #0
    mov     r10, #0
    mov     r11, #0
    mov     ip, #0
    sub     sp, #8
    str     r4, [sp, #0]
    str     r4, [sp, #4]
    // ####

    // Input arguments:
    // r0 = a_ptr, r1 = b_ptr
    mov     r5, r0
    mov     r6, r1

    // Setup random number r
    bl      bf64_rand
    ldr     r3, =0x80808081
    add     ip, r0, r1
    adc     ip, ip, #0
    umull   r2, r3, r3, ip
    lsr     r3, r3, #7
    rsb     r3, r3, r3, lsl #8
    sub     ip, ip, r3
    add     ip, #1
    uxtb    r4, ip

    // ####
    // Residue clean up
    mov     r0, #0
    mov     r1, #0
    mov     r2, #0
    mov     r3, #0
    mov     ip, #0
    uxtb    r3, ip
    // Note that the uxtb instruction can perform a rotation
    // This is different from other ALU instruction, so it might go
    // trough the pipeline differently than other, hence have a temporary
    // register that store the mask, and can cause a leak when uxtb is
    // called before the bf8_inv table lookup
    // ####

    // setup random number r1
    mov     r1, #1
    mov     r0, sp
    bl      rand_mask
    ldrb    r10, [sp]
    
    // ####
    // Residue clean up
    mov     r0, #0
    mov     r1, #0
    mov     r2, #0
    mov     r3, #0
    // ####

    // r * r1
    mov     r0, r4
    mov     r1, r10
    bl      bf8_mul_asm
    mov     r11, r0
    mov     r0, #0          // CLEAR
    // Note that the bf8_mul_asm, cleans all registers it touched
    // and does not use the stack, hence we assume we are clean after
    // the return, except for the return value in r0


    // a * r
    ldrb    r0, [r5, #0]
    mov     r1, r4
    bl      bf8_mul_asm
    mov     r7, r0
    mov     r0, #0          // CLEAR

    // ####
    // Memory pipeline cleaning
    ldrb    r0, [sp, #4]
    ldr     r1, [sp, #4]
    strb    r1, [sp, #4]
    str     r0, [sp, #4]
    // Note that the secret have been loaded, hence traversed the memory pipeline
    // in some way. To ensure value is left, the pipeline is filled with 0's and
    // the stack pointer.
    // ####
// 456
    // b * r
    ldrb    r0, [r6, #0]
    mov     r1, r4
    bl      bf8_mul_asm
    mov     r8, r0
    mov     r0, #0          // CLEAR
// 523
    // ####
    // Memory pipeline cleaning
    ldrb    r0, [sp, #4]
    ldr     r1, [sp, #4]
    strb    r1, [sp, #4]
    str     r0, [sp, #4]
    // ####

    // Inverse
    // Inv + r1
// 528
    eor     r0, r8, r7
    uxtb    r0, r0
    bl      bf8_inv
    eor     r0, r0, r10
// 540    
    // Note the the instruction would not require clearning next to eachother
    // as the values is masked by the product with r1.
    // However as the inverse is a table lookup, the memory pipeline must be cleaned

    // ####
    // Residue clean up
    mov     r1, #0
    mov     r2, #0
    mov     r3, #0
    // ####

    // ####
    // Memory pipeline cleaning
    ldrb    r2, [sp, #4]
    ldr     r1, [sp, #4]
    strb    r1, [sp, #4]
    str     r2, [sp, #4]
    // ####
// 548
    // (Inv + r1) * r
    mov     r1, r4
    bl      bf8_mul_asm
    strb    r0, [r5, #0]
    mov     r0, #0          // CLEAR

    // ####
    // Memory pipeline cleaning
    ldrb    r2, [sp, #4]
    ldr     r1, [sp, #4]
    strb    r1, [sp, #4]
    str     r2, [sp, #4]
    // ####

    strb    r11, [r6, #0]
    mov     r11, #0         // CLEAR

    // ####
    // Memory pipeline cleaning
    ldrb    r2, [sp, #4]
    ldr     r1, [sp, #4]
    strb    r1, [sp, #4]
    str     r2, [sp, #4]
    // ####

    // Exit cleaning
    mov     r0, #0
    mov     r1, #0
    mov     r2, #0
    mov     r3, #0
    mov     r4, #0
    mov     r5, #0
    mov     r6, #0
    mov     r7, #0
    mov     r8, #0
    mov     r10, #0
    mov     r11, #0
    mov     ip, #0
    str     r4, [sp, #0]
    str     r4, [sp, #4]

    add     sp, #8
    pop     {r4, r5, r6, r7, r8, r10, r11, pc}