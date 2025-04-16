.syntax unified
.thumb

/*
    R0: transpose_date_0
    R1: transpose_date_1
    R2: mod index + 25
 */
.global sca_transpose
sca_transpose:
    ror     r0, r0, r2
    ror     r1, r1, r2
    ldr     r2, =0xffffff0f

    ror     r0, r0, #7
    rrxs    r0, r0
    ror     r0, r0, #6
    rrx     r0, r0
    ror     r0, r0, #8
    rrxs    r0, r0
    ror     r0, r0, #24
    rrx     r0, r0
    ror     r0, r0, #14
    rrxs    r0, r0
    ror     r0, r0, #17
    rrx     r0, r0
    ror     r0, r0, #28

    ror     r1, r1, #7
    rrxs    r1, r1
    ror     r1, r1, #6
    rrx     r1, r1
    ror     r1, r1, #8
    rrxs    r1, r1
    ror     r1, r1, #24
    rrx     r1, r1
    ror     r1, r1, #14
    rrxs    r1, r1
    ror     r1, r1, #17
    rrx     r1, r1
    ror     r1, r1, #28

    and     r1, r1, r2
    eor     r0, r1, r0, lsl #4
    rev     r0, r0
    rbit    r0, r0

    bx      lr

/*
    R0: vole_cache
    R1: v_buf
    R2: idx
 */
.global transpose_vole_asm
transpose_vole_asm:
    push    {r4, r5, r6, r7, r8, r9, r10, lr}
    sub     sp, #8
    mov     r6, r0
    mov     r7, r1
    // R2 = idx % 8
    and     r5, r2, #0x07
    // Offset vole_cache
    add     r6, r6, r2, lsr #3

    mov     r8, #0
    mov     r9, #0

    mov     r0, sp
    mov     r1, #4
    bl      rand_mask
    ldr     r8, [sp]


transpose_vole_loop:
    // Load vole_cache
    mov     r0, #0
    ldrb    r1, [r6]
    ldrb    r2, [r6, #234]
    ldrb    r3, [r6, #468]
    ldrb    r4, [r6, #702]

    orr     r0, r0, r1, lsl #24
    orr     r0, r0, r2, lsl #16
    orr     r0, r0, r3, lsl #8
    orr     r0, r0, r4
    //add     r6, r6, #936

    mov     r1, #0
    ldrb    r2, [r6, #936]
    ldrb    r3, [r6, #1170]
    ldrb    r4, [r6, #1404]
    ldrb    r10, [r6, #1638]
    orr     r1, r1, r2, lsl #24
    orr     r1, r1, r3, lsl #16
    orr     r1, r1, r4, lsl #8
    orr     r1, r1, r10
    add     r6, r6, #1872

    add     r2, r5, #25
    bl      sca_transpose

    lsr     r8, r8, #8
    orr     r8, r8, r0, lsl #24
    add     r9, r9, #1

    cmp     r9, #4
    bne     transpose_vole_loop
    str     r8, [r7]
    mov     r8, #0
    
    add     sp, #8
    pop     {r4, r5, r6, r7, r8, r9, r10, pc}