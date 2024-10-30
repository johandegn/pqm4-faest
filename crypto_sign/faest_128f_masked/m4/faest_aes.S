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

.global aes_key_schedule_backward_128_vbb_vk_round_share
aes_key_schedule_backward_128_vbb_vk_round_share:
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
    clear_stack 236
    add     sp, #236        @ 0xec
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}

.global aes_enc_forward_backward_128_share
aes_enc_forward_backward_128_share:
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
    clear_stack 772
    add     sp, sp, #772    @ 0x304
    ldmia   sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}


.global get_vk_128
get_vk_128:
    push    {r4, r5, r6, lr}
    cmp     r1, #127        @ 0x7f
    sub     sp, #48 @ 0x30
    mov     r4, r1
    mov     r5, r0
    bls     vk_128_l1
    tst     r1, #96 @ 0x60
    beq     vk_128_l2
    subs    r1, #128        @ 0x80
    bl      get_vk_128
    ldmia   r0, {r0, r1, r2, r3}
    stmia   sp, {r0, r1, r2, r3}
    sub     r1, r4, #32
    mov     r0, r5
    bl      get_vk_128
    ldmia   r0, {r0, r1, r2, r3}
    add     r4, sp, #16
    stmia   r4, {r0, r1, r2, r3}
    ldrd    r4, r0, [sp, #8]
    eors    r2, r4
    str     r2, [sp, #40]   @ 0x28
    ldr     r2, [sp, #28]
    ldr     r6, [r5, #56]   @ 0x38
    ldrd    r1, r3, [sp]
    eors    r2, r0
    str     r2, [sp, #44]   @ 0x2c
    ldr     r2, [sp, #16]
    eors    r2, r1
    ldr     r1, [sp, #20]
    str     r2, [sp, #32]
    eors    r3, r1
    str     r3, [sp, #36]   @ 0x24
    add     r4, sp, #32
    ldmia   r4!, {r0, r1, r2, r3}
    str     r0, [r6, #0]
    str     r1, [r6, #4]
    str     r2, [r6, #8]
    str     r3, [r6, #12]
    ldr     r0, [r5, #56]   @ 0x38
    add     sp, #48 @ 0x30
    pop     {r4, r5, r6, pc}
vk_128_l2:
    ldr     r2, [r0, #32]
    ldr     r0, [r0, #60]   @ 0x3c
    ldrh    r2, [r2, #0]
    and     r3, r1, #127    @ 0x7f
    lsrs    r4, r1, #7
    adds    r3, #96 @ 0x60
    lsrs    r2, r2, #3
    add     r3, r3, r4, lsl #5
    mul     r3, r2, r3
    ldr     r2, [r5, #56]   @ 0x38
    adds    r1, r0, r3
vk_128_l3:
    ldr     r4, [r0, r3]
    ldr     r0, [r1, #4]
    ldr     r3, [r1, #8]
    ldr     r1, [r1, #12]
    str     r1, [r2, #12]
    str     r4, [r2, #0]
    str     r0, [r2, #4]
    str     r3, [r2, #8]
    ldr     r0, [r5, #56]   @ 0x38
    mov     r4, r0
    clear_stack 48
    mov     r0, r4
    add     sp, #48 @ 0x30
    pop     {r4, r5, r6, pc}
vk_128_l1:
    ldr     r3, [r0, #32]
    ldr     r0, [r0, #60]   @ 0x3c
    ldrh    r3, [r3, #0]
    ldr     r2, [r5, #56]   @ 0x38
    lsrs    r3, r3, #3
    mul     r3, r1, r3
    adds    r1, r0, r3
    b       vk_128_l3
