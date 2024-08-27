# Hey Emacs, this is a -*- makefile -*-
#----------------------------------------------------------------------------
#
# Makefile for ChipWhisperer SimpleSerial-AES Program
#
#----------------------------------------------------------------------------
# On command line:
#
# make all = Make software.
#
# make clean = Clean out built project files.
#
# make coff = Convert ELF to AVR COFF.
#
# make extcoff = Convert ELF to AVR Extended COFF.
#
# make program = Download the hex file to the device, using avrdude.
#                Please customize the avrdude settings below first!
#
# make debug = Start either simulavr or avarice as specified for debugging,
#              with avr-gdb or avr-insight as the front end for debugging.
#
# make filename.s = Just compile filename.c into the assembler code only.
#
# make filename.i = Create a preprocessed source file for use in submitting
#                   bug reports to the GCC project.
#
# To rebuild project do "make clean" then "make all".
#----------------------------------------------------------------------------

# Target file name (without extension).
# This is the base name of the compiled .hex file.
TARGET = simpleserial-faest-128f-masked

# List C source files here.
# Header files (.h) are automatically pulled in.
SRC += main.c

SRC += aes.c
SRC += compat.c
SRC += crypto_sign.c
SRC += faest_128f.c
SRC += faest_aes.c
SRC += faest.c
SRC += fields.c
SRC += instances.c
SRC += owf.c
SRC += random_oracle.c
SRC += universal_hashing.c
SRC += vbb.c
SRC += vc.c
SRC += vole.c
SRC += randomness.c
SRC += aes-publicinputs.c
SRC += randombytes.c
SRC += rng_common_v1.c

# Masked keccak sources
SRC += KeccakHash.c
SRC += KeccakHashtimes4.c
SRC += KeccakMPCHash.c
SRC += KeccakMPCSponge.c
SRC += KeccakP-1600-mpc.c
#SRC += KeccakP-1600-times4-on1.c
SRC += KeccakSponge.c
#SRC += KeccakSpongetimes4.c

#SRC += fips202.c
SRC += keccakf1600.c


ASRC += aes-encrypt.S
ASRC += extended-witness.S
ASRC += aes-keyschedule.S
ASRC += aes-publicinputs_.S
ASRC += keccakf1600.S
ASRC += KeccakP-1600-mpc_.S


# -----------------------------------------------------------------------------
#EXTRA_OPTS = 
#CFLAGS += -D$(EXTRA_OPTS)

ASFLAGS += -DSTM32F4

OPT = 2

ifeq ($(CRYPTO_TARGET),)
  ${info No CRYPTO_TARGET passed - defaulting to TINYAES128C}
  CRYPTO_TARGET = NONE
endif

ifeq ($(CRYPTO_OPTIONS),)
CRYPTO_OPTIONS = NONE
endif

${info Building for platform ${PLATFORM} with CRYPTO_TARGET=$(CRYPTO_TARGET)}

#Add simpleserial project to build
include ../simpleserial/Makefile.simpleserial

FIRMWAREPATH = ../.
include $(FIRMWAREPATH)/Makefile.inc

