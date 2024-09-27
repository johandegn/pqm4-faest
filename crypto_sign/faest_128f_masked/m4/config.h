#define HAVE_RANDOMBYTES
#define PQCLEAN

#define XKCP_has_Sponge_Keccak
#define XKCP_has_FIPS202
#define XKCP_has_KeccakP1600
#define XKCP_has_KeccakP1600times4
#define XKCP_has_MaskedKeccakP1600

/* Configurations used in the paper. */
//#define CONFIG_MASK_NONE         // No Keccak masking.
//#define CONFIG_SNI_SECURE        // Refresh AND inputs/output and Keccak result, "All-SNI" in the paper.
//#define CONFIG_DOM_SECURE        // Randomize ANDs using Domain-Oriented Masking, "All-DOM" in the paper.
//#define CONFIG_IND_SECURE        // Mask all calls to Keccak using independent values, "All-IND" in the paper.
#define CONFIG_MASK_SEL            // Selective masking with independent values.
//#define CONFIG_SEL_HALF          // Selective half-masking with independent values.
//#define CONFIG_LEAK_TEST         // Reduced version for running leakage test.
//#define CONFIG_TEST_VECTORS      // Compatible with known-answer tests.

/* SHA-3 (Keccak) configuration. */

#ifdef CONFIG_MASK_NONE
#define KECCAK_MASK_NONE    // Mask no calls to SHA-3.
#endif

#ifdef CONFIG_SNI_SECURE
#define KECCAK_MASK_ALL     // Mask all calls to SHA-3.
#define KECCAK_SNI_SECURE   // Refresh AND inputs/output and Keccak result.
#endif

#ifdef CONFIG_DOM_SECURE
#define KECCAK_MASK_ALL     // Mask all calls to SHA-3.
#define KECCAK_DOM_SECURE   // Refresh AND operations.
#endif

#ifdef CONFIG_IND_SECURE
#define KECCAK_MASK_ALL     // Mask all calls to SHA-3.
#endif

#ifdef CONFIG_MASK_SEL
#undef KECCAK_MASK_ALL
#endif

#ifdef CONFIG_SEL_HALF
#undef KECCAK_MASK_ALL
#define KECCAK_MASK_HALF    // Mask first 12 rounds of SHA-3 only.
#endif

#ifdef CONFIG_LEAK_TEST
#define KECCAK_MASK_ALL       // Mask all calls to SHA-3.
#define KECCAK_MASK_PUBLIC    // Mask public values with random masks.
#define KECCAK_REDUCED_ROUNDS // Run only 4 rounds of Keccak.
#endif

#ifdef CONFIG_TEST_VECTORS
#define MASKED_IMPL_DETERMINISTIC // Do not randomize seeds and satisfy test vectors.
#endif
