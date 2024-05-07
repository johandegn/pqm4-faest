KEM_SCHEMES := \
	crypto_kem/bikel3/m4f \
	crypto_kem/kyber768-90s/m4fspeed \
	crypto_kem/kyber768-90s/m4fstack \
	crypto_kem/kyber1024/m4fspeed \
	crypto_kem/kyber1024/m4fstack \
	crypto_kem/bikel1/m4f \
	crypto_kem/kyber512/m4fspeed \
	crypto_kem/kyber512/m4fstack \
	crypto_kem/kyber768/m4fspeed \
	crypto_kem/kyber768/m4fstack \
	crypto_kem/kyber512-90s/m4fspeed \
	crypto_kem/kyber512-90s/m4fstack \
	crypto_kem/kyber1024-90s/m4fspeed \
	crypto_kem/kyber1024-90s/m4fstack \
	mupq/crypto_kem/bikel3/opt \
	mupq/crypto_kem/bikel1/opt \
	mupq/pqclean/crypto_kem/mceliece6688128f/clean \
	mupq/pqclean/crypto_kem/kyber768-90s/clean \
	mupq/pqclean/crypto_kem/kyber1024/clean \
	mupq/pqclean/crypto_kem/mceliece6960119/clean \
	mupq/pqclean/crypto_kem/mceliece348864/clean \
	mupq/pqclean/crypto_kem/kyber512/clean \
	mupq/pqclean/crypto_kem/kyber768/clean \
	mupq/pqclean/crypto_kem/mceliece8192128f/clean \
	mupq/pqclean/crypto_kem/mceliece348864f/clean \
	mupq/pqclean/crypto_kem/kyber512-90s/clean \
	mupq/pqclean/crypto_kem/hqc-rmrs-128/clean \
	mupq/pqclean/crypto_kem/mceliece6960119f/clean \
	mupq/pqclean/crypto_kem/kyber1024-90s/clean \
	mupq/pqclean/crypto_kem/hqc-rmrs-256/clean \
	mupq/pqclean/crypto_kem/hqc-rmrs-192/clean \
	mupq/pqclean/crypto_kem/mceliece460896/clean \
	mupq/pqclean/crypto_kem/mceliece6688128/clean \
	mupq/pqclean/crypto_kem/mceliece460896f/clean \
	mupq/pqclean/crypto_kem/mceliece8192128/clean

SIGN_SCHEMES := \
	crypto_sign/dilithium5/m4f \
	crypto_sign/mayo3/m4f \
	crypto_sign/faest_192f/ref \
	crypto_sign/haetae2/m4f \
	crypto_sign/dilithium3/m4f \
	crypto_sign/perk-192-fast-3/m4 \
	crypto_sign/perk-256-fast-5/m4 \
	crypto_sign/perk-256-short-5/m4 \
	crypto_sign/falcon-512-tree/m4-ct \
	crypto_sign/perk-128-short-5/m4 \
	crypto_sign/faest_128f/ref \
	crypto_sign/faest_128f/m4 \
	crypto_sign/ov-Ip-pkc/m4fspeed \
	crypto_sign/ov-Ip-pkc/m4fstack \
	crypto_sign/perk-192-short-5/m4 \
	crypto_sign/falcon-1024/m4-ct \
	crypto_sign/perk-192-short-3/m4 \
	crypto_sign/ov-Ip/m4f \
	crypto_sign/perk-128-short-3/m4 \
	crypto_sign/haetae3/m4f \
	crypto_sign/mayo1/m4f \
	crypto_sign/haetae5/m4f \
	crypto_sign/perk-192-fast-5/m4 \
	crypto_sign/perk-256-short-3/m4 \
	crypto_sign/ov-Ip-pkc-skc/m4fspeed \
	crypto_sign/ov-Ip-pkc-skc/m4fstack \
	crypto_sign/perk-128-fast-5/m4 \
	crypto_sign/falcon-512/m4-ct \
	crypto_sign/mayo2/m4f \
	crypto_sign/dilithium2/m4f \
	crypto_sign/faest_em_128f/ref \
	crypto_sign/perk-256-fast-3/m4 \
	crypto_sign/perk-128-fast-3/m4 \
	mupq/crypto_sign/mirith_Ib_short/ref \
	mupq/crypto_sign/mirith_Ia_fast/ref \
	mupq/crypto_sign/cross-sha2-r-sdpg-5-small/ref \
	mupq/crypto_sign/snova-37-8-16-4-esk/ref \
	mupq/crypto_sign/sphincs-a-shake-128s/ref \
	mupq/crypto_sign/meds167717/ref \
	mupq/crypto_sign/tuov_iii_pkc_skc/ref \
	mupq/crypto_sign/mirith_hypercube_Va_fast/ref \
	mupq/crypto_sign/mayo3/ref \
	mupq/crypto_sign/meds9923/ref \
	mupq/crypto_sign/sphincs-a-shake-256f/ref \
	mupq/crypto_sign/snova-61-33-16-2-ssk/ref \
	mupq/crypto_sign/tuov_iii/ref \
	mupq/crypto_sign/ascon-sign-192f-robust/ref \
	mupq/crypto_sign/snova-43-25-16-2-esk/ref \
	mupq/crypto_sign/sphincs-a-shake-192f/ref \
	mupq/crypto_sign/mirith_hypercube_Vb_short/ref \
	mupq/crypto_sign/snova-61-33-16-2-esk/ref \
	mupq/crypto_sign/snova-24-5-16-4-ssk/ref \
	mupq/crypto_sign/ascon-sign-128f-simple/ref \
	mupq/crypto_sign/haetae2/ref \
	mupq/crypto_sign/perk-192-fast-3/ref \
	mupq/crypto_sign/biscuit128f/ref \
	mupq/crypto_sign/biscuit256f/ref \
	mupq/crypto_sign/ascon-sign-128s-simple/ref \
	mupq/crypto_sign/cross-sha3-r-sdp-3-fast/ref \
	mupq/crypto_sign/cross-sha3-r-sdp-5-fast/ref \
	mupq/crypto_sign/mirith_hypercube_Va_short/ref \
	mupq/crypto_sign/sphincs-a-sha2-256f/ref \
	mupq/crypto_sign/perk-256-fast-5/ref \
	mupq/crypto_sign/cross-sha3-r-sdp-1-small/ref \
	mupq/crypto_sign/cross-sha2-r-sdp-1-fast/ref \
	mupq/crypto_sign/meds134180/ref \
	mupq/crypto_sign/aimer-l3-param2/ref \
	mupq/crypto_sign/mirith_hypercube_Ia_fast/ref \
	mupq/crypto_sign/mirith_hypercube_Ia_fast/opt \
	mupq/crypto_sign/snova-28-17-16-2-ssk/ref \
	mupq/crypto_sign/sphincs-a-sha2-128s/ref \
	mupq/crypto_sign/sphincs-a-sha2-128f/ref \
	mupq/crypto_sign/cross-sha3-r-sdpg-1-small/ref \
	mupq/crypto_sign/snova-49-11-16-3-esk/ref \
	mupq/crypto_sign/mirith_Vb_short/ref \
	mupq/crypto_sign/mqom_cat3_gf251_short/ref \
	mupq/crypto_sign/perk-256-short-5/ref \
	mupq/crypto_sign/meds41711/ref \
	mupq/crypto_sign/falcon-512-tree/opt-ct \
	mupq/crypto_sign/falcon-512-tree/opt-leaktime \
	mupq/crypto_sign/tuov_v_pkc_skc/ref \
	mupq/crypto_sign/perk-128-short-5/ref \
	mupq/crypto_sign/mqom_cat3_gf31_fast/ref \
	mupq/crypto_sign/cross-sha3-r-sdp-1-fast/ref \
	mupq/crypto_sign/biscuit128s/ref \
	mupq/crypto_sign/cross-sha3-r-sdp-5-small/ref \
	mupq/crypto_sign/cross-sha2-r-sdpg-5-fast/ref \
	mupq/crypto_sign/snova-37-8-16-4-ssk/ref \
	mupq/crypto_sign/mirith_hypercube_Vb_fast/ref \
	mupq/crypto_sign/aimer-l5-param1/ref \
	mupq/crypto_sign/snova-24-5-16-4-esk/ref \
	mupq/crypto_sign/snova-49-11-16-3-ssk/ref \
	mupq/crypto_sign/tuov_ip_pkc/ref \
	mupq/crypto_sign/ov-Ip-pkc/ref \
	mupq/crypto_sign/cross-sha3-r-sdpg-1-fast/ref \
	mupq/crypto_sign/perk-192-short-5/ref \
	mupq/crypto_sign/tuov_v/ref \
	mupq/crypto_sign/sphincs-a-shake-256s/ref \
	mupq/crypto_sign/mirith_Vb_fast/ref \
	mupq/crypto_sign/falcon-1024/opt-ct \
	mupq/crypto_sign/falcon-1024/opt-leaktime \
	mupq/crypto_sign/mirith_hypercube_Ia_shorter/ref \
	mupq/crypto_sign/perk-192-short-3/ref \
	mupq/crypto_sign/mirith_hypercube_IIIb_shorter/ref \
	mupq/crypto_sign/tuov_is_pkc_skc/ref \
	mupq/crypto_sign/ov-Ip/ref \
	mupq/crypto_sign/mqom_cat1_gf31_fast/ref \
	mupq/crypto_sign/perk-128-short-3/ref \
	mupq/crypto_sign/cross-sha2-r-sdp-5-small/ref \
	mupq/crypto_sign/snova-60-10-16-4-ssk/ref \
	mupq/crypto_sign/haetae3/ref \
	mupq/crypto_sign/ascon-sign-192f-simple/ref \
	mupq/crypto_sign/cross-sha2-r-sdpg-1-fast/ref \
	mupq/crypto_sign/mirith_IIIa_short/ref \
	mupq/crypto_sign/mqom_cat5_gf251_short/ref \
	mupq/crypto_sign/ascon-sign-128f-robust/ref \
	mupq/crypto_sign/aimer-l1-param2/ref \
	mupq/crypto_sign/mqom_cat1_gf251_fast/ref \
	mupq/crypto_sign/cross-sha3-r-sdp-3-small/ref \
	mupq/crypto_sign/hawk1024/ref \
	mupq/crypto_sign/mirith_hypercube_IIIa_short/ref \
	mupq/crypto_sign/mqom_cat3_gf31_short/ref \
	mupq/crypto_sign/cross-sha3-r-sdpg-3-fast/ref \
	mupq/crypto_sign/mirith_IIIb_fast/ref \
	mupq/crypto_sign/mayo1/ref \
	mupq/crypto_sign/cross-sha2-r-sdpg-3-small/ref \
	mupq/crypto_sign/ascon-sign-128s-robust/ref \
	mupq/crypto_sign/haetae5/ref \
	mupq/crypto_sign/cross-sha2-r-sdpg-3-fast/ref \
	mupq/crypto_sign/tuov_is/ref \
	mupq/crypto_sign/mirith_Va_short/ref \
	mupq/crypto_sign/snova-66-15-16-3-ssk/ref \
	mupq/crypto_sign/cross-sha2-r-sdpg-1-small/ref \
	mupq/crypto_sign/sphincs-a-shake-128f/ref \
	mupq/crypto_sign/aimer-l1-param3/ref \
	mupq/crypto_sign/mirith_hypercube_IIIa_shorter/ref \
	mupq/crypto_sign/mirith_hypercube_IIIa_fast/ref \
	mupq/crypto_sign/snova-43-25-16-2-ssk/ref \
	mupq/crypto_sign/snova-25-8-16-3-esk/ref \
	mupq/crypto_sign/mirith_hypercube_Ia_short/ref \
	mupq/crypto_sign/sphincs-a-sha2-192f/ref \
	mupq/crypto_sign/sphincs-a-sha2-192s/ref \
	mupq/crypto_sign/biscuit256s/ref \
	mupq/crypto_sign/cross-sha2-r-sdp-1-small/ref \
	mupq/crypto_sign/tuov_ip_pkc_skc/ref \
	mupq/crypto_sign/cross-sha2-r-sdp-3-small/ref \
	mupq/crypto_sign/perk-192-fast-5/ref \
	mupq/crypto_sign/cross-sha3-r-sdpg-5-fast/ref \
	mupq/crypto_sign/mirith_hypercube_Ib_fast/ref \
	mupq/crypto_sign/mirith_hypercube_Ib_fast/opt \
	mupq/crypto_sign/snova-66-15-16-3-esk/ref \
	mupq/crypto_sign/mqom_cat3_gf251_fast/ref \
	mupq/crypto_sign/biscuit192s/ref \
	mupq/crypto_sign/mirith_IIIb_short/ref \
	mupq/crypto_sign/ascon-sign-192s-simple/ref \
	mupq/crypto_sign/mirith_IIIa_fast/ref \
	mupq/crypto_sign/snova-60-10-16-4-esk/ref \
	mupq/crypto_sign/cross-sha3-r-sdpg-3-small/ref \
	mupq/crypto_sign/perk-256-short-3/ref \
	mupq/crypto_sign/snova-28-17-16-2-esk/ref \
	mupq/crypto_sign/mqom_cat5_gf251_fast/ref \
	mupq/crypto_sign/hawk512/ref \
	mupq/crypto_sign/cross-sha2-r-sdp-5-fast/ref \
	mupq/crypto_sign/sphincs-a-sha2-256s/ref \
	mupq/crypto_sign/tuov_v_pkc/ref \
	mupq/crypto_sign/meds55604/ref \
	mupq/crypto_sign/ov-Ip-pkc-skc/ref \
	mupq/crypto_sign/perk-128-fast-5/ref \
	mupq/crypto_sign/mirith_hypercube_Ib_short/ref \
	mupq/crypto_sign/aimer-l5-param2/ref \
	mupq/crypto_sign/mirith_Va_fast/ref \
	mupq/crypto_sign/falcon-512/opt-ct \
	mupq/crypto_sign/falcon-512/opt-leaktime \
	mupq/crypto_sign/tuov_ip/ref \
	mupq/crypto_sign/mayo2/ref \
	mupq/crypto_sign/mirith_Ib_fast/ref \
	mupq/crypto_sign/mqom_cat1_gf251_short/ref \
	mupq/crypto_sign/mirith_Ia_short/ref \
	mupq/crypto_sign/tuov_iii_pkc/ref \
	mupq/crypto_sign/mqom_cat1_gf31_short/ref \
	mupq/crypto_sign/aimer-l1-param1/ref \
	mupq/crypto_sign/mirith_hypercube_Vb_shorter/ref \
	mupq/crypto_sign/snova-25-8-16-3-ssk/ref \
	mupq/crypto_sign/cross-sha3-r-sdpg-5-small/ref \
	mupq/crypto_sign/falcon-1024-tree/opt-ct \
	mupq/crypto_sign/falcon-1024-tree/opt-leaktime \
	mupq/crypto_sign/cross-sha2-r-sdp-3-fast/ref \
	mupq/crypto_sign/perk-256-fast-3/ref \
	mupq/crypto_sign/biscuit192f/ref \
	mupq/crypto_sign/hawk256/ref \
	mupq/crypto_sign/mirith_hypercube_IIIb_fast/ref \
	mupq/crypto_sign/aimer-l3-param1/ref \
	mupq/crypto_sign/meds13220/ref \
	mupq/crypto_sign/sphincs-a-shake-192s/ref \
	mupq/crypto_sign/ascon-sign-192s-robust/ref \
	mupq/crypto_sign/perk-128-fast-3/ref \
	mupq/crypto_sign/mirith_hypercube_Ib_shorter/ref \
	mupq/crypto_sign/mirith_hypercube_Va_shorter/ref \
	mupq/crypto_sign/tuov_is_pkc/ref \
	mupq/crypto_sign/mirith_hypercube_IIIb_short/ref \
	mupq/pqclean/crypto_sign/sphincs-sha256-192s-simple/clean \
	mupq/pqclean/crypto_sign/dilithium5/clean \
	mupq/pqclean/crypto_sign/dilithium5aes/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-192s-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-256f-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-192s-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-256f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-128s-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-128s-simple/clean \
	mupq/pqclean/crypto_sign/dilithium3/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-192f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-256s-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-192f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-256f-robust/clean \
	mupq/pqclean/crypto_sign/dilithium3aes/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-256s-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-128f-robust/clean \
	mupq/pqclean/crypto_sign/falcon-1024/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-192f-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-256s-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-128f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-256s-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-256f-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-128s-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-256f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-192s-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-192f-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-256f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-128s-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-192f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-128f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-192f-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-256s-robust/clean \
	mupq/pqclean/crypto_sign/falcon-512/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-192s-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-192s-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-shake256-128f-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-128f-robust/clean \
	mupq/pqclean/crypto_sign/dilithium2/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-128f-robust/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-128s-simple/clean \
	mupq/pqclean/crypto_sign/sphincs-haraka-256s-robust/clean \
	mupq/pqclean/crypto_sign/dilithium2aes/clean \
	mupq/pqclean/crypto_sign/sphincs-sha256-128s-robust/clean