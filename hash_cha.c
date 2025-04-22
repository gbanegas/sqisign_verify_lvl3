//
// Created by gustavo on 18/04/25.
//
#include "fp.h"
#include "fp2.h"
#include "verification.h"
#include "encoded_sizes.h"
#include "structures.h"
#include "ec_params.h"
#include "fips202.h"
#include "mp.h"



void
ec_j_inv(fp2_t *j_inv, const ec_curve_t *curve)
{ // j-invariant computation for Montgommery coefficient A2=(A+2C:4C)
    fp2_t t0, t1;

    fp2_sqr(&t1, &curve->C);
    fp2_sqr(j_inv, &curve->A);
    fp2_add(&t0, &t1, &t1);
    fp2_sub(&t0, j_inv, &t0);
    fp2_sub(&t0, &t0, &t1);
    fp2_sub(j_inv, &t0, &t1);
    fp2_sqr(&t1, &t1);
    fp2_mul(j_inv, j_inv, &t1);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&t0, &t0, &t0);
    fp2_sqr(&t1, &t0);
    fp2_mul(&t0, &t0, &t1);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&t0, &t0, &t0);
    fp2_inv(j_inv);
    fp2_mul(j_inv, &t0, j_inv);
}

// compute the challenge as the hash of the message and the commitment curve and public key
void
hash_to_challenge(scalar_t *scalar,
                  const public_key_t *pk,
                  const ec_curve_t *com_curve,
                  const unsigned char *message,
                  size_t length)
{
    unsigned char buf[2 * FP2_ENCODED_BYTES];
    {
        fp2_t j1, j2;
        ec_j_inv(&j1, &pk->curve);
        ec_j_inv(&j2, com_curve);
        fp2_encode(buf, &j1);
        fp2_encode(buf + FP2_ENCODED_BYTES, &j2);
    }

    {
        // The type scalar_t represents an element of GF(p), which is about
        // 2*lambda bits, where lambda = 128, 192 or 256, according to the
        // security level. Thus, the variable scalar should have enough memory
        // for the values produced by SHAKE256 in the intermediate iterations.

        shake256incctx ctx;

        size_t hash_bytes = ((2 * SECURITY_BITS) + 7) / 8;
        size_t limbs = (hash_bytes + sizeof(digit_t) - 1) / sizeof(digit_t);
        size_t bits = (2 * SECURITY_BITS) % RADIX;
        digit_t mask = ((digit_t)-1) >> ((RADIX - bits) % RADIX);
#ifdef TARGET_BIG_ENDIAN
        mask = BSWAP_DIGIT(mask);
#endif

        shake256_inc_init(&ctx);
        shake256_inc_absorb(&ctx, buf, 2 * FP2_ENCODED_BYTES);
        shake256_inc_absorb(&ctx, message, length);
        shake256_inc_finalize(&ctx);
        shake256_inc_squeeze((void *)(*scalar), hash_bytes, &ctx);
        (*scalar)[limbs - 1] &= mask;
        for (int i = 2; i < HASH_ITERATIONS; i++) {
            shake256_inc_init(&ctx);
            shake256_inc_absorb(&ctx, (void *)(*scalar), hash_bytes);
            shake256_inc_finalize(&ctx);
            shake256_inc_squeeze((void *)(*scalar), hash_bytes, &ctx);
            (*scalar)[limbs - 1] &= mask;
        }
        shake256_inc_init(&ctx);
        shake256_inc_absorb(&ctx, (void *)(*scalar), hash_bytes);
        shake256_inc_finalize(&ctx);

        hash_bytes = ((TORSION_EVEN_POWER - SQIsign_response_length) + 7) / 8;
        limbs = (hash_bytes + sizeof(digit_t) - 1) / sizeof(digit_t);
        bits = (TORSION_EVEN_POWER - SQIsign_response_length) % RADIX;
        mask = ((digit_t)-1) >> ((RADIX - bits) % RADIX);
#ifdef TARGET_BIG_ENDIAN
        mask = BSWAP_DIGIT(mask);
#endif

        memset(*scalar, 0, NWORDS_ORDER * sizeof(digit_t));
        shake256_inc_squeeze((void *)(*scalar), hash_bytes, &ctx);
        (*scalar)[limbs - 1] &= mask;

#ifdef TARGET_BIG_ENDIAN
        for (int i = 0; i < NWORDS_ORDER; i++)
            (*scalar)[i] = BSWAP_DIGIT((*scalar)[i]);
#endif

        mp_mod_2exp(*scalar, SECURITY_BITS, NWORDS_ORDER);
    }
}