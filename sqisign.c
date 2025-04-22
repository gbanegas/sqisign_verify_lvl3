//
// Created by gustavo on 18/04/25.
//

#include "sqisign.h"
#include "ec_params.h"
#include "e0_basis.h"


int
sqisign_verify(const unsigned char *m,
               unsigned long long mlen,
               const unsigned char *sig,
               unsigned long long siglen,
               const unsigned char *pk) {
    int ret = 0;
    public_key_t pkt = {0};
    signature_t sigt;

    public_key_from_bytes(&pkt, pk);
    signature_from_bytes(&sigt, sig);

    ret = !protocols_verify(&sigt, &pkt, m, mlen);

    return ret;
}


void
xDBL_E0(ec_point_t *Q, const ec_point_t *P) {
    // Doubling of a Montgomery point in projective coordinates (X:Z) on the curve E0 with (A:C) = (0:1).
    // Input: projective Montgomery x-coordinates P = (XP:ZP), where xP=XP/ZP, and Montgomery curve constants (A:C) = (0:1).
    // Output: projective Montgomery x-coordinates Q <- 2*P = (XQ:ZQ) such that x(2P)=XQ/ZQ.
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t1, &t1, &t1);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_add(&Q->z, &t1, &t2);
    fp2_mul(&Q->z, &Q->z, &t2);
}

void
xDBL(ec_point_t *Q, const ec_point_t *P, const ec_point_t *AC) {
    // Doubling of a Montgomery point in projective coordinates (X:Z). Computation of coefficient values A+2C and 4C
    // on-the-fly.
    // Input: projective Montgomery x-coordinates P = (XP:ZP), where xP=XP/ZP, and Montgomery curve constants (A:C).
    // Output: projective Montgomery x-coordinates Q <- 2*P = (XQ:ZQ) such that x(2P)=XQ/ZQ.
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t3, &AC->z, &AC->z);
    fp2_mul(&t1, &t1, &t3);
    fp2_add(&t1, &t1, &t1);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_add(&t0, &t3, &AC->x);
    fp2_mul(&t0, &t0, &t2);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void
xDBL_A24(ec_point_t *Q, const ec_point_t *P, const ec_point_t *A24, const bool A24_normalized) {
    // Doubling of a Montgomery point in projective coordinates (X:Z).
    // Input: projective Montgomery x-coordinates P = (XP:ZP), where xP=XP/ZP, and
    //        the Montgomery curve constants A24 = (A+2C:4C) (or A24 = (A+2C/4C:1) if normalized).
    // Output: projective Montgomery x-coordinates Q <- 2*P = (XQ:ZQ) such that x(2P)=XQ/ZQ.
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    if (!A24_normalized)
        fp2_mul(&t1, &t1, &A24->z);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &A24->x);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void
xADD(ec_point_t *R, const ec_point_t *P, const ec_point_t *Q, const ec_point_t *PQ) {
    // Differential addition of Montgomery points in projective coordinates (X:Z).
    // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, and difference
    //        PQ=P-Q=(XPQ:ZPQ).
    // Output: projective Montgomery point R <- P+Q = (XR:ZR) such that x(P+Q)=XR/ZR.
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&t2, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &t2);
    fp2_add(&t2, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&t2, &t2);
    fp2_sqr(&t3, &t3);
    fp2_mul(&t2, &PQ->z, &t2);
    fp2_mul(&R->z, &PQ->x, &t3);
    fp2_copy(&R->x, &t2);
}

void
xDBLADD(ec_point_t *R,
        ec_point_t *S,
        const ec_point_t *P,
        const ec_point_t *Q,
        const ec_point_t *PQ,
        const ec_point_t *A24,
        const bool A24_normalized) {
    // Simultaneous doubling and differential addition.
    // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, the difference
    //         PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants A24 = (A+2C:4C) (or A24 = (A+2C/4C:1) if normalized).
    // Output: projective Montgomery points R <- 2*P = (XR:ZR) such that x(2P)=XR/ZR, and S <- P+Q = (XS:ZS) such that =
    //         x(Q+P)=XS/ZS.
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&R->x, &t0);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_add(&S->x, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t2);
    fp2_sqr(&R->z, &t1);
    fp2_mul(&t1, &t1, &S->x);
    fp2_sub(&t2, &R->x, &R->z);
    if (!A24_normalized)
        fp2_mul(&R->z, &R->z, &A24->z);
    fp2_mul(&R->x, &R->x, &R->z);
    fp2_mul(&S->x, &A24->x, &t2);
    fp2_sub(&S->z, &t0, &t1);
    fp2_add(&R->z, &R->z, &S->x);
    fp2_add(&S->x, &t0, &t1);
    fp2_mul(&R->z, &R->z, &t2);
    fp2_sqr(&S->z, &S->z);
    fp2_sqr(&S->x, &S->x);
    fp2_mul(&S->z, &S->z, &PQ->x);
    fp2_mul(&S->x, &S->x, &PQ->z);
}


void
ec_normalize_curve(ec_curve_t *E) {
    fp2_inv(&E->C);
    fp2_mul(&E->A, &E->A, &E->C);
    fp2_set_one(&E->C);
}

void
ec_normalize_curve_and_A24(ec_curve_t *E) {
    // Neither the curve or A24 are guaranteed to be normalized.
    // First we normalize (A/C : 1) and conditionally compute
    if (!fp2_is_one(&E->C)) {
        ec_normalize_curve(E);
    }

    if (!E->is_A24_computed_and_normalized) {
        // Now compute A24 = ((A + 2) / 4 : 1)
        fp2_add_one(&E->A24.x, &E->A); // re(A24.x) = re(A) + 1
        fp2_add_one(&E->A24.x, &E->A24.x); // re(A24.x) = re(A) + 2
        fp_copy(&E->A24.x.im, &E->A.im); // im(A24.x) = im(A)

        fp2_half(&E->A24.x, &E->A24.x); // (A + 2) / 2
        fp2_half(&E->A24.x, &E->A24.x); // (A + 2) / 4
        fp2_set_one(&E->A24.z);

        E->is_A24_computed_and_normalized = true;
    }
}

static void
difference_point(ec_point_t *PQ, const ec_point_t *P, const ec_point_t *Q, const ec_curve_t *curve) {
    // Given P,Q in projective x-only, computes a deterministic choice for (P-Q)
    // Based on Proposition 3 of https://eprint.iacr.org/2017/518.pdf

    fp2_t Bxx, Bxz, Bzz, t0, t1;

    fp2_mul(&t0, &P->x, &Q->x);
    fp2_mul(&t1, &P->z, &Q->z);
    fp2_sub(&Bxx, &t0, &t1);
    fp2_sqr(&Bxx, &Bxx);
    fp2_mul(&Bxx, &Bxx, &curve->C); // C*(P.x*Q.x-P.z*Q.z)^2
    fp2_add(&Bxz, &t0, &t1);
    fp2_mul(&t0, &P->x, &Q->z);
    fp2_mul(&t1, &P->z, &Q->x);
    fp2_add(&Bzz, &t0, &t1);
    fp2_mul(&Bxz, &Bxz, &Bzz); // (P.x*Q.x+P.z*Q.z)(P.x*Q.z+P.z*Q.x)
    fp2_sub(&Bzz, &t0, &t1);
    fp2_sqr(&Bzz, &Bzz);
    fp2_mul(&Bzz, &Bzz, &curve->C); // C*(P.x*Q.z-P.z*Q.x)^2
    fp2_mul(&Bxz, &Bxz, &curve->C); // C*(P.x*Q.x+P.z*Q.z)(P.x*Q.z+P.z*Q.x)
    fp2_mul(&t0, &t0, &t1);
    fp2_mul(&t0, &t0, &curve->A);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&Bxz, &Bxz, &t0); // C*(P.x*Q.x+P.z*Q.z)(P.x*Q.z+P.z*Q.x) + 2*A*P.x*Q.z*P.z*Q.x

    // To ensure that the denominator is a fourth power in Fp, we normalize by
    // C*C_bar^2*(P.z)_bar^2*(Q.z)_bar^2
    fp_copy(&t0.re, &curve->C.re);
    fp_neg(&t0.im, &curve->C.im);
    fp2_sqr(&t0, &t0);
    fp2_mul(&t0, &t0, &curve->C);
    fp_copy(&t1.re, &P->z.re);
    fp_neg(&t1.im, &P->z.im);
    fp2_sqr(&t1, &t1);
    fp2_mul(&t0, &t0, &t1);
    fp_copy(&t1.re, &Q->z.re);
    fp_neg(&t1.im, &Q->z.im);
    fp2_sqr(&t1, &t1);
    fp2_mul(&t0, &t0, &t1);
    fp2_mul(&Bxx, &Bxx, &t0);
    fp2_mul(&Bxz, &Bxz, &t0);
    fp2_mul(&Bzz, &Bzz, &t0);

    // Solving quadratic equation
    fp2_sqr(&t0, &Bxz);
    fp2_mul(&t1, &Bxx, &Bzz);
    fp2_sub(&t0, &t0, &t1);
    // No need to check if t0 is square, as per the entangled basis algorithm.
    fp2_sqrt(&t0);
    fp2_add(&PQ->x, &Bxz, &t0);
    fp2_copy(&PQ->z, &Bzz);
}


// Copying points, bases and curves
static inline void
copy_point(ec_point_t *P, const ec_point_t *Q) {
    fp2_copy(&P->x, &Q->x);
    fp2_copy(&P->z, &Q->z);
}

static inline void
copy_basis(ec_basis_t *B1, const ec_basis_t *B0) {
    copy_point(&B1->P, &B0->P);
    copy_point(&B1->Q, &B0->Q);
    copy_point(&B1->PmQ, &B0->PmQ);
}

static inline void
copy_curve(ec_curve_t *E1, const ec_curve_t *E2) {
    fp2_copy(&(E1->A), &(E2->A));
    fp2_copy(&(E1->C), &(E2->C));
    E1->is_A24_computed_and_normalized = E2->is_A24_computed_and_normalized;
    copy_point(&E1->A24, &E2->A24);
}


// The entangled basis generation does not allow A = 0
// so we simply return the one we have already precomputed
static void
ec_basis_E0_2f(ec_basis_t *PQ2, ec_curve_t *curve, int f) {
    ec_point_t P, Q;

    // Set P, Q to precomputed (X : 1) values
    fp2_copy(&P.x, &BASIS_E0_PX);
    fp2_copy(&Q.x, &BASIS_E0_QX);
    fp2_set_one(&P.z);
    fp2_set_one(&Q.z);

    // clear the power of two to get a point of order 2^f
    for (int i = 0; i < TORSION_EVEN_POWER - f; i++) {
        xDBL_E0(&P, &P);
        xDBL_E0(&Q, &Q);
    }

    // Set P, Q in the basis and compute x(P - Q)
    copy_point(&PQ2->P, &P);
    copy_point(&PQ2->Q, &Q);
    difference_point(&PQ2->PmQ, &P, &Q, curve);
}

// Given an x-coordinate, determines if this is a valid
// point on the curve. Assumes C=1.
static uint32_t
is_on_curve(const fp2_t *x, const ec_curve_t *curve) {
    fp2_t t0;

    fp2_add(&t0, x, &curve->A); // x + (A/C)
    fp2_mul(&t0, &t0, x); // x^2 + (A/C)*x
    fp2_add_one(&t0, &t0); // x^2 + (A/C)*x + 1
    fp2_mul(&t0, &t0, x); // x^3 + (A/C)*x^2 + x

    return fp2_is_square(&t0);
}

// Helper function which finds a point x(P) = n * A
static uint8_t
find_nA_x_coord(fp2_t *x, ec_curve_t *curve, const uint8_t start) {
    // when A is NQR we allow x(P) to be a multiple n*A of A
    uint8_t n = start;
    if (n == 1) {
        fp2_copy(x, &curve->A);
    } else {
        fp2_mul_small(x, &curve->A, n);
    }

    while (!is_on_curve(x, curve)) {
        fp2_add(x, x, &curve->A);
        n++;
    }

    /*
     * With very low probability (1/2^128), n will not fit in 7 bits.
     * In this case, we set hint = 0 which signals failure and the need
     * to generate a value on the fly during verification
     */
    uint8_t hint = n < 128 ? n : 0;
    return hint;
}

// Helper function which finds an NQR -1 / (1 + i*b) for entangled basis generation
static uint8_t
find_nqr_factor(fp2_t *x, ec_curve_t *curve, const uint8_t start) {
    // factor = -1/(1 + i*b) for b in Fp will be NQR whenever 1 + b^2 is NQR
    // in Fp, so we find one of these and then invert (1 + i*b). We store b
    // as a u8 hint to save time in verification.

    // We return the hint as a u8, but use (uint16_t)n to give 2^16 - 1
    // to make failure cryptographically negligible, with a fallback when
    // n > 128 is required.
    uint8_t hint;
    uint32_t found = 0;
    uint16_t n = start;

    bool qr_b = 1;
    fp_t b, tmp;
    fp2_t z, t0, t1;

    do {
        while (qr_b) {
            // find b with 1 + b^2 a non-quadratic residue
            fp_set_small(&tmp, (uint32_t) n * n + 1);
            qr_b = fp_is_square(&tmp);
            n++; // keeps track of b = n - 1
        }

        // for Px := -A/(1 + i*b) to be on the curve
        // is equivalent to A^2*(z-1) - z^2 NQR for z = 1 + i*b
        // thus prevents unnecessary inversion pre-check

        // t0 = z - 1 = i*b
        // t1 = z = 1 + i*b
        fp_set_small(&b, (uint32_t) n - 1);
        fp2_set_zero(&t0);
        fp2_set_one(&z);
        fp_copy(&z.im, &b);
        fp_copy(&t0.im, &b);

        // A^2*(z-1) - z^2
        fp2_sqr(&t1, &curve->A);
        fp2_mul(&t0, &t0, &t1); // A^2 * (z - 1)
        fp2_sqr(&t1, &z);
        fp2_sub(&t0, &t0, &t1); // A^2 * (z - 1) - z^2
        found = !fp2_is_square(&t0);

        qr_b = 1;
    } while (!found);

    // set Px to -A/(1 + i*b)
    fp2_copy(x, &z);
    fp2_inv(x);
    fp2_mul(x, x, &curve->A);
    fp2_neg(x, x);

    /*
     * With very low probability n will not fit in 7 bits.
     * We set hint = 0 which signals failure and the need
     * to generate a value on the fly during verification
     */
    hint = n <= 128 ? n - 1 : 0;

    return hint;
}

static inline void
AC_to_A24(ec_point_t *A24, const ec_curve_t *E) {
    // Maybe we already have this computed
    if (E->is_A24_computed_and_normalized) {
        copy_point(A24, &E->A24);
        return;
    }

    // A24 = (A+2C : 4C)
    fp2_add(&A24->z, &E->C, &E->C);
    fp2_add(&A24->x, &E->A, &A24->z);
    fp2_add(&A24->z, &A24->z, &A24->z);
}

void
ec_normalize_point(ec_point_t *P) {
    fp2_inv(&P->z);
    fp2_mul(&P->x, &P->x, &P->z);
    fp2_set_one(&(P->z));
}


void
ec_curve_normalize_A24(ec_curve_t *E) {
    if (!E->is_A24_computed_and_normalized) {
        AC_to_A24(&E->A24, E);
        ec_normalize_point(&E->A24);
        E->is_A24_computed_and_normalized = true;
    }
}

void
cswap_points(ec_point_t *P, ec_point_t *Q, const digit_t option) {
    // Swap points in constant time
    // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    fp2_cswap(&(P->x), &(Q->x), option);
    fp2_cswap(&(P->z), &(Q->z), option);
}


void
ec_point_init(ec_point_t *P) {
    // Initialize point as identity element (1:0)
    fp2_set_one(&(P->x));
    fp2_set_zero(&(P->z));
}

void
xMUL(ec_point_t *Q, const ec_point_t *P, const digit_t *k, const int kbits, const ec_curve_t *curve) {
    // The Montgomery ladder
    // Input: projective Montgomery point P=(XP:ZP) such that xP=XP/ZP, a scalar k of bitlength kbits, and
    //        the Montgomery curve constants (A:C) (or A24 = (A+2C/4C:1) if normalized).
    // Output: projective Montgomery points Q <- k*P = (XQ:ZQ) such that x(k*P)=XQ/ZQ.
    ec_point_t R0, R1, A24;
    digit_t mask;
    unsigned int bit, prevbit = 0, swap;

    if (!curve->is_A24_computed_and_normalized) {
        // Computation of A24=(A+2C:4C)
        fp2_add(&A24.x, &curve->C, &curve->C);
        fp2_add(&A24.z, &A24.x, &A24.x);
        fp2_add(&A24.x, &A24.x, &curve->A);
    } else {
        fp2_copy(&A24.x, &curve->A24.x);
        fp2_copy(&A24.z, &curve->A24.z);
        // Assert A24 has been normalised
    }

    // R0 <- (1:0), R1 <- P
    ec_point_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    for (int i = kbits - 1; i >= 0; i--) {
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX - 1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t) swap;

        cswap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, &A24, true);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t) swap;
    cswap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

void
ec_mul(ec_point_t *res, const digit_t *scalar, const int kbits, const ec_point_t *P, ec_curve_t *curve) {
    // For large scalars it's worth normalising anyway
    if (kbits > 50) {
        ec_curve_normalize_A24(curve);
    }

    // When A24 is computed and normalized we save some Fp2 multiplications
    xMUL(res, P, scalar, kbits, curve);
}


// Helper function which given a point of order k*2^n with n maximal
// and k odd, computes a point of order 2^f
static inline void
clear_cofactor_for_maximal_even_order(ec_point_t *P, ec_curve_t *curve, int f) {
    // clear out the odd cofactor to get a point of order 2^n
    ec_mul(P, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, P, curve);

    // clear the power of two to get a point of order 2^f
    for (int i = 0; i < TORSION_EVEN_POWER - f; i++) {
        xDBL_A24(P, P, &curve->A24, curve->is_A24_computed_and_normalized);
    }
}

uint32_t
ec_is_zero(const ec_point_t *P) {
    return fp2_is_zero(&P->z);
}


void
ec_dbl_iter(ec_point_t *res, int n, const ec_point_t *P, ec_curve_t *curve) {
    if (n == 0) {
        copy_point(res, P);
        return;
    }

    // When the chain is long enough, we should normalise A24
    if (n > 50) {
        ec_curve_normalize_A24(curve);
    }

    // When A24 is normalized we can save some multiplications
    if (curve->is_A24_computed_and_normalized) {
        xDBL_A24(res, P, &curve->A24, true);
        for (int i = 0; i < n - 1; i++) {
            xDBL_A24(res, res, &curve->A24, true);
        }
    } else {
        // Otherwise we do normal doubling
        xDBL(res, P, (const ec_point_t *) curve);
        for (int i = 0; i < n - 1; i++) {
            xDBL(res, res, (const ec_point_t *) curve);
        }
    }
}


void
ec_dbl(ec_point_t *res, const ec_point_t *P, const ec_curve_t *curve) {
    // If A24 = ((A+2)/4 : 1) we save multiplications
    if (curve->is_A24_computed_and_normalized) {
        xDBL_A24(res, P, &curve->A24, true);
    } else {
        // Otherwise we compute A24 on the fly for doubling
        xDBL(res, P, (const ec_point_t *) curve);
    }
}

/**
 * @brief Check if a point (X : Z) has order exactly 2^t
 *
 * @param P: a point
 * @param E: an elliptic curve
 * @param t: an integer
 *
 * @return 0xFFFFFFFF if the order is correct, 0 otherwise
 */
static int
test_point_order_twof(const ec_point_t *P, const ec_curve_t *E, int t) {
    ec_point_t test;
    ec_curve_t curve;
    test = *P;
    copy_curve(&curve, E);

    if (ec_is_zero(&test))
        return 0;
    // Scale point by 2^(t-1)
    ec_dbl_iter(&test, t - 1, &test, &curve);
    // If it's zero now, it doesnt have order 2^t
    if (ec_is_zero(&test))
        return 0;
    // Ensure [2^t] P = 0
    ec_dbl(&test, &test, &curve);
    return ec_is_zero(&test);
}


/**
 * @brief Check if basis points (P, Q, PmQ) all have order exactly 2^t
 *
 * @param B: a basis
 * @param E: an elliptic curve
 * @param t: an integer
 *
 * @return 0xFFFFFFFF if the order is correct, 0 otherwise
 */
static int
test_basis_order_twof(const ec_basis_t *B, const ec_curve_t *E, int t) {
    int check_P = test_point_order_twof(&B->P, E, t);
    int check_Q = test_point_order_twof(&B->Q, E, t);
    int check_PmQ = test_point_order_twof(&B->PmQ, E, t);

    return check_P & check_Q & check_PmQ;
}


// Computes a basis E[2^f] = <P, Q> where the point Q is above (0 : 0)
// given the hints as an array for faster basis computation
int
ec_curve_to_basis_2f_from_hint(ec_basis_t *PQ2, ec_curve_t *curve, int f, const uint8_t hint) {
    // Normalise (A/C : 1) and ((A + 2)/4 : 1)
    ec_normalize_curve_and_A24(curve);

    if (fp2_is_zero(&curve->A)) {
        ec_basis_E0_2f(PQ2, curve, f);
        return 1;
    }

    // The LSB of hint encodes whether A is a QR
    // The remaining 7-bits are used to find a valid x(P)
    bool hint_A = hint & 1;
    uint8_t hint_P = hint >> 1;

    // Compute the points P, Q
    ec_point_t P, Q;

    if (!hint_P) {
        // When hint_P = 0 it means we did not find a point in 128 attempts
        // this is very rare and we almost never expect to need this fallback
        // In either case, we can start with b = 128 to skip testing the known
        // values which will not work
        if (!hint_A) {
            find_nA_x_coord(&P.x, curve, 128);
        } else {
            find_nqr_factor(&P.x, curve, 128);
        }
    } else {
        // Otherwise we use the hint to directly find x(P) based on hint_A
        if (!hint_A) {
            // when A is NQR, we have found n such that x(P) = n*A
            fp2_mul_small(&P.x, &curve->A, hint_P);
        } else {
            // when A is QR we have found b such that (1 + b^2) is a NQR in
            // Fp, so we must compute x(P) = -A / (1 + i*b)
            fp_set_one(&P.x.re);
            fp_set_small(&P.x.im, hint_P);
            fp2_inv(&P.x);
            fp2_mul(&P.x, &P.x, &curve->A);
            fp2_neg(&P.x, &P.x);
        }
    }
    fp2_set_one(&P.z);

#ifndef NDEBUG
    int passed = 1;
    passed = is_on_curve(&P.x, curve);
    passed &= !fp2_is_square(&P.x);

    if (!passed)
        return 0;
#endif

    // set xQ to -xP - A
    fp2_add(&Q.x, &curve->A, &P.x);
    fp2_neg(&Q.x, &Q.x);
    fp2_set_one(&Q.z);

    // clear out the odd cofactor to get a point of order 2^f
    clear_cofactor_for_maximal_even_order(&P, curve, f);
    clear_cofactor_for_maximal_even_order(&Q, curve, f);

    // compute PmQ, set PmQ to Q to ensure Q above (0,0)
    difference_point(&PQ2->Q, &P, &Q, curve);
    copy_point(&PQ2->P, &P);
    copy_point(&PQ2->PmQ, &Q);

#ifndef NDEBUG
    passed &= test_basis_order_twof(PQ2, curve, f);

    if (!passed)
        return 0;
#endif

    return 1;
}


// Check that the basis change matrix elements are canonical
// representatives modulo 2^(SQIsign_response_length + 2).
static int
check_canonical_basis_change_matrix(const signature_t *sig) {
    // This works as long as all values in sig->mat_Bchall_can_to_B_chall are
    // positive integers.
    int ret = 1;
    scalar_t aux;

    memset(aux, 0, NWORDS_ORDER * sizeof(digit_t));
    aux[0] = 0x1;
    multiple_mp_shiftl(aux, SQIsign_response_length + HD_extra_torsion - (int) sig->backtracking, NWORDS_ORDER);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if (mp_compare(aux, sig->mat_Bchall_can_to_B_chall[i][j], NWORDS_ORDER) <= 0) {
                ret = 0;
            }
        }
    }

    return ret;
}

/*
// Copying points, bases and curves
static inline void
copy_point(ec_point_t *P, const ec_point_t *Q)
{
    fp2_copy(&P->x, &Q->x);
    fp2_copy(&P->z, &Q->z);
}

static inline void
copy_curve(ec_curve_t *E1, const ec_curve_t *E2)
{
    fp2_copy(&(E1->A), &(E2->A));
    fp2_copy(&(E1->C), &(E2->C));
    E1->is_A24_computed_and_normalized = E2->is_A24_computed_and_normalized;
    copy_point(&E1->A24, &E2->A24);
}
*/
uint32_t
ec_has_zero_coordinate(const ec_point_t *P) {
    return fp2_is_zero(&P->x) | fp2_is_zero(&P->z);
}


int
ec_ladder3pt(ec_point_t *R,
             const digit_t *m,
             const ec_point_t *P,
             const ec_point_t *Q,
             const ec_point_t *PQ,
             const ec_curve_t *E) {
    // The 3-point Montgomery ladder
    // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, a scalar k of
    //         bitlength kbits, the difference PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants A24 = (A+2C/4C:1).
    // Output: projective Montgomery point R <- P + m*Q = (XR:ZR) such that x(P + m*Q)=XR/ZR.
    //assert(E->is_A24_computed_and_normalized);
    if (!fp2_is_one(&E->A24.z)) {
        return 0;
    }
    // Formulas are not valid in that case
    if (ec_has_zero_coordinate(PQ)) {
        return 0;
    }

    ec_point_t X0, X1, X2;
    copy_point(&X0, Q);
    copy_point(&X1, P);
    copy_point(&X2, PQ);

    int i, j;
    digit_t t;
    for (i = 0; i < NWORDS_ORDER; i++) {
        t = 1;
        for (j = 0; j < RADIX; j++) {
            cswap_points(&X1, &X2, -((t & m[i]) == 0));
            xDBLADD(&X0, &X1, &X0, &X1, &X2, &E->A24, true);
            cswap_points(&X1, &X2, -((t & m[i]) == 0));
            t <<= 1;
        };
    };
    copy_point(R, &X1);
    return 1;
}

uint32_t
ec_is_two_torsion(const ec_point_t *P, const ec_curve_t *E) {
    if (ec_is_zero(P))
        return 0;

    uint32_t x_is_zero, tmp_is_zero;
    fp2_t t0, t1, t2;
    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t1, &t0, &t1);
    fp2_mul(&t2, &t2, &E->A);
    fp2_mul(&t1, &t1, &E->C);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t0, &t1, &t2); // 4 (CX^2+CZ^2+AXZ)

    x_is_zero = fp2_is_zero(&P->x);
    tmp_is_zero = fp2_is_zero(&t0);

    // two torsion if x or x^2 + Ax + 1 is zero
    return x_is_zero | tmp_is_zero;
}


uint32_t
ec_is_four_torsion(const ec_point_t *P, const ec_curve_t *E) {
    ec_point_t test;
    xDBL_A24(&test, P, &E->A24, E->is_A24_computed_and_normalized);
    return ec_is_two_torsion(&test, E);
}

// Degree-4 isogeny with kernel generated by P such that [2]P != (0 ,0)
// Outputs the curve coefficient in the form A24=(A+2C:4C)
void
xisog_4(ec_kps4_t *kps, ec_point_t *B, const ec_point_t P) {
    ec_point_t *K = kps->K;

    fp2_sqr(&K[0].x, &P.x);
    fp2_sqr(&K[0].z, &P.z);
    fp2_add(&K[1].x, &K[0].z, &K[0].x);
    fp2_sub(&K[1].z, &K[0].z, &K[0].x);
    fp2_mul(&B->x, &K[1].x, &K[1].z);
    fp2_sqr(&B->z, &K[0].z);

    // Constants for xeval_4
    fp2_add(&K[2].x, &P.x, &P.z);
    fp2_sub(&K[1].x, &P.x, &P.z);
    fp2_add(&K[0].x, &K[0].z, &K[0].z);
    fp2_add(&K[0].x, &K[0].x, &K[0].x);
}


// Degree-4 isogeny evaluation with kenerl generated by P such that [2]P != (0, 0)
void
xeval_4(ec_point_t *R, const ec_point_t *Q, const int lenQ, const ec_kps4_t *kps) {
    const ec_point_t *K = kps->K;

    fp2_t t0, t1;

    for (int i = 0; i < lenQ; i++) {
        fp2_add(&t0, &Q[i].x, &Q[i].z);
        fp2_sub(&t1, &Q[i].x, &Q[i].z);
        fp2_mul(&(R[i].x), &t0, &K[1].x);
        fp2_mul(&(R[i].z), &t1, &K[2].x);
        fp2_mul(&t0, &t0, &t1);
        fp2_mul(&t0, &t0, &K[0].x);
        fp2_add(&t1, &(R[i].x), &(R[i].z));
        fp2_sub(&(R[i].z), &(R[i].x), &(R[i].z));
        fp2_sqr(&t1, &t1);
        fp2_sqr(&(R[i].z), &(R[i].z));
        fp2_add(&(R[i].x), &t0, &t1);
        fp2_sub(&t0, &t0, &(R[i].z));
        fp2_mul(&(R[i].x), &(R[i].x), &t1);
        fp2_mul(&(R[i].z), &(R[i].z), &t0);
    }
}

// Degree-2 isogeny with kernel generated by P != (0 ,0)
// Outputs the curve coefficient in the form A24=(A+2C:4C)
void
xisog_2(ec_kps2_t *kps, ec_point_t *B, const ec_point_t P) {
    fp2_sqr(&B->x, &P.x);
    fp2_sqr(&B->z, &P.z);
    fp2_sub(&B->x, &B->z, &B->x);
    fp2_add(&kps->K.x, &P.x, &P.z);
    fp2_sub(&kps->K.z, &P.x, &P.z);
}

void
xeval_2(ec_point_t *R, ec_point_t *const Q, const int lenQ, const ec_kps2_t *kps) {
    fp2_t t0, t1, t2;
    for (int j = 0; j < lenQ; j++) {
        fp2_add(&t0, &Q[j].x, &Q[j].z);
        fp2_sub(&t1, &Q[j].x, &Q[j].z);
        fp2_mul(&t2, &kps->K.x, &t1);
        fp2_mul(&t1, &kps->K.z, &t0);
        fp2_add(&t0, &t2, &t1);
        fp2_sub(&t1, &t2, &t1);
        fp2_mul(&R[j].x, &Q[j].x, &t0);
        fp2_mul(&R[j].z, &Q[j].z, &t1);
    }
}

/**
 * @brief Given a curve the point (A+2 : 4C) compute the curve coefficients (A : C)
 *
 * @param E a curve to compute
 * @param A24 the value (A+2 : 4C)
 */
static inline void
A24_to_AC(ec_curve_t *E, const ec_point_t *A24) {
    // (A:C) = ((A+2C)*2-4C : 4C)
    fp2_add(&E->A, &A24->x, &A24->x);
    fp2_sub(&E->A, &E->A, &A24->z);
    fp2_add(&E->A, &E->A, &E->A);
    fp2_copy(&E->C, &A24->z);
}


// since we use degree 4 isogeny steps, we need to handle the odd case with care
static uint32_t
ec_eval_even_strategy(ec_curve_t *curve,
                      ec_point_t *points,
                      unsigned len_points,
                      const ec_point_t *kernel,
                      const int isog_len) {
    ec_curve_normalize_A24(curve);
    ec_point_t A24;
    copy_point(&A24, &curve->A24);

    int space = 1;
    for (int i = 1; i < isog_len; i *= 2)
        ++space;

    // Stack of remaining kernel points and their associated orders
    ec_point_t splits[space];
    uint16_t todo[space];
    splits[0] = *kernel;
    todo[0] = isog_len;

    int current = 0; // Pointer to current top of stack

    // Chain of 4-isogenies
    for (int j = 0; j < isog_len / 2; ++j) {
        // assert(current >= 0);
        // assert(todo[current] >= 1);
        // Get the next point of order 4
        while (todo[current] != 2) {
            //assert(todo[current] >= 3);
            // A new split will be added
            ++current;
            //assert(current < space);
            // We set the seed of the new split to be computed and saved
            copy_point(&splits[current], &splits[current - 1]);
            // if we copied from the very first element, then we perform one additional doubling
            unsigned num_dbls = todo[current - 1] / 4 * 2 + todo[current - 1] % 2;
            todo[current] = todo[current - 1] - num_dbls;
            while (num_dbls--)
                xDBL_A24(&splits[current], &splits[current], &A24, false);
        }

        if (j == 0) {
            // assert(fp2_is_one(&A24.z));
            if (!ec_is_four_torsion(&splits[current], curve))
                return -1;

            ec_point_t T;
            xDBL_A24(&T, &splits[current], &A24, false);
            if (fp2_is_zero(&T.x))
                return -1; // special isogenies not allowed
        } else {
            // assert(todo[current] == 2);
#ifndef NDEBUG
            if (fp2_is_zero(&splits[current].z))
                debug_print("splitting point z coordinate is unexpectedly zero");

            ec_point_t test;
            xDBL_A24(&test, &splits[current], &A24, false);
            if (fp2_is_zero(&test.z))
                debug_print("z coordinate is unexpectedly zero before doubling");
            xDBL_A24(&test, &test, &A24, false);
            if (!fp2_is_zero(&test.z))
                debug_print("z coordinate is unexpectedly not zero after doubling");
#endif
        }

        // Evaluate 4-isogeny
        ec_kps4_t kps4;
        xisog_4(&kps4, &A24, splits[current]);
        xeval_4(splits, splits, current, &kps4);
        for (int i = 0; i < current; ++i)
            todo[i] -= 2;
        xeval_4(points, points, len_points, &kps4);

        --current;
    }
    //  assert(isog_len % 2 ? !current : current == -1);

    // Final 2-isogeny
    if (isog_len % 2) {
#ifndef NDEBUG
        if (fp2_is_zero(&splits[0].z))
            debug_print("splitting point z coordinate is unexpectedly zero");
        ec_point_t test;
        copy_point(&test, &splits[0]);
        xDBL_A24(&test, &test, &A24, false);
        if (!fp2_is_zero(&test.z))
            debug_print("z coordinate is unexpectedly not zero after doubling");
#endif

        // We need to check the order of this point in case there were no 4-isogenies
        if (isog_len == 1 && !ec_is_two_torsion(&splits[0], curve))
            return -1;
        if (fp2_is_zero(&splits[0].x)) {
            // special isogenies not allowed
            // this case can only happen if isog_len == 1; otherwise the
            // previous 4-isogenies we computed ensure that $T=(0:1)$ is put
            // as the kernel of the dual isogeny
            return -1;
        }

        ec_kps2_t kps2;
        xisog_2(&kps2, &A24, splits[0]);
        xeval_2(points, points, len_points, &kps2);
    }

    // Output curve in the form (A:C)
    A24_to_AC(curve, &A24);

    curve->is_A24_computed_and_normalized = false;

    return 0;
}

uint32_t
ec_eval_even(ec_curve_t *image, ec_isog_even_t *phi, ec_point_t *points, unsigned len_points) {
    copy_curve(image, &phi->curve);
    return ec_eval_even_strategy(image, points, len_points, &phi->kernel, phi->length);
}

// Compute the 2^n isogeny from the signature with kernel
// P + [chall_coeff]Q and store the codomain in E_chall
static int
compute_challenge_verify(ec_curve_t *E_chall, const signature_t *sig, const ec_curve_t *Epk, const uint8_t hint_pk) {
    ec_basis_t bas_EA;
    ec_isog_even_t phi_chall;

    // Set domain and length of 2^n isogeny
    copy_curve(&phi_chall.curve, Epk);
    phi_chall.length = TORSION_EVEN_POWER - sig->backtracking;

    // Compute the basis from the supplied hint
    if (!ec_curve_to_basis_2f_from_hint(&bas_EA, &phi_chall.curve, TORSION_EVEN_POWER, hint_pk)) // canonical
        return 0;

    // recovering the exact challenge
    {
        if (!ec_ladder3pt(&phi_chall.kernel, sig->chall_coeff, &bas_EA.P, &bas_EA.Q, &bas_EA.PmQ, &phi_chall.curve)) {
            return 0;
        };
    }

    // Double the kernel until is has the correct order
    ec_dbl_iter(&phi_chall.kernel, sig->backtracking, &phi_chall.kernel, &phi_chall.curve);

    // Compute the codomain
    copy_curve(E_chall, &phi_chall.curve);
    if (ec_eval_even(E_chall, &phi_chall, NULL, 0))
        return 0;
    return 1;
}

void
select_point(ec_point_t *Q, const ec_point_t *P1, const ec_point_t *P2, const digit_t option)
{ // Select points in constant time
    // If option = 0 then Q <- P1, else if option = 0xFF...FF then Q <- P2
    fp2_select(&(Q->x), &(P1->x), &(P2->x), option);
    fp2_select(&(Q->z), &(P1->z), &(P2->z), option);
}




int
xDBLMUL(ec_point_t *S,
        const ec_point_t *P,
        const digit_t *k,
        const ec_point_t *Q,
        const digit_t *l,
        const ec_point_t *PQ,
        const int kbits,
        const ec_curve_t *curve)
{ // The Montgomery biladder
  // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, scalars k and l of
  //         bitlength kbits, the difference PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants (A:C).
  // Output: projective Montgomery point S <- k*P + l*Q = (XS:ZS) such that x(k*P + l*Q)=XS/ZS.

    int i, A_is_zero;
    digit_t evens, mevens, bitk0, bitl0, maskk, maskl, temp, bs1_ip1, bs2_ip1, bs1_i, bs2_i, h;
    digit_t sigma[2] = { 0 }, pre_sigma = 0;
    digit_t k_t[NWORDS_ORDER], l_t[NWORDS_ORDER], one[NWORDS_ORDER] = { 0 }, r[2 * BITS] = { 0 };
    ec_point_t DIFF1a, DIFF1b, DIFF2a, DIFF2b, R[3] = { 0 }, T[3];

    // differential additions formulas are invalid in this case
    if (ec_has_zero_coordinate(P) | ec_has_zero_coordinate(Q) | ec_has_zero_coordinate(PQ))
        return 0;

    // Derive sigma according to parity
    bitk0 = (k[0] & 1);
    bitl0 = (l[0] & 1);
    maskk = 0 - bitk0; // Parity masks: 0 if even, otherwise 1...1
    maskl = 0 - bitl0;
    sigma[0] = (bitk0 ^ 1);
    sigma[1] = (bitl0 ^ 1);
    evens = sigma[0] + sigma[1]; // Count number of even scalars
    mevens = 0 - (evens & 1);    // Mask mevens <- 0 if # even of scalars = 0 or 2, otherwise mevens = 1...1

    // If k and l are both even or both odd, pick sigma = (0,1)
    sigma[0] = (sigma[0] & mevens);
    sigma[1] = (sigma[1] & mevens) | (1 & ~mevens);

    // Convert even scalars to odd
    one[0] = 1;
    mp_sub(k_t, k, one, NWORDS_ORDER);
    mp_sub(l_t, l, one, NWORDS_ORDER);
    select_ct(k_t, k_t, k, maskk, NWORDS_ORDER);
    select_ct(l_t, l_t, l, maskl, NWORDS_ORDER);

    // Scalar recoding
    for (i = 0; i < kbits; i++) {
        // If sigma[0] = 1 swap k_t and l_t
        maskk = 0 - (sigma[0] ^ pre_sigma);
        swap_ct(k_t, l_t, maskk, NWORDS_ORDER);

        if (i == kbits - 1) {
            bs1_ip1 = 0;
            bs2_ip1 = 0;
        } else {
            bs1_ip1 = mp_shiftr(k_t, 1, NWORDS_ORDER);
            bs2_ip1 = mp_shiftr(l_t, 1, NWORDS_ORDER);
        }
        bs1_i = k_t[0] & 1;
        bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        // Revert sigma if second bit, r_(2i+1), is 1
        pre_sigma = sigma[0];
        maskk = 0 - r[2 * i + 1];
        select_ct(&temp, &sigma[0], &sigma[1], maskk, 1);
        select_ct(&sigma[1], &sigma[1], &sigma[0], maskk, 1);
        sigma[0] = temp;
    }

    // Point initialization
    ec_point_init(&R[0]);
    maskk = 0 - sigma[0];
    select_point(&R[1], P, Q, maskk);
    select_point(&R[2], Q, P, maskk);

    fp2_copy(&DIFF1a.x, &R[1].x);
    fp2_copy(&DIFF1a.z, &R[1].z);
    fp2_copy(&DIFF1b.x, &R[2].x);
    fp2_copy(&DIFF1b.z, &R[2].z);

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    xADD(&R[2], &R[1], &R[2], PQ);
    if (ec_has_zero_coordinate(&R[2]))
        return 0; // non valid formulas

    fp2_copy(&DIFF2a.x, &R[2].x);
    fp2_copy(&DIFF2a.z, &R[2].z);
    fp2_copy(&DIFF2b.x, &PQ->x);
    fp2_copy(&DIFF2b.z, &PQ->z);

    A_is_zero = fp2_is_zero(&curve->A);

    // Main loop
    for (i = kbits - 1; i >= 0; i--) {
        h = r[2 * i] + r[2 * i + 1]; // in {0, 1, 2}
        maskk = 0 - (h & 1);
        select_point(&T[0], &R[0], &R[1], maskk);
        maskk = 0 - (h >> 1);
        select_point(&T[0], &T[0], &R[2], maskk);
        if (A_is_zero) {
            xDBL_E0(&T[0], &T[0]);
        } else {
            //assert(fp2_is_one(&curve->A24.z));
            xDBL_A24(&T[0], &T[0], &curve->A24, true);
        }

        maskk = 0 - r[2 * i + 1]; // in {0, 1}
        select_point(&T[1], &R[0], &R[1], maskk);
        select_point(&T[2], &R[1], &R[2], maskk);

        cswap_points(&DIFF1a, &DIFF1b, maskk);
        xADD(&T[1], &T[1], &T[2], &DIFF1a);
        xADD(&T[2], &R[0], &R[2], &DIFF2a);

        // If hw (mod 2) = 1 then swap DIFF2a and DIFF2b
        maskk = 0 - (h & 1);
        cswap_points(&DIFF2a, &DIFF2b, maskk);

        // R <- T
        copy_point(&R[0], &T[0]);
        copy_point(&R[1], &T[1]);
        copy_point(&R[2], &T[2]);
    }

    // Output R[evens]
    select_point(S, &R[0], &R[1], mevens);

    maskk = 0 - (bitk0 & bitl0);
    select_point(S, S, &R[2], maskk);
    return 1;
}

int
ec_biscalar_mul(ec_point_t *res,
                const digit_t *scalarP,
                const digit_t *scalarQ,
                const int kbits,
                const ec_basis_t *PQ,
                const ec_curve_t *curve) {
    if (fp2_is_zero(&PQ->PmQ.z))
        return 0;

    /* Differential additions behave badly when PmQ = (0:1), so we need to
     * treat this case specifically. Since we assume P, Q are a basis, this
     * can happen only if kbits==1 */
    if (kbits == 1) {
        // Sanity check: our basis should be given by 2-torsion points
        if (!ec_is_two_torsion(&PQ->P, curve) || !ec_is_two_torsion(&PQ->Q, curve) ||
            !ec_is_two_torsion(&PQ->PmQ, curve))
            return 0;
        digit_t bP, bQ;
        bP = (scalarP[0] & 1);
        bQ = (scalarQ[0] & 1);
        if (bP == 0 && bQ == 0) {
            ec_point_init(res); //(1: 0)
        } else {
            if (bP == 1 && bQ == 0) {
                copy_point(res, &PQ->P);
            } else {
                if (bP == 0 && bQ == 1) {
                    copy_point(res, &PQ->Q);
                } else {
                    if (bP == 1 && bQ == 1) {
                        copy_point(res, &PQ->PmQ);
                    }
                }
            }
            //else // should never happen
            //assert(0);
        }
        return 1;
    } else {
        ec_curve_t E;
        copy_curve(&E, curve);

        if (!fp2_is_zero(&curve->A)) {
            // If A is not zero normalize
            ec_curve_normalize_A24(&E);
        }
        return xDBLMUL(res, &PQ->P, scalarP, &PQ->Q, scalarQ, &PQ->PmQ, kbits, (const ec_curve_t *) &E);
    }
}

// same as matrix_application_even_basis() in id2iso.c, with some modifications:
// - this version works with a matrix of scalars (not ibz_t).
// - reduction modulo 2^f of matrix elements is removed here, because it is
//   assumed that the elements are already cannonical representatives modulo
//   2^f; this is ensured by calling check_canonical_basis_change_matrix() at
//   the beginning of protocols_verify().
static int
matrix_scalar_application_even_basis(ec_basis_t *bas, const ec_curve_t *E, scalar_mtx_2x2_t *mat, int f) {
    scalar_t scalar0, scalar1;
    memset(scalar0, 0, NWORDS_ORDER * sizeof(digit_t));
    memset(scalar1, 0, NWORDS_ORDER * sizeof(digit_t));

    ec_basis_t tmp_bas;
    copy_basis(&tmp_bas, bas);

    // For a matrix [[a, c], [b, d]] we compute:
    //
    // first basis element R = [a]P + [b]Q
    if (!ec_biscalar_mul(&bas->P, (*mat)[0][0], (*mat)[1][0], f, &tmp_bas, E))
        return 0;
    // second basis element S = [c]P + [d]Q
    if (!ec_biscalar_mul(&bas->Q, (*mat)[0][1], (*mat)[1][1], f, &tmp_bas, E))
        return 0;
    // Their difference R - S = [a - c]P + [b - d]Q
    mp_sub(scalar0, (*mat)[0][0], (*mat)[0][1], NWORDS_ORDER);
    mp_mod_2exp(scalar0, f, NWORDS_ORDER);
    mp_sub(scalar1, (*mat)[1][0], (*mat)[1][1], NWORDS_ORDER);
    mp_mod_2exp(scalar1, f, NWORDS_ORDER);
    return ec_biscalar_mul(&bas->PmQ, scalar0, scalar1, f, &tmp_bas, E);
}

void
ec_dbl_iter_basis(ec_basis_t *res, int n, const ec_basis_t *B, ec_curve_t *curve)
{
    ec_dbl_iter(&res->P, n, &B->P, curve);
    ec_dbl_iter(&res->Q, n, &B->Q, curve);
    ec_dbl_iter(&res->PmQ, n, &B->PmQ, curve);
}
// Compute the bases for the challenge and auxillary curve from
// the canonical bases. Challenge basis is reconstructed from the
// compressed scalars within the challenge.
static int
challenge_and_aux_basis_verify(ec_basis_t *B_chall_can,
                               ec_basis_t *B_aux_can,
                               ec_curve_t *E_chall,
                               ec_curve_t *E_aux,
                               signature_t *sig,
                               const int pow_dim2_deg_resp) {
    // recovering the canonical basis as TORSION_EVEN_POWER for consistency with signing
    if (!ec_curve_to_basis_2f_from_hint(B_chall_can, E_chall, TORSION_EVEN_POWER, sig->hint_chall))
        return 0;

    // setting to the right order
    ec_dbl_iter_basis(B_chall_can,
                      TORSION_EVEN_POWER - pow_dim2_deg_resp - HD_extra_torsion - sig->two_resp_length,
                      B_chall_can,
                      E_chall);

    if (!ec_curve_to_basis_2f_from_hint(B_aux_can, E_aux, TORSION_EVEN_POWER, sig->hint_aux))
        return 0;

    // setting to the right order
    ec_dbl_iter_basis(B_aux_can, TORSION_EVEN_POWER - pow_dim2_deg_resp - HD_extra_torsion, B_aux_can, E_aux);

#ifndef NDEBUG
    if (!test_basis_order_twof(B_chall_can, E_chall, HD_extra_torsion + pow_dim2_deg_resp + sig->two_resp_length))
        debug_print("canonical basis has wrong order, expect something to fail");
#endif

    // applying the change matrix on the basis of E_chall
    return matrix_scalar_application_even_basis(B_chall_can,
                                                E_chall,
                                                &sig->mat_Bchall_can_to_B_chall,
                                                pow_dim2_deg_resp + HD_extra_torsion + sig->two_resp_length);
}



void
xisog_2_singular(ec_kps2_t *kps, ec_point_t *B24, ec_point_t A24)
{
    // No need to check the square root, only used for signing.
    fp2_t t0, four;
    fp2_set_small(&four, 4);
    fp2_add(&t0, &A24.x, &A24.x);
    fp2_sub(&t0, &t0, &A24.z);
    fp2_add(&t0, &t0, &t0);
    fp2_inv(&A24.z);
    fp2_mul(&t0, &t0, &A24.z);
    fp2_copy(&kps->K.x, &t0);
    fp2_add(&B24->x, &t0, &t0);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t0, &t0, &four);
    fp2_sqrt(&t0);
    fp2_neg(&kps->K.z, &t0);
    fp2_add(&B24->z, &t0, &t0);
    fp2_add(&B24->x, &B24->x, &B24->z);
    fp2_add(&B24->z, &B24->z, &B24->z);
}

void
xeval_2_singular(ec_point_t *R, const ec_point_t *Q, const int lenQ, const ec_kps2_t *kps)
{
    fp2_t t0, t1;
    for (int i = 0; i < lenQ; i++) {
        fp2_mul(&t0, &Q[i].x, &Q[i].z);
        fp2_mul(&t1, &kps->K.x, &Q[i].z);
        fp2_add(&t1, &t1, &Q[i].x);
        fp2_mul(&t1, &t1, &Q[i].x);
        fp2_sqr(&R[i].x, &Q[i].z);
        fp2_add(&R[i].x, &R[i].x, &t1);
        fp2_mul(&R[i].z, &t0, &kps->K.z);
    }
}



// naive implementation
uint32_t
ec_eval_small_chain(ec_curve_t *curve,
                    const ec_point_t *kernel,
                    int len,
                    ec_point_t *points,
                    unsigned len_points,
                    bool special) // do we allow special isogenies?
{

    ec_point_t A24;
    AC_to_A24(&A24, curve);

    ec_kps2_t kps;
    ec_point_t small_K, big_K;
    copy_point(&big_K, kernel);

    for (int i = 0; i < len; i++) {
        copy_point(&small_K, &big_K);
        // small_K = big_K;
        for (int j = 0; j < len - i - 1; j++) {
            xDBL_A24(&small_K, &small_K, &A24, false);
        }
        // Check the order of the point before the first isogeny step
        if (i == 0 && !ec_is_two_torsion(&small_K, curve))
            return (uint32_t)-1;
        // Perform isogeny step
        if (fp2_is_zero(&small_K.x)) {
            if (special) {
                ec_point_t B24;
                xisog_2_singular(&kps, &B24, A24);
                xeval_2_singular(&big_K, &big_K, 1, &kps);
                xeval_2_singular(points, points, len_points, &kps);
                copy_point(&A24, &B24);
            } else {
                return (uint32_t)-1;
            }
        } else {
            xisog_2(&kps, &A24, small_K);
            xeval_2(&big_K, &big_K, 1, &kps);
            xeval_2(points, points, len_points, &kps);
        }
    }
    A24_to_AC(curve, &A24);

    curve->is_A24_computed_and_normalized = false;
    return 0;
}

// When two_resp_length is non-zero, we must compute a small 2^n-isogeny
// updating E_chall as the codomain as well as push the basis on E_chall
// through this isogeny
static int
two_response_isogeny_verify(ec_curve_t *E_chall, ec_basis_t *B_chall_can, const signature_t *sig,
                            int pow_dim2_deg_resp) {
    ec_point_t ker, points[3];

    // choosing the right point for the small two_isogenies
    if (mp_is_even(sig->mat_Bchall_can_to_B_chall[0][0], NWORDS_ORDER) &&
        mp_is_even(sig->mat_Bchall_can_to_B_chall[1][0], NWORDS_ORDER)) {
        copy_point(&ker, &B_chall_can->Q);
    } else {
        copy_point(&ker, &B_chall_can->P);
    }

    copy_point(&points[0], &B_chall_can->P);
    copy_point(&points[1], &B_chall_can->Q);
    copy_point(&points[2], &B_chall_can->PmQ);

    ec_dbl_iter(&ker, pow_dim2_deg_resp + HD_extra_torsion, &ker, E_chall);

#ifndef NDEBUG
    if (!test_point_order_twof(&ker, E_chall, sig->two_resp_length))
        debug_print("kernel does not have order 2^(two_resp_length");
#endif

    if (ec_eval_small_chain(E_chall, &ker, sig->two_resp_length, points, 3, false)) {
        return 0;
    }

#ifndef NDEBUG
    if (!test_point_order_twof(&points[0], E_chall, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("points[0] does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
    if (!test_point_order_twof(&points[1], E_chall, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("points[1] does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
    if (!test_point_order_twof(&points[2], E_chall, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("points[2] does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
#endif

    copy_point(&B_chall_can->P, &points[0]);
    copy_point(&B_chall_can->Q, &points[1]);
    copy_point(&B_chall_can->PmQ, &points[2]);
    return 1;
}


void
copy_bases_to_kernel(theta_kernel_couple_points_t *ker, const ec_basis_t *B1, const ec_basis_t *B2)
{
    // Copy the basis on E1 to (P, _) on T1, T2 and T1 - T2
    copy_point(&ker->T1.P1, &B1->P);
    copy_point(&ker->T2.P1, &B1->Q);
    copy_point(&ker->T1m2.P1, &B1->PmQ);

    // Copy the basis on E2 to (_, P) on T1, T2 and T1 - T2
    copy_point(&ker->T1.P2, &B2->P);
    copy_point(&ker->T2.P2, &B2->Q);
    copy_point(&ker->T1m2.P2, &B2->PmQ);
}

void
ec_curve_init(ec_curve_t *E)
{ // Initialize the curve struct
    // Initialize the constants
    fp2_set_zero(&(E->A));
    fp2_set_one(&(E->C));

    // Initialize the point (A+2 : 4C)
    ec_point_init(&(E->A24));

    // Set the bool to be false by default
    E->is_A24_computed_and_normalized = false;
}

uint32_t
ec_is_equal(const ec_point_t *P, const ec_point_t *Q)
{ // Evaluate if two points in Montgomery coordinates (X:Z) are equal
    // Returns 0xFFFFFFFF (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1;

    // Check if P, Q are the points at infinity
    uint32_t l_zero = ec_is_zero(P);
    uint32_t r_zero = ec_is_zero(Q);

    // Check if PX * QZ = QX * PZ
    fp2_mul(&t0, &P->x, &Q->z);
    fp2_mul(&t1, &P->z, &Q->x);
    uint32_t lr_equal = fp2_is_equal(&t0, &t1);

    // Points are equal if
    // - Both are zero, or
    // - neither are zero AND PX * QZ = QX * PZ
    return (l_zero & r_zero) | (~l_zero & ~r_zero * lr_equal);
}

uint32_t
ec_is_basis_four_torsion(const ec_basis_t *B, const ec_curve_t *E)
{ // Check if basis points (P, Q) form a full 2^t-basis
    ec_point_t P2, Q2;
    xDBL_A24(&P2, &B->P, &E->A24, E->is_A24_computed_and_normalized);
    xDBL_A24(&Q2, &B->Q, &E->A24, E->is_A24_computed_and_normalized);
    return (ec_is_two_torsion(&P2, E) & ec_is_two_torsion(&Q2, E) & ~ec_is_equal(&P2, &Q2));
}



uint32_t
ec_recover_y(fp2_t *y, const fp2_t *Px, const ec_curve_t *curve)
{ // Recover y-coordinate of a point on the Montgomery curve y^2 = x^3 + Ax^2 + x
    fp2_t t0;

    fp2_sqr(&t0, Px);
    fp2_mul(y, &t0, &curve->A); // Ax^2
    fp2_add(y, y, Px);          // Ax^2 + x
    fp2_mul(&t0, &t0, Px);
    fp2_add(y, y, &t0); // x^3 + Ax^2 + x
    // This is required, because we do not yet know that our curves are
    // supersingular so our points live on the twist with B = 1.
    return fp2_sqrt_verify(y);
}


// Lifts a basis x(P), x(Q), x(P-Q) assuming the curve has (A/C : 1) and the point
// P = (X/Z : 1). For generic implementation see lift_basis()
uint32_t
lift_basis_normalized(jac_point_t *P, jac_point_t *Q, ec_basis_t *B, ec_curve_t *E)
{
    /*assert(fp2_is_one(&B->P.z));
    assert(fp2_is_one(&E->C));*/

    fp2_copy(&P->x, &B->P.x);
    fp2_copy(&Q->x, &B->Q.x);
    fp2_copy(&Q->z, &B->Q.z);
    fp2_set_one(&P->z);
    uint32_t ret = ec_recover_y(&P->y, &P->x, E);

    // Algorithm of Okeya-Sakurai to recover y.Q in the montgomery model
    fp2_t v1, v2, v3, v4;
    fp2_mul(&v1, &P->x, &Q->z);
    fp2_add(&v2, &Q->x, &v1);
    fp2_sub(&v3, &Q->x, &v1);
    fp2_sqr(&v3, &v3);
    fp2_mul(&v3, &v3, &B->PmQ.x);
    fp2_add(&v1, &E->A, &E->A);
    fp2_mul(&v1, &v1, &Q->z);
    fp2_add(&v2, &v2, &v1);
    fp2_mul(&v4, &P->x, &Q->x);
    fp2_add(&v4, &v4, &Q->z);
    fp2_mul(&v2, &v2, &v4);
    fp2_mul(&v1, &v1, &Q->z);
    fp2_sub(&v2, &v2, &v1);
    fp2_mul(&v2, &v2, &B->PmQ.z);
    fp2_sub(&Q->y, &v3, &v2);
    fp2_add(&v1, &P->y, &P->y);
    fp2_mul(&v1, &v1, &Q->z);
    fp2_mul(&v1, &v1, &B->PmQ.z);
    fp2_mul(&Q->x, &Q->x, &v1);
    fp2_mul(&Q->z, &Q->z, &v1);

    // Transforming to a jacobian coordinate
    fp2_sqr(&v1, &Q->z);
    fp2_mul(&Q->y, &Q->y, &v1);
    fp2_mul(&Q->x, &Q->x, &Q->z);
    return ret;
}

uint32_t
lift_basis(jac_point_t *P, jac_point_t *Q, ec_basis_t *B, ec_curve_t *E)
{
    // Normalise the curve E such that (A : C) is (A/C : 1)
    // and the point x(P) = (X/Z : 1).
    fp2_t inverses[2];
    fp2_copy(&inverses[0], &B->P.z);
    fp2_copy(&inverses[1], &E->C);

    fp2_batched_inv(inverses, 2);
    fp2_set_one(&B->P.z);
    fp2_set_one(&E->C);

    fp2_mul(&B->P.x, &B->P.x, &inverses[0]);
    fp2_mul(&E->A, &E->A, &inverses[1]);

    // Lift the basis to Jacobian points P, Q
    return lift_basis_normalized(P, Q, B, E);
}

void
DBL(jac_point_t *Q, const jac_point_t *P, const ec_curve_t *AC)
{ // Cost of 6M + 6S.
    // Doubling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to
    // (X/Z^2,Y/Z^3) This version receives the coefficient value A
    fp2_t t0, t1, t2, t3;

    uint32_t flag = fp2_is_zero(&P->x) & fp2_is_zero(&P->z);

    fp2_sqr(&t0, &P->x); // t0 = x1^2
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t0, &t0, &t1); // t0 = 3x1^2
    fp2_sqr(&t1, &P->z);    // t1 = z1^2
    fp2_mul(&t2, &P->x, &AC->A);
    fp2_add(&t2, &t2, &t2); // t2 = 2Ax1
    fp2_add(&t2, &t1, &t2); // t2 = 2Ax1+z1^2
    fp2_mul(&t2, &t1, &t2); // t2 = z1^2(2Ax1+z1^2)
    fp2_add(&t2, &t0, &t2); // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)
    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z); // z2 = 2y1z1
    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t0, &t0, &AC->A); // t0 = 4Ay1^2z1^2
    fp2_sqr(&t1, &P->y);
    fp2_add(&t1, &t1, &t1);     // t1 = 2y1^2
    fp2_add(&t3, &P->x, &P->x); // t3 = 2x1
    fp2_mul(&t3, &t1, &t3);     // t3 = 4x1y1^2
    fp2_sqr(&Q->x, &t2);        // x2 = alpha^2
    fp2_sub(&Q->x, &Q->x, &t0); // x2 = alpha^2 - 4Ay1^2z1^2
    fp2_sub(&Q->x, &Q->x, &t3);
    fp2_sub(&Q->x, &Q->x, &t3); // x2 = alpha^2 - 4Ay1^2z1^2 - 8x1y1^2
    fp2_sub(&Q->y, &t3, &Q->x); // y2 = 4x1y1^2 - x2
    fp2_mul(&Q->y, &Q->y, &t2); // y2 = alpha(4x1y1^2 - x2)
    fp2_sqr(&t1, &t1);          // t1 = 4y1^4
    fp2_sub(&Q->y, &Q->y, &t1);
    fp2_sub(&Q->y, &Q->y, &t1); // y2 = alpha(4x1y1^2 - x2) - 8y1^4

    fp2_select(&Q->x, &Q->x, &P->x, -flag);
    fp2_select(&Q->z, &Q->z, &P->z, -flag);
}

/**
 * @brief Check if a Jacobian point (X : Y : Z) has order exactly 2^f
 *
 * @param P: a point
 * @param E: an elliptic curve
 * @param t: an integer
 *
 * @return 0xFFFFFFFF if the order is correct, 0 otherwise
 */
static int
test_jac_order_twof(const jac_point_t *P, const ec_curve_t *E, int t)
{
    jac_point_t test;
    test = *P;
    if (fp2_is_zero(&test.z))
        return 0;
    for (int i = 0; i < t - 1; i++) {
        DBL(&test, &test, E);
    }
    if (fp2_is_zero(&test.z))
        return 0;
    DBL(&test, &test, E);
    return (fp2_is_zero(&test.z));
}


void
double_couple_jac_point(theta_couple_jac_point_t *out,
                        const theta_couple_jac_point_t *in,
                        const theta_couple_curve_t *E1E2)
{
    DBL(&out->P1, &in->P1, &E1E2->E1);
    DBL(&out->P2, &in->P2, &E1E2->E2);
}

void
jac_to_ws(jac_point_t *Q, fp2_t *t, fp2_t *ao3, const jac_point_t *P, const ec_curve_t *curve)
{
    // Cost of 3M + 2S when A != 0.
    fp_t one;
    fp2_t a;
    /* a = 1 - A^2/3, U = X + (A*Z^2)/3, V = Y, W = Z, T = a*Z^4*/
    fp_set_one(&one);
    if (!fp2_is_zero(&(curve->A))) {
        fp_div3(&(ao3->re), &(curve->A.re));
        fp_div3(&(ao3->im), &(curve->A.im));
        fp2_sqr(t, &P->z);
        fp2_mul(&Q->x, ao3, t);
        fp2_add(&Q->x, &Q->x, &P->x);
        fp2_sqr(t, t);
        fp2_mul(&a, ao3, &(curve->A));
        fp_sub(&(a.re), &one, &(a.re));
        fp_neg(&(a.im), &(a.im));
        fp2_mul(t, t, &a);
    } else {
        fp2_copy(&Q->x, &P->x);
        fp2_sqr(t, &P->z);
        fp2_sqr(t, t);
    }
    fp2_copy(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
}

void
jac_from_ws(jac_point_t *Q, const jac_point_t *P, const fp2_t *ao3, const ec_curve_t *curve)
{
    // Cost of 1M + 1S when A != 0.
    fp2_t t;
    /* X = U - (A*W^2)/3, Y = V, Z = W. */
    if (!fp2_is_zero(&(curve->A))) {
        fp2_sqr(&t, &P->z);
        fp2_mul(&t, &t, ao3);
        fp2_sub(&Q->x, &P->x, &t);
    }
    fp2_copy(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
}


void
DBLW(jac_point_t *Q, fp2_t *u, const jac_point_t *P, const fp2_t *t)
{ // Cost of 3M + 5S.
    // Doubling on a Weierstrass curve, representation in modified Jacobian coordinates
    // (X:Y:Z:T=a*Z^4) corresponding to (X/Z^2,Y/Z^3), where a is the curve coefficient.
    // Formula from https://hyperelliptic.org/EFD/g1p/auto-shortw-modified.html

    uint32_t flag = fp2_is_zero(&P->x) & fp2_is_zero(&P->z);

    fp2_t xx, c, cc, r, s, m;
    // XX = X^2
    fp2_sqr(&xx, &P->x);
    // A = 2*Y^2
    fp2_sqr(&c, &P->y);
    fp2_add(&c, &c, &c);
    // AA = A^2
    fp2_sqr(&cc, &c);
    // R = 2*AA
    fp2_add(&r, &cc, &cc);
    // S = (X+A)^2-XX-AA
    fp2_add(&s, &P->x, &c);
    fp2_sqr(&s, &s);
    fp2_sub(&s, &s, &xx);
    fp2_sub(&s, &s, &cc);
    // M = 3*XX+T1
    fp2_add(&m, &xx, &xx);
    fp2_add(&m, &m, &xx);
    fp2_add(&m, &m, t);
    // X3 = M^2-2*S
    fp2_sqr(&Q->x, &m);
    fp2_sub(&Q->x, &Q->x, &s);
    fp2_sub(&Q->x, &Q->x, &s);
    // Z3 = 2*Y*Z
    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z);
    // Y3 = M*(S-X3)-R
    fp2_sub(&Q->y, &s, &Q->x);
    fp2_mul(&Q->y, &Q->y, &m);
    fp2_sub(&Q->y, &Q->y, &r);
    // T3 = 2*R*T1
    fp2_mul(u, t, &r);
    fp2_add(u, u, u);

    fp2_select(&Q->x, &Q->x, &P->x, -flag);
    fp2_select(&Q->z, &Q->z, &P->z, -flag);
}

void
double_couple_jac_point_iter(theta_couple_jac_point_t *out,
                             unsigned n,
                             const theta_couple_jac_point_t *in,
                             const theta_couple_curve_t *E1E2)
{
    if (n == 0) {
        *out = *in;
    } else if (n == 1) {
        double_couple_jac_point(out, in, E1E2);
    } else {
        fp2_t a1, a2, t1, t2;

        jac_to_ws(&out->P1, &t1, &a1, &in->P1, &E1E2->E1);
        jac_to_ws(&out->P2, &t2, &a2, &in->P2, &E1E2->E2);

        DBLW(&out->P1, &t1, &out->P1, &t1);
        DBLW(&out->P2, &t2, &out->P2, &t2);
        for (unsigned i = 0; i < n - 1; i++) {
            DBLW(&out->P1, &t1, &out->P1, &t1);
            DBLW(&out->P2, &t2, &out->P2, &t2);
        }

        jac_from_ws(&out->P1, &out->P1, &a1, &E1E2->E1);
        jac_from_ws(&out->P2, &out->P2, &a2, &E1E2->E2);
    }
}

void
jac_to_xz(ec_point_t *P, const jac_point_t *xyP)
{
    fp2_copy(&P->x, &xyP->x);
    fp2_copy(&P->z, &xyP->z);
    fp2_sqr(&P->z, &P->z);

    // If xyP = (0:1:0), we currently have P=(0 : 0) but we want to set P=(1:0)
    uint32_t c1, c2;
    fp2_t one;
    fp2_set_one(&one);

    c1 = fp2_is_zero(&P->x);
    c2 = fp2_is_zero(&P->z);
    fp2_select(&P->x, &P->x, &one, c1 & c2);
}

void
couple_jac_to_xz(theta_couple_point_t *P, const theta_couple_jac_point_t *xyP)
{
    jac_to_xz(&P->P1, &xyP->P1);
    jac_to_xz(&P->P2, &xyP->P2);
}



void
double_couple_point(theta_couple_point_t *out, const theta_couple_point_t *in, const theta_couple_curve_t *E1E2)
{
    ec_dbl(&out->P1, &in->P1, &E1E2->E1);
    ec_dbl(&out->P2, &in->P2, &E1E2->E2);
}


// Returns 1 if the basis is as expected and 0 otherwise
// We only expect this to fail for malformed signatures, so
// do not require this to run in constant time.
static int
verify_two_torsion(const theta_couple_point_t *K1_2, const theta_couple_point_t *K2_2, const theta_couple_curve_t *E12)
{
    // First check if any point in K1_2 or K2_2 is zero, if they are then the points did not have
    // order 8 when we started gluing
    if (ec_is_zero(&K1_2->P1) | ec_is_zero(&K1_2->P2) | ec_is_zero(&K2_2->P1) | ec_is_zero(&K2_2->P2)) {
        return 0;
    }

    // Now ensure that P1, Q1 and P2, Q2 are independent. For points of order two this means
    // that they're not the same
    if (ec_is_equal(&K1_2->P1, &K2_2->P1) | ec_is_equal(&K1_2->P2, &K2_2->P2)) {
        return 0;
    }

    // Finally, double points to ensure all points have order exactly 0
    theta_couple_point_t O1, O2;
    double_couple_point(&O1, K1_2, E12);
    double_couple_point(&O2, K2_2, E12);
    // If this check fails then the points had order 2*f for some f, and the kernel is malformed.
    if (!(ec_is_zero(&O1.P1) & ec_is_zero(&O1.P2) & ec_is_zero(&O2.P1) & ec_is_zero(&O2.P2))) {
        return 0;
    }

    return 1;
}


static void
action_by_translation_z_and_det(fp2_t *z_inv, fp2_t *det_inv, const ec_point_t *P4, const ec_point_t *P2)
{
    // Store the Z-coordinate to invert
    fp2_copy(z_inv, &P4->z);

    // Then collect detij = xij wij - uij zij
    fp2_t tmp;
    fp2_mul(det_inv, &P4->x, &P2->z);
    fp2_mul(&tmp, &P4->z, &P2->x);
    fp2_sub(det_inv, det_inv, &tmp);
}

static void
action_by_translation_compute_matrix(translation_matrix_t *G,
                                     const ec_point_t *P4,
                                     const ec_point_t *P2,
                                     const fp2_t *z_inv,
                                     const fp2_t *det_inv)
{
    fp2_t tmp;

    // Gi.g10 = uij xij /detij - xij/zij
    fp2_mul(&tmp, &P4->x, z_inv);
    fp2_mul(&G->g10, &P4->x, &P2->x);
    fp2_mul(&G->g10, &G->g10, det_inv);
    fp2_sub(&G->g10, &G->g10, &tmp);

    // Gi.g11 = uij zij * detij
    fp2_mul(&G->g11, &P2->x, det_inv);
    fp2_mul(&G->g11, &G->g11, &P4->z);

    // Gi.g00 = -Gi.g11
    fp2_neg(&G->g00, &G->g11);

    // Gi.g01 = - wij zij detij
    fp2_mul(&G->g01, &P2->z, det_inv);
    fp2_mul(&G->g01, &G->g01, &P4->z);
    fp2_neg(&G->g01, &G->g01);
}


// Computes the action by translation for four points
// (P1, P2) and (Q1, Q2) on E1 x E2 simultaneously to
// save on inversions.
// Returns 0 if any of Pi or Qi does not have order 2
// and 1 otherwise
static int
action_by_translation(translation_matrix_t *Gi,
                      const theta_couple_point_t *K1_4,
                      const theta_couple_point_t *K2_4,
                      const theta_couple_curve_t *E12)
{
    // Compute points of order 2 from Ki_4
    theta_couple_point_t K1_2, K2_2;
    double_couple_point(&K1_2, K1_4, E12);
    double_couple_point(&K2_2, K2_4, E12);

    if (!verify_two_torsion(&K1_2, &K2_2, E12)) {
        return 0;
    }

    // We need to invert four Z coordinates and
    // four determinants which we do with batched
    // inversion
    fp2_t inverses[8];
    action_by_translation_z_and_det(&inverses[0], &inverses[4], &K1_4->P1, &K1_2.P1);
    action_by_translation_z_and_det(&inverses[1], &inverses[5], &K1_4->P2, &K1_2.P2);
    action_by_translation_z_and_det(&inverses[2], &inverses[6], &K2_4->P1, &K2_2.P1);
    action_by_translation_z_and_det(&inverses[3], &inverses[7], &K2_4->P2, &K2_2.P2);

    fp2_batched_inv(inverses, 8);
    if (fp2_is_zero(&inverses[0]))
        return 0; // something was wrong with our input (which somehow was not caught by
    // verify_two_torsion)

    action_by_translation_compute_matrix(&Gi[0], &K1_4->P1, &K1_2.P1, &inverses[0], &inverses[4]);
    action_by_translation_compute_matrix(&Gi[1], &K1_4->P2, &K1_2.P2, &inverses[1], &inverses[5]);
    action_by_translation_compute_matrix(&Gi[2], &K2_4->P1, &K2_2.P1, &inverses[2], &inverses[6]);
    action_by_translation_compute_matrix(&Gi[3], &K2_4->P2, &K2_2.P2, &inverses[3], &inverses[7]);

    return 1;
}

// Given the appropriate four torsion, computes the
// change of basis to compute the correct theta null
// point.
// Returns 0 if the order of K1_4 or K2_4 is not 4
static int
gluing_change_of_basis(basis_change_matrix_t *M,
                       const theta_couple_point_t *K1_4,
                       const theta_couple_point_t *K2_4,
                       const theta_couple_curve_t *E12)
{
    // Compute the four 2x2 matrices for the action by translation
    // on the four points:
    translation_matrix_t Gi[4];
    if (!action_by_translation(Gi, K1_4, K2_4, E12))
        return 0;

    // Computation of the 4x4 matrix from Mij
    // t001, t101 (resp t002, t102) first column of M11 * M21 (resp M12 * M22)
    fp2_t t001, t101, t002, t102, tmp;

    fp2_mul(&t001, &Gi[0].g00, &Gi[2].g00);
    fp2_mul(&tmp, &Gi[0].g01, &Gi[2].g10);
    fp2_add(&t001, &t001, &tmp);

    fp2_mul(&t101, &Gi[0].g10, &Gi[2].g00);
    fp2_mul(&tmp, &Gi[0].g11, &Gi[2].g10);
    fp2_add(&t101, &t101, &tmp);

    fp2_mul(&t002, &Gi[1].g00, &Gi[3].g00);
    fp2_mul(&tmp, &Gi[1].g01, &Gi[3].g10);
    fp2_add(&t002, &t002, &tmp);

    fp2_mul(&t102, &Gi[1].g10, &Gi[3].g00);
    fp2_mul(&tmp, &Gi[1].g11, &Gi[3].g10);
    fp2_add(&t102, &t102, &tmp);

    // trace for the first row
    fp2_set_one(&M->m[0][0]);
    fp2_mul(&tmp, &t001, &t002);
    fp2_add(&M->m[0][0], &M->m[0][0], &tmp);
    fp2_mul(&tmp, &Gi[2].g00, &Gi[3].g00);
    fp2_add(&M->m[0][0], &M->m[0][0], &tmp);
    fp2_mul(&tmp, &Gi[0].g00, &Gi[1].g00);
    fp2_add(&M->m[0][0], &M->m[0][0], &tmp);

    fp2_mul(&M->m[0][1], &t001, &t102);
    fp2_mul(&tmp, &Gi[2].g00, &Gi[3].g10);
    fp2_add(&M->m[0][1], &M->m[0][1], &tmp);
    fp2_mul(&tmp, &Gi[0].g00, &Gi[1].g10);
    fp2_add(&M->m[0][1], &M->m[0][1], &tmp);

    fp2_mul(&M->m[0][2], &t101, &t002);
    fp2_mul(&tmp, &Gi[2].g10, &Gi[3].g00);
    fp2_add(&M->m[0][2], &M->m[0][2], &tmp);
    fp2_mul(&tmp, &Gi[0].g10, &Gi[1].g00);
    fp2_add(&M->m[0][2], &M->m[0][2], &tmp);

    fp2_mul(&M->m[0][3], &t101, &t102);
    fp2_mul(&tmp, &Gi[2].g10, &Gi[3].g10);
    fp2_add(&M->m[0][3], &M->m[0][3], &tmp);
    fp2_mul(&tmp, &Gi[0].g10, &Gi[1].g10);
    fp2_add(&M->m[0][3], &M->m[0][3], &tmp);

    // Compute the action of (0,out.K2_4.P2) for the second row
    fp2_mul(&tmp, &Gi[3].g01, &M->m[0][1]);
    fp2_mul(&M->m[1][0], &Gi[3].g00, &M->m[0][0]);
    fp2_add(&M->m[1][0], &M->m[1][0], &tmp);

    fp2_mul(&tmp, &Gi[3].g11, &M->m[0][1]);
    fp2_mul(&M->m[1][1], &Gi[3].g10, &M->m[0][0]);
    fp2_add(&M->m[1][1], &M->m[1][1], &tmp);

    fp2_mul(&tmp, &Gi[3].g01, &M->m[0][3]);
    fp2_mul(&M->m[1][2], &Gi[3].g00, &M->m[0][2]);
    fp2_add(&M->m[1][2], &M->m[1][2], &tmp);

    fp2_mul(&tmp, &Gi[3].g11, &M->m[0][3]);
    fp2_mul(&M->m[1][3], &Gi[3].g10, &M->m[0][2]);
    fp2_add(&M->m[1][3], &M->m[1][3], &tmp);

    // compute the action of (K1_4.P1,0) for the third row
    fp2_mul(&tmp, &Gi[0].g01, &M->m[0][2]);
    fp2_mul(&M->m[2][0], &Gi[0].g00, &M->m[0][0]);
    fp2_add(&M->m[2][0], &M->m[2][0], &tmp);

    fp2_mul(&tmp, &Gi[0].g01, &M->m[0][3]);
    fp2_mul(&M->m[2][1], &Gi[0].g00, &M->m[0][1]);
    fp2_add(&M->m[2][1], &M->m[2][1], &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m[0][2]);
    fp2_mul(&M->m[2][2], &Gi[0].g10, &M->m[0][0]);
    fp2_add(&M->m[2][2], &M->m[2][2], &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m[0][3]);
    fp2_mul(&M->m[2][3], &Gi[0].g10, &M->m[0][1]);
    fp2_add(&M->m[2][3], &M->m[2][3], &tmp);

    // compute the action of (K1_4.P1,K2_4.P2) for the final row
    fp2_mul(&tmp, &Gi[0].g01, &M->m[1][2]);
    fp2_mul(&M->m[3][0], &Gi[0].g00, &M->m[1][0]);
    fp2_add(&M->m[3][0], &M->m[3][0], &tmp);

    fp2_mul(&tmp, &Gi[0].g01, &M->m[1][3]);
    fp2_mul(&M->m[3][1], &Gi[0].g00, &M->m[1][1]);
    fp2_add(&M->m[3][1], &M->m[3][1], &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m[1][2]);
    fp2_mul(&M->m[3][2], &Gi[0].g10, &M->m[1][0]);
    fp2_add(&M->m[3][2], &M->m[3][2], &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m[1][3]);
    fp2_mul(&M->m[3][3], &Gi[0].g10, &M->m[1][1]);
    fp2_add(&M->m[3][3], &M->m[3][3], &tmp);

    return 1;
}

// same as apply_isomorphism method but more efficient when the t component of P is zero.
static void
apply_isomorphism_general(theta_point_t *res,
                          const basis_change_matrix_t *M,
                          const theta_point_t *P,
                          const bool Pt_not_zero)
{
    fp2_t x1;
    theta_point_t temp;

    fp2_mul(&temp.x, &P->x, &M->m[0][0]);
    fp2_mul(&x1, &P->y, &M->m[0][1]);
    fp2_add(&temp.x, &temp.x, &x1);
    fp2_mul(&x1, &P->z, &M->m[0][2]);
    fp2_add(&temp.x, &temp.x, &x1);

    fp2_mul(&temp.y, &P->x, &M->m[1][0]);
    fp2_mul(&x1, &P->y, &M->m[1][1]);
    fp2_add(&temp.y, &temp.y, &x1);
    fp2_mul(&x1, &P->z, &M->m[1][2]);
    fp2_add(&temp.y, &temp.y, &x1);

    fp2_mul(&temp.z, &P->x, &M->m[2][0]);
    fp2_mul(&x1, &P->y, &M->m[2][1]);
    fp2_add(&temp.z, &temp.z, &x1);
    fp2_mul(&x1, &P->z, &M->m[2][2]);
    fp2_add(&temp.z, &temp.z, &x1);

    fp2_mul(&temp.t, &P->x, &M->m[3][0]);
    fp2_mul(&x1, &P->y, &M->m[3][1]);
    fp2_add(&temp.t, &temp.t, &x1);
    fp2_mul(&x1, &P->z, &M->m[3][2]);
    fp2_add(&temp.t, &temp.t, &x1);

    if (Pt_not_zero) {
        fp2_mul(&x1, &P->t, &M->m[0][3]);
        fp2_add(&temp.x, &temp.x, &x1);

        fp2_mul(&x1, &P->t, &M->m[1][3]);
        fp2_add(&temp.y, &temp.y, &x1);

        fp2_mul(&x1, &P->t, &M->m[2][3]);
        fp2_add(&temp.z, &temp.z, &x1);

        fp2_mul(&x1, &P->t, &M->m[3][3]);
        fp2_add(&temp.t, &temp.t, &x1);
    }

    fp2_copy(&res->x, &temp.x);
    fp2_copy(&res->y, &temp.y);
    fp2_copy(&res->z, &temp.z);
    fp2_copy(&res->t, &temp.t);
}

static void
apply_isomorphism(theta_point_t *res, const basis_change_matrix_t *M, const theta_point_t *P)
{
    apply_isomorphism_general(res, M, P, true);
}

// compute the theta_point corresponding to the couple of point T on an elliptic product
static void
base_change(theta_point_t *out, const theta_gluing_t *phi, const theta_couple_point_t *T)
{
    theta_point_t null_point;

    // null_point = (a : b : c : d)
    // a = P1.x P2.x, b = P1.x P2.z, c = P1.z P2.x, d = P1.z P2.z
    fp2_mul(&null_point.x, &T->P1.x, &T->P2.x);
    fp2_mul(&null_point.y, &T->P1.x, &T->P2.z);
    fp2_mul(&null_point.z, &T->P2.x, &T->P1.z);
    fp2_mul(&null_point.t, &T->P1.z, &T->P2.z);

    // Apply the basis change
    apply_isomorphism(out, &phi->M, &null_point);
}
/**
 * @brief Square the coordinates of a theta point
 * @param out Output: the theta_point
 * @param in a theta point*
 * in = (x,y,z,t)
 * out = (x^2, y^2, z^2, t^2)
 *
 */
static inline void
pointwise_square(theta_point_t *out, const theta_point_t *in)
{
    fp2_sqr(&out->x, &in->x);
    fp2_sqr(&out->y, &in->y);
    fp2_sqr(&out->z, &in->z);
    fp2_sqr(&out->t, &in->t);
}

/**
 * @brief Perform the hadamard transform on a theta point
 *
 * @param out Output: the theta_point
 * @param in a theta point*
 * in = (x,y,z,t)
 * out = (x+y+z+t, x-y+z-t, x+y-z-t, x-y-z+t)
 *
 */
static inline void
hadamard(theta_point_t *out, const theta_point_t *in)
{
    fp2_t t1, t2, t3, t4;

    // t1 = x + y
    fp2_add(&t1, &in->x, &in->y);
    // t2 = x - y
    fp2_sub(&t2, &in->x, &in->y);
    // t3 = z + t
    fp2_add(&t3, &in->z, &in->t);
    // t4 = z - t
    fp2_sub(&t4, &in->z, &in->t);

    fp2_add(&out->x, &t1, &t3);
    fp2_add(&out->y, &t2, &t4);
    fp2_sub(&out->z, &t1, &t3);
    fp2_sub(&out->t, &t2, &t4);
}


/**
 * @brief Square the coordinates and then perform the hadamard transform
 *
 * @param out Output: the theta_point
 * @param in a theta point*
 * in = (x,y,z,t)
 * out = (x^2+y^2+z^2+t^2, x^2-y^2+z^2-t^2, x^2+y^2-z^2-t^2, x^2-y^2-z^2+t^2)
 *
 */
static inline void
to_squared_theta(theta_point_t *out, const theta_point_t *in)
{
    pointwise_square(out, in);
    hadamard(out, out);
}


/**
 * @brief Compute the gluing isogeny from an elliptic product
 *
 * @param out Output: the theta_gluing
 * @param K1_8 a couple point
 * @param E12 an elliptic curve product
 * @param K2_8 a point in E2[8]
 *
 * out : E1xE2 -> A of kernel [4](K1_8,K2_8)
 * if the kernel supplied has the incorrect order, or gluing seems malformed,
 * returns 0, otherwise returns 1.
 */
static int
gluing_compute(theta_gluing_t *out,
               const theta_couple_curve_t *E12,
               const theta_couple_jac_point_t *xyK1_8,
               const theta_couple_jac_point_t *xyK2_8,
               bool verify)
{
    // Ensure that we have been given the eight torsion
#ifndef NDEBUG
    {
        int check = test_jac_order_twof(&xyK1_8->P1, &E12->E1, 3);
        if (!check)
            debug_print("xyK1_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK2_8->P1, &E12->E1, 3);
        if (!check)
            debug_print("xyK2_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK1_8->P2, &E12->E2, 3);
        if (!check)
            debug_print("xyK2_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK2_8->P2, &E12->E2, 3);
        if (!check)
            debug_print("xyK2_8->P2 does not have order 8");
    }
#endif

    out->xyK1_8 = *xyK1_8;
    out->domain = *E12;

    // Given points in E[8] x E[8] we need the four torsion below
    theta_couple_jac_point_t xyK1_4, xyK2_4;

    double_couple_jac_point(&xyK1_4, xyK1_8, E12);
    double_couple_jac_point(&xyK2_4, xyK2_8, E12);

    // Convert from (X:Y:Z) coordinates to (X:Z)
    theta_couple_point_t K1_8, K2_8;
    theta_couple_point_t K1_4, K2_4;

    couple_jac_to_xz(&K1_8, xyK1_8);
    couple_jac_to_xz(&K2_8, xyK2_8);
    couple_jac_to_xz(&K1_4, &xyK1_4);
    couple_jac_to_xz(&K2_4, &xyK2_4);

    // Set the basis change matrix, if we have not been given a valid K[8] for this computation
    // gluing_change_of_basis will detect this and return 0
    if (!gluing_change_of_basis(&out->M, &K1_4, &K2_4, E12)) {
        debug_print("gluing failed as kernel does not have correct order");
        return 0;
    }

    // apply the base change to the kernel
    theta_point_t TT1, TT2;

    base_change(&TT1, out, &K1_8);
    base_change(&TT2, out, &K2_8);

    // compute the codomain
    to_squared_theta(&TT1, &TT1);
    to_squared_theta(&TT2, &TT2);

    // If the kernel is well formed then TT1.t and TT2.t are zero
    // if they are not, we exit early as the signature we are validating
    // is probably malformed
    if (!(fp2_is_zero(&TT1.t) & fp2_is_zero(&TT2.t))) {
        debug_print("gluing failed TT1.t or TT2.t is not zero");
        return 0;
    }
    // Test our projective factors are non zero
    if (fp2_is_zero(&TT1.x) | fp2_is_zero(&TT2.x) | fp2_is_zero(&TT1.y) | fp2_is_zero(&TT2.z) | fp2_is_zero(&TT1.z))
        return 0; // invalid input

    // Projective factor: Ax
    fp2_mul(&out->codomain.x, &TT1.x, &TT2.x);
    fp2_mul(&out->codomain.y, &TT1.y, &TT2.x);
    fp2_mul(&out->codomain.z, &TT1.x, &TT2.z);
    fp2_set_zero(&out->codomain.t);
    // Projective factor: ABCxz
    fp2_mul(&out->precomputation.x, &TT1.y, &TT2.z);
    fp2_copy(&out->precomputation.y, &out->codomain.z);
    fp2_copy(&out->precomputation.z, &out->codomain.y);
    fp2_set_zero(&out->precomputation.t);

    // Compute the two components of phi(K1_8) = (x:x:y:y).
    fp2_mul(&out->imageK1_8.x, &TT1.x, &out->precomputation.x);
    fp2_mul(&out->imageK1_8.y, &TT1.z, &out->precomputation.z);

    // If K1_8 and K2_8 are our 8-torsion points, this ensures that the
    // 4-torsion points [2]K1_8 and [2]K2_8 are isotropic.
    if (verify) {
        fp2_t t1, t2;
        fp2_mul(&t1, &TT1.y, &out->precomputation.y);
        if (!fp2_is_equal(&out->imageK1_8.x, &t1))
            return 0;
        fp2_mul(&t1, &TT2.x, &out->precomputation.x);
        fp2_mul(&t2, &TT2.z, &out->precomputation.z);
        if (!fp2_is_equal(&t2, &t1))
            return 0;
    }

    // compute the final codomain
    hadamard(&out->codomain, &out->codomain);
    return 1;
}



// Same as gluing_eval_point but in the very special case where we already know that the point will
// have a zero coordinate at the place where the zero coordinate of the dual_theta_nullpoint would
// have made the computation difficult
static int
gluing_eval_point_special_case(theta_point_t *image, const theta_couple_point_t *P, const theta_gluing_t *phi)
{
    theta_point_t T;

    // Apply the basis change
    base_change(&T, phi, P);

    // Apply the to_squared_theta transform
    to_squared_theta(&T, &T);

    // This coordinate should always be 0 in a gluing because D=0.
    // If this is not the case, something went very wrong, so reject
    if (!fp2_is_zero(&T.t))
        return 0;

    // Compute (x, y, z, t)
    fp2_mul(&image->x, &T.x, &phi->precomputation.x);
    fp2_mul(&image->y, &T.y, &phi->precomputation.y);
    fp2_mul(&image->z, &T.z, &phi->precomputation.z);
    fp2_set_zero(&image->t);

    hadamard(image, image);
    return 1;
}


void
jac_to_xz_add_components(add_components_t *add_comp, const jac_point_t *P, const jac_point_t *Q, const ec_curve_t *AC)
{
    // Take P and Q in E distinct, two jac_point_t, return three components u,v and w in Fp2 such
    // that the xz coordinates of P+Q are (u-v:w) and of P-Q are (u+v:w)

    fp2_t t0, t1, t2, t3, t4, t5, t6;

    fp2_sqr(&t0, &P->z);             // t0 = z1^2
    fp2_sqr(&t1, &Q->z);             // t1 = z2^2
    fp2_mul(&t2, &P->x, &t1);        // t2 = x1z2^2
    fp2_mul(&t3, &t0, &Q->x);        // t3 = z1^2x2
    fp2_mul(&t4, &P->y, &Q->z);      // t4 = y1z2
    fp2_mul(&t4, &t4, &t1);          // t4 = y1z2^3
    fp2_mul(&t5, &P->z, &Q->y);      // t5 = z1y2
    fp2_mul(&t5, &t5, &t0);          // t5 = z1^3y2
    fp2_mul(&t0, &t0, &t1);          // t0 = (z1z2)^2
    fp2_mul(&t6, &t4, &t5);          // t6 = (z1z_2)^3y1y2
    fp2_add(&add_comp->v, &t6, &t6); // v  = 2(z1z_2)^3y1y2
    fp2_sqr(&t4, &t4);               // t4 = y1^2z2^6
    fp2_sqr(&t5, &t5);               // t5 = z1^6y_2^2
    fp2_add(&t4, &t4, &t5);          // t4 = z1^6y_2^2 + y1^2z2^6
    fp2_add(&t5, &t2, &t3);          // t5 = x1z2^2 +z_1^2x2
    fp2_add(&t6, &t3, &t3);          // t6 = 2z_1^2x2
    fp2_sub(&t6, &t5, &t6);          // t6 = lambda = x1z2^2 - z_1^2x2
    fp2_sqr(&t6, &t6);               // t6 = lambda^2 = (x1z2^2 - z_1^2x2)^2
    fp2_mul(&t1, &AC->A, &t0);       // t1 = A*(z1z2)^2
    fp2_add(&t1, &t5, &t1);          // t1 = gamma =A*(z1z2)^2 + x1z2^2 +z_1^2x2
    fp2_mul(&t1, &t1, &t6);          // t1 = gamma*lambda^2
    fp2_sub(&add_comp->u, &t4, &t1); // u  = z1^6y_2^2 + y1^2z2^6 - gamma*lambda^2
    fp2_mul(&add_comp->w, &t6, &t0); // w  = (z1z2)^2(lambda)^2
}

// sub routine of the gluing eval
static void
gluing_eval_point(theta_point_t *image, const theta_couple_jac_point_t *P, const theta_gluing_t *phi)
{
    theta_point_t T1, T2;
    add_components_t add_comp1, add_comp2;

    // Compute the cross addition components of P1+Q1 and P2+Q2
    jac_to_xz_add_components(&add_comp1, &P->P1, &phi->xyK1_8.P1, &phi->domain.E1);
    jac_to_xz_add_components(&add_comp2, &P->P2, &phi->xyK1_8.P2, &phi->domain.E2);

    // Compute T1 and T2 derived from the cross addition components.
    fp2_mul(&T1.x, &add_comp1.u, &add_comp2.u); // T1x = u1u2
    fp2_mul(&T2.t, &add_comp1.v, &add_comp2.v); // T2t = v1v2
    fp2_add(&T1.x, &T1.x, &T2.t);               // T1x = u1u2 + v1v2
    fp2_mul(&T1.y, &add_comp1.u, &add_comp2.w); // T1y = u1w2
    fp2_mul(&T1.z, &add_comp1.w, &add_comp2.u); // T1z = w1u2
    fp2_mul(&T1.t, &add_comp1.w, &add_comp2.w); // T1t = w1w2
    fp2_add(&T2.x, &add_comp1.u, &add_comp1.v); // T2x = (u1+v1)
    fp2_add(&T2.y, &add_comp2.u, &add_comp2.v); // T2y = (u2+v2)
    fp2_mul(&T2.x, &T2.x, &T2.y);               // T2x = (u1+v1)(u2+v2)
    fp2_sub(&T2.x, &T2.x, &T1.x);               // T1x = v1u2 + u1v2
    fp2_mul(&T2.y, &add_comp1.v, &add_comp2.w); // T2y = v1w2
    fp2_mul(&T2.z, &add_comp1.w, &add_comp2.v); // T2z = w1v2
    fp2_set_zero(&T2.t);                        // T2t = 0

    // Apply the basis change and compute their respective square
    // theta(P+Q) = M.T1 - M.T2 and theta(P-Q) = M.T1 + M.T2
    apply_isomorphism_general(&T1, &phi->M, &T1, true);
    apply_isomorphism_general(&T2, &phi->M, &T2, false);
    pointwise_square(&T1, &T1);
    pointwise_square(&T2, &T2);

    // the difference between the two is therefore theta(P+Q)theta(P-Q)
    // whose hadamard transform is then the product of the dual
    // theta_points of phi(P) and phi(Q).
    fp2_sub(&T1.x, &T1.x, &T2.x);
    fp2_sub(&T1.y, &T1.y, &T2.y);
    fp2_sub(&T1.z, &T1.z, &T2.z);
    fp2_sub(&T1.t, &T1.t, &T2.t);
    hadamard(&T1, &T1);

    // Compute (x, y, z, t)
    // As imageK1_8 = (x:x:y:y), its inverse is (y:y:x:x).
    fp2_mul(&image->x, &T1.x, &phi->imageK1_8.y);
    fp2_mul(&image->y, &T1.y, &phi->imageK1_8.y);
    fp2_mul(&image->z, &T1.z, &phi->imageK1_8.x);
    fp2_mul(&image->t, &T1.t, &phi->imageK1_8.x);

    hadamard(image, image);
}

/**
 * @brief Evaluate a gluing isogeny from an elliptic product on a basis
 *
 * @param image1 Output: the theta_point of the image of the first couple of points
 * @param image2 Output : the theta point of the image of the second couple of points
 * @param xyT1: A pair of points (X : Y : Z) on E1E2 to glue using phi
 * @param xyT2: A pair of points (X : Y : Z) on E1E2 to glue using phi
 * @param phi : a gluing isogeny E1 x E2 -> A
 *
 **/
static void
gluing_eval_basis(theta_point_t *image1,
                  theta_point_t *image2,
                  const theta_couple_jac_point_t *xyT1,
                  const theta_couple_jac_point_t *xyT2,
                  const theta_gluing_t *phi)
{
    gluing_eval_point(image1, xyT1, phi);
    gluing_eval_point(image2, xyT2, phi);
}


void
theta_precomputation(theta_structure_t *A)
{

    if (A->precomputation) {
        return;
    }

    theta_point_t A_dual;
    to_squared_theta(&A_dual, &A->null_point);

    fp2_t t1, t2;
    fp2_mul(&t1, &A_dual.x, &A_dual.y);
    fp2_mul(&t2, &A_dual.z, &A_dual.t);
    fp2_mul(&A->XYZ0, &t1, &A_dual.z);
    fp2_mul(&A->XYT0, &t1, &A_dual.t);
    fp2_mul(&A->YZT0, &t2, &A_dual.y);
    fp2_mul(&A->XZT0, &t2, &A_dual.x);

    fp2_mul(&t1, &A->null_point.x, &A->null_point.y);
    fp2_mul(&t2, &A->null_point.z, &A->null_point.t);
    fp2_mul(&A->xyz0, &t1, &A->null_point.z);
    fp2_mul(&A->xyt0, &t1, &A->null_point.t);
    fp2_mul(&A->yzt0, &t2, &A->null_point.y);
    fp2_mul(&A->xzt0, &t2, &A->null_point.x);

    A->precomputation = true;
}


void
double_point(theta_point_t *out, theta_structure_t *A, const theta_point_t *in)
{
    to_squared_theta(out, in);
    fp2_sqr(&out->x, &out->x);
    fp2_sqr(&out->y, &out->y);
    fp2_sqr(&out->z, &out->z);
    fp2_sqr(&out->t, &out->t);

    if (!A->precomputation) {
        theta_precomputation(A);
    }
    fp2_mul(&out->x, &out->x, &A->YZT0);
    fp2_mul(&out->y, &out->y, &A->XZT0);
    fp2_mul(&out->z, &out->z, &A->XYT0);
    fp2_mul(&out->t, &out->t, &A->XYZ0);

    hadamard(out, out);

    fp2_mul(&out->x, &out->x, &A->yzt0);
    fp2_mul(&out->y, &out->y, &A->xzt0);
    fp2_mul(&out->z, &out->z, &A->xyt0);
    fp2_mul(&out->t, &out->t, &A->xyz0);
}

void
double_iter(theta_point_t *out, theta_structure_t *A, const theta_point_t *in, int exp)
{
    if (exp == 0) {
        *out = *in;
    } else {
        double_point(out, A, in);
        for (int i = 1; i < exp; i++) {
            double_point(out, A, out);
        }
    }
}



/**
 * @brief Compute a (2,2) isogeny in dimension 2 in the theta_model
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_8 a point in A[8]
 * @param T2_8 a point in A[8]
 * @param hadamard_bool_1 a boolean used for the last two steps of the chain
 * @param hadamard_bool_2 a boolean used for the last two steps of the chain
 *
 * out : A -> B of kernel [4](T1_8,T2_8)
 * hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 * verify: add extra sanity check to ensure our 8-torsion points are coherent with the isogeny
 *
 */
static int
theta_isogeny_compute(theta_isogeny_t *out,
                      const theta_structure_t *A,
                      const theta_point_t *T1_8,
                      const theta_point_t *T2_8,
                      bool hadamard_bool_1,
                      bool hadamard_bool_2,
                      bool verify)
{
    out->hadamard_bool_1 = hadamard_bool_1;
    out->hadamard_bool_2 = hadamard_bool_2;
    out->domain = *A;
    out->T1_8 = *T1_8;
    out->T2_8 = *T2_8;
    out->codomain.precomputation = false;

    theta_point_t TT1, TT2;

    if (hadamard_bool_1) {
        hadamard(&TT1, T1_8);
        to_squared_theta(&TT1, &TT1);
        hadamard(&TT2, T2_8);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT1, T1_8);
        to_squared_theta(&TT2, T2_8);
    }

    fp2_t t1, t2;

    // Test that our projective factor ABCDxzw is non zero, where
    // TT1=(Ax, Bx, Cy, Dy), TT2=(Az, Bw, Cz, Dw)
    // But ABCDxzw=0 can only happen if we had an unexpected splitting in
    // the isogeny chain.
    // In either case reject
    // (this is not strictly necessary, we could just return (0:0:0:0))
    if (fp2_is_zero(&TT2.x) | fp2_is_zero(&TT2.y) | fp2_is_zero(&TT2.z) | fp2_is_zero(&TT2.t) | fp2_is_zero(&TT1.x) |
        fp2_is_zero(&TT1.y))
        return 0;

    fp2_mul(&t1, &TT1.x, &TT2.y);
    fp2_mul(&t2, &TT1.y, &TT2.x);
    fp2_mul(&out->codomain.null_point.x, &TT2.x, &t1);
    fp2_mul(&out->codomain.null_point.y, &TT2.y, &t2);
    fp2_mul(&out->codomain.null_point.z, &TT2.z, &t1);
    fp2_mul(&out->codomain.null_point.t, &TT2.t, &t2);
    fp2_t t3;
    fp2_mul(&t3, &TT2.z, &TT2.t);
    fp2_mul(&out->precomputation.x, &t3, &TT1.y);
    fp2_mul(&out->precomputation.y, &t3, &TT1.x);
    fp2_copy(&out->precomputation.z, &out->codomain.null_point.t);
    fp2_copy(&out->precomputation.t, &out->codomain.null_point.z);

    // If T1_8 and T2_8 are our 8-torsion points, this ensures that the
    // 4-torsion points 2T1_8 and 2T2_8 are isotropic.
    if (verify) {
        fp2_mul(&t1, &TT1.x, &out->precomputation.x);
        fp2_mul(&t2, &TT1.y, &out->precomputation.y);
        if (!fp2_is_equal(&t1, &t2))
            return 0;
        fp2_mul(&t1, &TT1.z, &out->precomputation.z);
        fp2_mul(&t2, &TT1.t, &out->precomputation.t);
        if (!fp2_is_equal(&t1, &t2))
            return 0;
        fp2_mul(&t1, &TT2.x, &out->precomputation.x);
        fp2_mul(&t2, &TT2.z, &out->precomputation.z);
        if (!fp2_is_equal(&t1, &t2))
            return 0;
        fp2_mul(&t1, &TT2.y, &out->precomputation.y);
        fp2_mul(&t2, &TT2.t, &out->precomputation.t);
        if (!fp2_is_equal(&t1, &t2))
            return 0;
    }

    if (hadamard_bool_2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
    return 1;
}

static void
theta_isogeny_eval(theta_point_t *out, const theta_isogeny_t *phi, const theta_point_t *P)
{
    if (phi->hadamard_bool_1) {
        hadamard(out, P);
        to_squared_theta(out, out);
    } else {
        to_squared_theta(out, P);
    }
    fp2_mul(&out->x, &out->x, &phi->precomputation.x);
    fp2_mul(&out->y, &out->y, &phi->precomputation.y);
    fp2_mul(&out->z, &out->z, &phi->precomputation.z);
    fp2_mul(&out->t, &out->t, &phi->precomputation.t);

    if (phi->hadamard_bool_2) {
        hadamard(out, out);
    }
}



/**
 * @brief Compute a (2,2) isogeny when only the 4 torsion above the kernel is known and not the 8
 * torsion
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_4 a point in A[4]
 * @param T2_4 a point in A[4]
 * @param hadamard_bool_1 a boolean
 * @param hadamard_bool_2 a boolean
 *
 * out : A -> B of kernel [2](T1_4,T2_4)
 * hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 *
 */
static void
theta_isogeny_compute_4(theta_isogeny_t *out,
                        const theta_structure_t *A,
                        const theta_point_t *T1_4,
                        const theta_point_t *T2_4,
                        bool hadamard_bool_1,
                        bool hadamard_bool_2)
{
    out->hadamard_bool_1 = hadamard_bool_1;
    out->hadamard_bool_2 = hadamard_bool_2;
    out->domain = *A;
    out->T1_8 = *T1_4;
    out->T2_8 = *T2_4;
    out->codomain.precomputation = false;

    theta_point_t TT1, TT2;
    // we will compute:
    // TT1 = (xAB, _ , xCD, _)
    // TT2 = (AA,BB,CC,DD)

    // fp2_t xA_inv,zA_inv,tB_inv;

    if (hadamard_bool_1) {
        hadamard(&TT1, T1_4);
        to_squared_theta(&TT1, &TT1);

        hadamard(&TT2, &A->null_point);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT1, T1_4);
        to_squared_theta(&TT2, &A->null_point);
    }

    fp2_t sqaabb, sqaacc;
    fp2_mul(&sqaabb, &TT2.x, &TT2.y);
    fp2_mul(&sqaacc, &TT2.x, &TT2.z);
    // No need to check the square roots, only used for signing.
    // sqaabb = sqrt(AA*BB)
    fp2_sqrt(&sqaabb);
    // sqaacc = sqrt(AA*CC)
    fp2_sqrt(&sqaacc);

    // we compute out->codomain.null_point = (xAB * sqaacc * AA, xAB *sqaabb *sqaacc, xCD*sqaabb *
    // AA) out->precomputation = (xAB * BB * CC *DD , sqaabb * CC * DD * xAB , sqaacc * BB* DD * xAB
    // , xCD * sqaabb *sqaacc * BB)

    fp2_mul(&out->codomain.null_point.y, &sqaabb, &sqaacc);
    fp2_mul(&out->precomputation.t, &out->codomain.null_point.y, &TT1.z);
    fp2_mul(&out->codomain.null_point.y, &out->codomain.null_point.y,
            &TT1.x); // done for out->codomain.null_point.y

    fp2_mul(&out->codomain.null_point.t, &TT1.z, &sqaabb);
    fp2_mul(&out->codomain.null_point.t, &out->codomain.null_point.t,
            &TT2.x); // done for out->codomain.null_point.t

    fp2_mul(&out->codomain.null_point.x, &TT1.x, &TT2.x);
    fp2_mul(&out->codomain.null_point.z, &out->codomain.null_point.x,
            &TT2.z); // done for out->codomain.null_point.z
    fp2_mul(&out->codomain.null_point.x, &out->codomain.null_point.x,
            &sqaacc); // done for out->codomain.null_point.x

    fp2_mul(&out->precomputation.x, &TT1.x, &TT2.t);
    fp2_mul(&out->precomputation.z, &out->precomputation.x, &TT2.y);
    fp2_mul(&out->precomputation.x, &out->precomputation.x, &TT2.z);
    fp2_mul(&out->precomputation.y, &out->precomputation.x, &sqaabb); // done for out->precomputation.y
    fp2_mul(&out->precomputation.x, &out->precomputation.x, &TT2.y);  // done for out->precomputation.x
    fp2_mul(&out->precomputation.z, &out->precomputation.z, &sqaacc); // done for out->precomputation.z
    fp2_mul(&out->precomputation.t, &out->precomputation.t, &TT2.y);  // done for out->precomputation.t

    if (hadamard_bool_2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}


/*
 * @brief Compute a (2,2) isogeny when only the kernel is known and not the 8 or 4 torsion above
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_2 a point in A[2]
 * @param T2_2 a point in A[2]
 * @param hadamard_bool_1 a boolean
 * @param boo2 a boolean
 *
 * out : A -> B of kernel (T1_2,T2_2)
 * hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 *
 */
static void
theta_isogeny_compute_2(theta_isogeny_t *out,
                        const theta_structure_t *A,
                        const theta_point_t *T1_2,
                        const theta_point_t *T2_2,
                        bool hadamard_bool_1,
                        bool hadamard_bool_2)
{
    out->hadamard_bool_1 = hadamard_bool_1;
    out->hadamard_bool_2 = hadamard_bool_2;
    out->domain = *A;
    out->T1_8 = *T1_2;
    out->T2_8 = *T2_2;
    out->codomain.precomputation = false;

    theta_point_t TT2;
    // we will compute:
    // TT2 = (AA,BB,CC,DD)

    if (hadamard_bool_1) {
        hadamard(&TT2, &A->null_point);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT2, &A->null_point);
    }

    // we compute out->codomain.null_point = (AA,sqaabb, sqaacc, sqaadd)
    // out->precomputation = (  BB * CC *DD , sqaabb * CC * DD , sqaacc * BB* DD , sqaadd * BB * CC)
    fp2_copy(&out->codomain.null_point.x, &TT2.x);
    fp2_mul(&out->codomain.null_point.y, &TT2.x, &TT2.y);
    fp2_mul(&out->codomain.null_point.z, &TT2.x, &TT2.z);
    fp2_mul(&out->codomain.null_point.t, &TT2.x, &TT2.t);
    // No need to check the square roots, only used for signing.
    fp2_sqrt(&out->codomain.null_point.y);
    fp2_sqrt(&out->codomain.null_point.z);
    fp2_sqrt(&out->codomain.null_point.t);

    fp2_mul(&out->precomputation.x, &TT2.z, &TT2.t);
    fp2_mul(&out->precomputation.y,
            &out->precomputation.x,
            &out->codomain.null_point.y);                            // done for out->precomputation.y
    fp2_mul(&out->precomputation.x, &out->precomputation.x, &TT2.y); // done for out->precomputation.x
    fp2_mul(&out->precomputation.z, &TT2.t, &out->codomain.null_point.z);
    fp2_mul(&out->precomputation.z, &out->precomputation.z, &TT2.y); // done for out->precomputation.z
    fp2_mul(&out->precomputation.t, &TT2.z, &out->codomain.null_point.t);
    fp2_mul(&out->precomputation.t, &out->precomputation.t, &TT2.y); // done for out->precomputation.t

    if (hadamard_bool_2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}


static inline void
choose_index_theta_point(fp2_t *res, int ind, const theta_point_t *T)
{
    const fp2_t *src = NULL;
    switch (ind % 4) {
        case 0:
            src = &T->x;
            break;
        case 1:
            src = &T->y;
            break;
        case 2:
            src = &T->z;
            break;
        case 3:
            src = &T->t;
            break;
        default:
            //assert(0);
            break;
    }
    fp2_copy(res, src);
}


// Select a base change matrix in constant time, with M1 a regular
// base change matrix and M2 a precomputed base change matrix
// If option = 0 then M <- M1, else if option = 0xFF...FF then M <- M2
static inline void
select_base_change_matrix(basis_change_matrix_t *M,
                          const basis_change_matrix_t *M1,
                          const precomp_basis_change_matrix_t *M2,
                          const uint32_t option)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            fp2_select(&M->m[i][j], &M1->m[i][j], &FP2_CONSTANTS[M2->m[i][j]], option);
}


static bool
splitting_compute(theta_splitting_t *out, const theta_structure_t *A, int zero_index, bool randomize)

{
    // init
    uint32_t ctl;
    uint32_t count = 0;
    fp2_t U_cst, t1, t2;

    memset(&out->M, 0, sizeof(basis_change_matrix_t));

    // enumerate through all indices
    for (int i = 0; i < 10; i++) {
        fp2_set_zero(&U_cst);
        for (int t = 0; t < 4; t++) {
            // Iterate through the null point
            choose_index_theta_point(&t2, t, &A->null_point);
            choose_index_theta_point(&t1, t ^ EVEN_INDEX[i][1], &A->null_point);

            // Compute t1 * t2
            fp2_mul(&t1, &t1, &t2);
            // If CHI_EVAL(i,t) is +1 we want ctl to be 0 and
            // If CHI_EVAL(i,t) is -1 we want ctl to be 0xFF..FF
            ctl = (uint32_t)(CHI_EVAL[EVEN_INDEX[i][0]][t] >> 1);
            //assert(ctl == 0 || ctl == 0xffffffff);

            fp2_neg(&t2, &t1);
            fp2_select(&t1, &t1, &t2, ctl);

            // Then we compute U_cst ± (t1 * t2)
            fp2_add(&U_cst, &U_cst, &t1);
        }

        // If U_cst is 0 then update the splitting matrix
        ctl = fp2_is_zero(&U_cst);
        count -= ctl;
        select_base_change_matrix(&out->M, &out->M, &SPLITTING_TRANSFORMS[i], ctl);
        if (zero_index != -1 && i == zero_index &&
            !ctl) { // extra checks if we know exactly where the 0 index should be
            return 0;
        }
    }

#if defined(ENABLE_SIGN)
    // Pick a random normalization matrix
    if (randomize) {
        unsigned char secret_index = sample_random_index();
        basis_change_matrix_t Mrandom;

        set_base_change_matrix_from_precomp(&Mrandom, &NORMALIZATION_TRANSFORMS[0]);

        // Use a constant time selection to pick the index we want
        for (unsigned char i = 1; i < 6; i++) {
            // When i == secret_index, mask == 0 and 0xFF..FF otherwise
            int32_t mask = i - secret_index;
            mask = (mask | -mask) >> 31;
            select_base_change_matrix(&Mrandom, &Mrandom, &NORMALIZATION_TRANSFORMS[i], ~mask);
        }
        base_change_matrix_multiplication(&out->M, &Mrandom, &out->M);
    }
#else
    //assert(!randomize);
#endif

    // apply the isomorphism to ensure the null point is compatible with splitting
    apply_isomorphism(&out->B.null_point, &out->M, &A->null_point);

    // splitting was successful only if exactly one zero was identified
    return count == 1;
}

uint32_t
is_product_theta_point(const theta_point_t *P)
{
    fp2_t t1, t2;
    fp2_mul(&t1, &P->x, &P->t);
    fp2_mul(&t2, &P->y, &P->z);
    return fp2_is_equal(&t1, &t2);
}



static int
theta_product_structure_to_elliptic_product(theta_couple_curve_t *E12, theta_structure_t *A)
{
    fp2_t xx, yy;

    // This should be true from our computations in splitting_compute
    // but still check this for sanity
    if (!is_product_theta_point(&A->null_point))
        return 0;

    ec_curve_init(&(E12->E1));
    ec_curve_init(&(E12->E2));

    // A valid elliptic theta null point has no zero coordinate
    if (fp2_is_zero(&A->null_point.x) | fp2_is_zero(&A->null_point.y) | fp2_is_zero(&A->null_point.z))
        return 0;

    // xx = x², yy = y²
    fp2_sqr(&xx, &A->null_point.x);
    fp2_sqr(&yy, &A->null_point.y);
    // xx = x^4, yy = y^4
    fp2_sqr(&xx, &xx);
    fp2_sqr(&yy, &yy);

    // A2 = -2(x^4+y^4)/(x^4-y^4)
    fp2_add(&E12->E2.A, &xx, &yy);
    fp2_sub(&E12->E2.C, &xx, &yy);
    fp2_add(&E12->E2.A, &E12->E2.A, &E12->E2.A);
    fp2_neg(&E12->E2.A, &E12->E2.A);

    // same with x,z
    fp2_sqr(&xx, &A->null_point.x);
    fp2_sqr(&yy, &A->null_point.z);
    fp2_sqr(&xx, &xx);
    fp2_sqr(&yy, &yy);

    // A1 = -2(x^4+z^4)/(x^4-z^4)
    fp2_add(&E12->E1.A, &xx, &yy);
    fp2_sub(&E12->E1.C, &xx, &yy);
    fp2_add(&E12->E1.A, &E12->E1.A, &E12->E1.A);
    fp2_neg(&E12->E1.A, &E12->E1.A);

    if (fp2_is_zero(&E12->E1.C) | fp2_is_zero(&E12->E2.C))
        return 0;

    return 1;
}



static int
theta_point_to_montgomery_point(theta_couple_point_t *P12, const theta_point_t *P, const theta_structure_t *A)
{
    fp2_t temp;
    const fp2_t *x, *z;

    if (!is_product_theta_point(P))
        return 0;

    x = &P->x;
    z = &P->y;
    if (fp2_is_zero(x) & fp2_is_zero(z)) {
        x = &P->z;
        z = &P->t;
    }
    if (fp2_is_zero(x) & fp2_is_zero(z)) {
        return 0; // at this point P=(0:0:0:0) so is invalid
    }
    // P2.X = A.null_point.y * P.x + A.null_point.x * P.y
    // P2.Z = - A.null_point.y * P.x + A.null_point.x * P.y
    fp2_mul(&P12->P2.x, &A->null_point.y, x);
    fp2_mul(&temp, &A->null_point.x, z);
    fp2_sub(&P12->P2.z, &temp, &P12->P2.x);
    fp2_add(&P12->P2.x, &P12->P2.x, &temp);

    x = &P->x;
    z = &P->z;
    if (fp2_is_zero(x) & fp2_is_zero(z)) {
        x = &P->y;
        z = &P->t;
    }
    // P1.X = A.null_point.z * P.x + A.null_point.x * P.z
    // P1.Z = -A.null_point.z * P.x + A.null_point.x * P.z
    fp2_mul(&P12->P1.x, &A->null_point.z, x);
    fp2_mul(&temp, &A->null_point.x, z);
    fp2_sub(&P12->P1.z, &temp, &P12->P1.x);
    fp2_add(&P12->P1.x, &P12->P1.x, &temp);
    return 1;
}


static int
_theta_chain_compute_impl(unsigned n,
                          theta_couple_curve_t *E12,
                          const theta_kernel_couple_points_t *ker,
                          bool extra_torsion,
                          theta_couple_curve_t *E34,
                          theta_couple_point_t *P12,
                          size_t numP,
                          bool verify,
                          bool randomize)
{
    theta_structure_t theta;

    // lift the basis
    theta_couple_jac_point_t xyT1, xyT2;

    ec_basis_t bas1 = { .P = ker->T1.P1, .Q = ker->T2.P1, .PmQ = ker->T1m2.P1 };
    ec_basis_t bas2 = { .P = ker->T1.P2, .Q = ker->T2.P2, .PmQ = ker->T1m2.P2 };
    if (!lift_basis(&xyT1.P1, &xyT2.P1, &bas1, &E12->E1))
        return 0;
    if (!lift_basis(&xyT1.P2, &xyT2.P2, &bas2, &E12->E2))
        return 0;

    const unsigned extra = HD_extra_torsion * extra_torsion;

#ifndef NDEBUG
    //assert(extra == 0 || extra == 2); // only cases implemented
    if (!test_point_order_twof(&bas2.P, &E12->E2, n + extra))
        debug_print("bas2.P does not have correct order");

    if (!test_jac_order_twof(&xyT2.P2, &E12->E2, n + extra))
        debug_print("xyT2.P2 does not have correct order");
#endif

    theta_point_t pts[numP ? numP : 1];

    int space = 1;
    for (unsigned i = 1; i < n; i *= 2)
        ++space;

    uint16_t todo[space];
    todo[0] = n - 2 + extra;

    int current = 0;

    // kernel points for the gluing isogeny
    theta_couple_jac_point_t jacQ1[space], jacQ2[space];
    jacQ1[0] = xyT1;
    jacQ2[0] = xyT2;
    while (todo[current] != 1) {
        //assert(todo[current] >= 2);
        ++current;
        //assert(current < space);
        // the gluing isogeny is quite a bit more expensive than the others,
        // so we adjust the usual splitting rule here a little bit: towards
        // the end of the doubling chain it will be cheaper to recompute the
        // doublings after evaluation than to push the intermediate points.
        const unsigned num_dbls = todo[current - 1] >= 16 ? todo[current - 1] / 2 : todo[current - 1] - 1;
        //assert(num_dbls && num_dbls < todo[current - 1]);
        double_couple_jac_point_iter(&jacQ1[current], num_dbls, &jacQ1[current - 1], E12);
        double_couple_jac_point_iter(&jacQ2[current], num_dbls, &jacQ2[current - 1], E12);
        todo[current] = todo[current - 1] - num_dbls;
    }

    // kernel points for the remaining isogeny steps
    theta_point_t thetaQ1[space], thetaQ2[space];

    // the gluing step
    theta_gluing_t first_step;
    {
       // assert(todo[current] == 1);

        // compute the gluing isogeny
        if (!gluing_compute(&first_step, E12, &jacQ1[current], &jacQ2[current], verify))
            return 0;

        // evaluate
        for (unsigned j = 0; j < numP; ++j) {
            //assert(ec_is_zero(&P12[j].P1) || ec_is_zero(&P12[j].P2));
            if (!gluing_eval_point_special_case(&pts[j], &P12[j], &first_step))
                return 0;
        }

        // push kernel points through gluing isogeny
        for (int j = 0; j < current; ++j) {
            gluing_eval_basis(&thetaQ1[j], &thetaQ2[j], &jacQ1[j], &jacQ2[j], &first_step);
            --todo[j];
        }

        --current;
    }

    // set-up the theta_structure for the first codomain
    theta.null_point = first_step.codomain;
    theta.precomputation = 0;
    theta_precomputation(&theta);

    theta_isogeny_t step;

    // and now we do the remaining steps
    for (unsigned i = 1; current >= 0 && todo[current]; ++i) {
        //assert(current < space);
        while (todo[current] != 1) {
            //assert(todo[current] >= 2);
            ++current;
            //assert(current < space);
            const unsigned num_dbls = todo[current - 1] / 2;
            //assert(num_dbls && num_dbls < todo[current - 1]);
            double_iter(&thetaQ1[current], &theta, &thetaQ1[current - 1], num_dbls);
            double_iter(&thetaQ2[current], &theta, &thetaQ2[current - 1], num_dbls);
            todo[current] = todo[current - 1] - num_dbls;
        }

        // computing the next step
        int ret;
        if (i == n - 2) // penultimate step
            ret = theta_isogeny_compute(&step, &theta, &thetaQ1[current], &thetaQ2[current], 0, 0, verify);
        else if (i == n - 1) // ultimate step
            ret = theta_isogeny_compute(&step, &theta, &thetaQ1[current], &thetaQ2[current], 1, 0, false);
        else
            ret = theta_isogeny_compute(&step, &theta, &thetaQ1[current], &thetaQ2[current], 0, 1, verify);
        if (!ret)
            return 0;

        for (unsigned j = 0; j < numP; ++j)
            theta_isogeny_eval(&pts[j], &step, &pts[j]);

        // updating the codomain
        theta = step.codomain;

        // pushing the kernel
        //assert(todo[current] == 1);
        for (int j = 0; j < current; ++j) {
            theta_isogeny_eval(&thetaQ1[j], &step, &thetaQ1[j]);
            theta_isogeny_eval(&thetaQ2[j], &step, &thetaQ2[j]);
            //assert(todo[j]);
            --todo[j];
        }

        --current;
    }

    //assert(current == -1);

    if (!extra_torsion) {
        if (n >= 3) {
            // in the last step we've skipped pushing the kernel since current was == 0, let's do it now
            theta_isogeny_eval(&thetaQ1[0], &step, &thetaQ1[0]);
            theta_isogeny_eval(&thetaQ2[0], &step, &thetaQ2[0]);
        }

        // penultimate step
        theta_isogeny_compute_4(&step, &theta, &thetaQ1[0], &thetaQ2[0], 0, 0);
        for (unsigned j = 0; j < numP; ++j)
            theta_isogeny_eval(&pts[j], &step, &pts[j]);
        theta = step.codomain;
        theta_isogeny_eval(&thetaQ1[0], &step, &thetaQ1[0]);
        theta_isogeny_eval(&thetaQ2[0], &step, &thetaQ2[0]);

        // ultimate step
        theta_isogeny_compute_2(&step, &theta, &thetaQ1[0], &thetaQ2[0], 1, 0);
        for (unsigned j = 0; j < numP; ++j)
            theta_isogeny_eval(&pts[j], &step, &pts[j]);
        theta = step.codomain;
    }

    // final splitting step
    theta_splitting_t last_step;

    bool is_split = splitting_compute(&last_step, &theta, extra_torsion ? 8 : -1, randomize);

    if (!is_split) {
        debug_print("kernel did not generate an isogeny between elliptic products");
        return 0;
    }

    if (!theta_product_structure_to_elliptic_product(E34, &last_step.B))
        return 0;

    // evaluate
    for (size_t j = 0; j < numP; ++j) {
        apply_isomorphism(&pts[j], &last_step.M, &pts[j]);
        if (!theta_point_to_montgomery_point(&P12[j], &pts[j], &last_step.B))
            return 0;
    }

    return 1;
}

// Like theta_chain_compute_and_eval, adding extra verification checks;
// used in the signature verification
int
theta_chain_compute_and_eval_verify(unsigned n,
                                    /*const*/ theta_couple_curve_t *E12,
                                    const theta_kernel_couple_points_t *ker,
                                    bool extra_torsion,
                                    theta_couple_curve_t *E34,
                                    theta_couple_point_t *P12,
                                    size_t numP)
{
    return _theta_chain_compute_impl(n, E12, ker, extra_torsion, E34, P12, numP, true, false);
}


// The commitment curve can be recovered from the codomain of the 2D
// isogeny built from the bases computed during verification.
static int
compute_commitment_curve_verify(ec_curve_t *E_com,
                                const ec_basis_t *B_chall_can,
                                const ec_basis_t *B_aux_can,
                                const ec_curve_t *E_chall,
                                const ec_curve_t *E_aux,
                                int pow_dim2_deg_resp) {
#ifndef NDEBUG
    // Check all the points are the correct order
    if (!test_basis_order_twof(B_chall_can, E_chall, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("B_chall_can does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");

    if (!test_basis_order_twof(B_aux_can, E_aux, HD_extra_torsion + pow_dim2_deg_resp))
        debug_print("B_aux_can does not have order 2^(HD_extra_torsion + pow_dim2_deg_resp");
#endif

    // now compute the dim2 isogeny from Echall x E_aux -> E_com x E_aux'
    // of kernel B_chall_can x B_aux_can

    // first we set-up the kernel
    theta_couple_curve_t EchallxEaux;
    copy_curve(&EchallxEaux.E1, E_chall);
    copy_curve(&EchallxEaux.E2, E_aux);

    theta_kernel_couple_points_t dim_two_ker;
    copy_bases_to_kernel(&dim_two_ker, B_chall_can, B_aux_can);

    // computing the isogeny
    theta_couple_curve_t codomain;
    int codomain_splits;
    ec_curve_init(&codomain.E1);
    ec_curve_init(&codomain.E2);
    // handling the special case where we don't need to perform any dim2 computation
    if (pow_dim2_deg_resp == 0) {
        codomain_splits = 1;
        copy_curve(&codomain.E1, &EchallxEaux.E1);
        copy_curve(&codomain.E2, &EchallxEaux.E2);
        // We still need to check that E_chall is supersingular
        // This assumes that HD_extra_torsion == 2
        if (!ec_is_basis_four_torsion(B_chall_can, E_chall)) {
            return 0;
        }
    } else {
        codomain_splits = theta_chain_compute_and_eval_verify(
            pow_dim2_deg_resp, &EchallxEaux, &dim_two_ker, true, &codomain, NULL, 0);
    }

    // computing the commitment curve
    // its always the first one because of our (2^n,2^n)-isogeny formulae
    copy_curve(E_com, &codomain.E1);

    return codomain_splits;
}


int
ec_curve_verify_A(const fp2_t *A)
{ // Verify the Montgomery coefficient A is valid (A^2-4 \ne 0)
    // Return 1 if curve is valid, 0 otherwise
    fp2_t t;
    fp2_set_one(&t);
    fp_add(&t.re, &t.re, &t.re); // t=2
    if (fp2_is_equal(A, &t))
        return 0;
    fp_neg(&t.re, &t.re); // t=-2
    if (fp2_is_equal(A, &t))
        return 0;
    return 1;
}


int
ec_curve_init_from_A(ec_curve_t *E, const fp2_t *A)
{ // Initialize the curve from the A coefficient and check it is valid
    // Return 1 if curve is valid, 0 otherwise
    ec_curve_init(E);
    fp2_copy(&E->A, A); // Set A
    return ec_curve_verify_A(A);
}

// SQIsign verification
int
protocols_verify(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l) {
    int verify;

    if (!check_canonical_basis_change_matrix(sig))
        return 0;

    // Computation of the length of the dim 2 2^n isogeny
    int pow_dim2_deg_resp = SQIsign_response_length - (int) sig->two_resp_length - (int) sig->backtracking;

    // basic sanity test: checking that the response is not too long
    if (pow_dim2_deg_resp < 0)
        return 0;
    // The dim 2 isogeny embeds a dim 1 isogeny of odd degree, so it can
    // never be of length 2.
    if (pow_dim2_deg_resp == 1)
        return 0;

    // check the public curve is valid
    if (!ec_curve_verify_A(&(pk->curve).A))
        return 0;

    // Set auxiliary curve from the A-coefficient within the signature
    ec_curve_t E_aux;
    if (!ec_curve_init_from_A(&E_aux, &sig->E_aux_A))
        return 0; // invalid curve

    // checking that we are given A-coefficients and no precomputation
   // assert(fp2_is_one(&pk->curve.C) == 0xFFFFFFFF && !pk->curve.is_A24_computed_and_normalized);

    // computation of the challenge
    ec_curve_t E_chall;
    if (!compute_challenge_verify(&E_chall, sig, &pk->curve, pk->hint_pk)) {
        return 0;
    }

    // Computation of the canonical bases for the challenge and aux curve
    ec_basis_t B_chall_can, B_aux_can;

    if (!challenge_and_aux_basis_verify(&B_chall_can, &B_aux_can, &E_chall, &E_aux, sig, pow_dim2_deg_resp)) {
        return 0;
    }

    // When two_resp_length != 0 we need to compute a second, short 2^r-isogeny
    if (sig->two_resp_length > 0) {
        if (!two_response_isogeny_verify(&E_chall, &B_chall_can, sig, pow_dim2_deg_resp)) {
            return 0;
        }
    }

    // We can recover the commitment curve with a 2D isogeny
    // The supplied signature did not compute an isogeny between eliptic products
    // and so definitely is an invalid signature.
    ec_curve_t E_com;
    if (!compute_commitment_curve_verify(&E_com, &B_chall_can, &B_aux_can, &E_chall, &E_aux, pow_dim2_deg_resp))
        return 0;

    scalar_t chk_chall;

    // recomputing the challenge vector
    hash_to_challenge(&chk_chall, pk, &E_com, m, l);

    // performing the final check
    verify = mp_compare(sig->chall_coeff, chk_chall, NWORDS_ORDER) == 0;

    return verify;
}
