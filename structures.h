//
// Created by gustavo on 18/04/25.
//

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "fp_constants.h"
#include "fp.h"

typedef struct fp2_t
{
    fp_t re, im;
} fp2_t;

typedef digit_t scalar_t[NWORDS_ORDER];
typedef scalar_t scalar_mtx_2x2_t[2][2];



/** @defgroup ec_t Data structures
 * @{
 */

/** @brief Projective point on the Kummer line E/pm 1 in Montgomery coordinates
 *
 * @typedef ec_point_t
 *
 * @struct ec_point_t
 *
 * A projective point in (X:Z) or (X:Y:Z) coordinates (tbd).
 */
typedef struct ec_point_t
{
    fp2_t x;
    fp2_t z;
} ec_point_t;

/** @brief An elliptic curve
 *
 * @typedef ec_curve_t
 *
 * @struct ec_curve_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct ec_curve_t
{
    fp2_t A;
    fp2_t C;                             ///< cannot be 0
    ec_point_t A24;                      // the point (A+2 : 4C)
    bool is_A24_computed_and_normalized; // says if A24 has been computed and normalized
} ec_curve_t;


/** @brief Type for the signature
 *
 * @typedef signature_t
 *
 * @struct signature
 *
 */
typedef struct signature
{
    fp2_t E_aux_A; // the Montgomery A-coefficient for the auxiliary curve
    uint8_t backtracking;
    uint8_t two_resp_length;
    scalar_mtx_2x2_t mat_Bchall_can_to_B_chall; // the matrix of the desired basis
    scalar_t chall_coeff;
    uint8_t hint_aux;
    uint8_t hint_chall;
} signature_t;

/** @brief Type for the public keys
 *
 * @typedef public_key_t
 *
 * @struct public_key
 *
 */
typedef struct public_key
{
    ec_curve_t curve; // the normalized A-coefficient of the Montgomery curve
    uint8_t hint_pk;
} public_key_t;


/** @brief A basis of a torsion subgroup
 *
 * @typedef ec_basis_t
 *
 * @struct ec_basis_t
 *
 * A pair of points (or a triplet, tbd) forming a basis of a torsion subgroup.
 */
typedef struct ec_basis_t
{
    ec_point_t P;
    ec_point_t Q;
    ec_point_t PmQ;
} ec_basis_t;

/** @brief An isogeny of degree a power of 2
 *
 * @typedef ec_isog_even_t
 *
 * @struct ec_isog_even_t
 */
typedef struct ec_isog_even_t
{
    ec_curve_t curve;  ///< The domain curve
    ec_point_t kernel; ///< A kernel generator
    unsigned length;   ///< The length as a 2-isogeny walk
} ec_isog_even_t;


/* KPS structure for isogenies of degree 2 or 4 */
typedef struct
{
    ec_point_t K;
} ec_kps2_t;
typedef struct
{
    ec_point_t K[3];
} ec_kps4_t;


/** @brief Type for couple curve *
 * @typedef theta_couple_curve_t
 *
 * @struct theta_couple_curve
 *
 * the  theta_couple_curve structure
 */
typedef struct theta_couple_curve
{
    ec_curve_t E1;
    ec_curve_t E2;
} theta_couple_curve_t;

/** @brief Type for a product E1 x E2 with corresponding bases
 * @typedef theta_couple_curve_with_basis_t
 *
 * @struct theta_couple_curve_with_basis
 *
 * tType for a product E1 x E2 with corresponding bases Ei[2^n]
 */
typedef struct theta_couple_curve_with_basis
{
    ec_curve_t E1;
    ec_curve_t E2;
    ec_basis_t B1;
    ec_basis_t B2;
} theta_couple_curve_with_basis_t;

/** @brief Type for couple point with XZ coordinates
 * @typedef theta_couple_point_t
 *
 * @struct theta_couple_point
 *
 * Structure for the couple point on an elliptic product
 * using XZ coordinates
 */
typedef struct theta_couple_point
{
    ec_point_t P1;
    ec_point_t P2;
} theta_couple_point_t;

/** @brief Type for three couple points T1, T2, T1-T2 with XZ coordinates
 * @typedef theta_kernel_couple_points_t
 *
 * @struct theta_kernel_couple_points
 *
 * Structure for a triple of theta couple points T1, T2 and T1 - T2
 */
typedef struct theta_kernel_couple_points
{
    theta_couple_point_t T1;
    theta_couple_point_t T2;
    theta_couple_point_t T1m2;
} theta_kernel_couple_points_t;

/** @brief Type for theta point with repeating components
 * @typedef theta_point_compact_t
 *
 * @struct theta_point_compact
 *
 * the  theta_point structure used for points with repeated components
 */
typedef struct theta_point_compact
{
    fp2_t x;
    fp2_t y;
} theta_point_compact_t;


/** @brief Type for theta point *
 * @typedef theta_point_t
 *
 * @struct theta_point
 *
 * the  theta_point structure used
 */
typedef struct theta_point
{
    fp2_t x;
    fp2_t y;
    fp2_t z;
    fp2_t t;
} theta_point_t;

/** @brief Type for theta structure *
 * @typedef theta_structure_t
 *
 * @struct theta_structure
 *
 * the  theta_structure structure used
 */
typedef struct theta_structure
{
    theta_point_t null_point;
    bool precomputation;

    // Eight precomputed values used for doubling and
    // (2,2)-isogenies.
    fp2_t XYZ0;
    fp2_t YZT0;
    fp2_t XZT0;
    fp2_t XYT0;

    fp2_t xyz0;
    fp2_t yzt0;
    fp2_t xzt0;
    fp2_t xyt0;
} theta_structure_t;

/** @brief Projective point in Montgomery coordinates
 *
 * @typedef jac_point_t
 *
 * @struct jac_point_t
 *
 * A projective point in (X:Y:Z) coordinates
 */
typedef struct jac_point_t
{
    fp2_t x;
    fp2_t y;
    fp2_t z;
} jac_point_t;


/** @brief Type for couple point with XYZ coordinates
 * @typedef theta_couple_jac_point_t
 *
 * @struct theta_couple_jac_point
 *
 * Structure for the couple point on an elliptic product
 * using XYZ coordinates
 */
typedef struct theta_couple_jac_point
{
    jac_point_t P1;
    jac_point_t P2;
} theta_couple_jac_point_t;





/** @brief A 4x4 matrix used for basis changes
 * @typedef basis_change_matrix_t
 *
 * @struct basis_change_matrix
 *
 * Structure to hold 16 elements representing a 4x4 matrix used for changing
 * the basis of a theta point.
 */
typedef struct basis_change_matrix
{
    fp2_t m[4][4];
} basis_change_matrix_t;

/** @brief Type for gluing (2,2) theta isogeny *
 * @typedef theta_gluing_t
 *
 * @struct theta_gluing
 *
 * the  theta_gluing structure
 */
typedef struct theta_gluing
{

    theta_couple_curve_t domain;
    theta_couple_jac_point_t xyK1_8;
    theta_point_compact_t imageK1_8;
    basis_change_matrix_t M;
    theta_point_t precomputation;
    theta_point_t codomain;

} theta_gluing_t;

/** @brief A 2x2 matrix used for action by translation
 * @typedef translation_matrix_t
 *
 * @struct translation_matrix
 *
 * Structure to hold 4 fp2_t elements representing a 2x2 matrix used when computing
 * a compatible theta structure during gluing.
 */
typedef struct translation_matrix
{
    fp2_t g00;
    fp2_t g01;
    fp2_t g10;
    fp2_t g11;
} translation_matrix_t;

/** @brief Addition components
 *
 * @typedef add_components_t
 *
 * @struct add_components_t
 *
 * 3 components u,v,w that define the (X:Z) coordinates of both
 * addition and substraction of two distinct points with
 * P+Q =(u-v:w) and P-Q = (u+v=w)
 */
typedef struct add_components_t
{
    fp2_t u;
    fp2_t v;
    fp2_t w;
} add_components_t;

/** @brief Type for standard (2,2) theta isogeny *
 * @typedef theta_isogeny_t
 *
 * @struct theta_isogeny
 *
 * the  theta_isogeny structure
 */
typedef struct theta_isogeny
{
    theta_point_t T1_8;
    theta_point_t T2_8;
    bool hadamard_bool_1;
    bool hadamard_bool_2;
    theta_structure_t domain;
    theta_point_t precomputation;
    theta_structure_t codomain;
} theta_isogeny_t;

/** @brief Type for splitting isomorphism *
 * @typedef theta_splitting_t
 *
 * @struct theta_splitting
 *
 * the theta_splitting structure
 */
typedef struct theta_splitting
{
    basis_change_matrix_t M;
    theta_structure_t B;

} theta_splitting_t;

typedef struct precomp_basis_change_matrix {
    uint8_t m[4][4];
} precomp_basis_change_matrix_t;


#define FP2_ZERO 0
#define FP2_ONE 1
#define FP2_I 2
#define FP2_MINUS_ONE 3
#define FP2_MINUS_I 4

extern const fp2_t FP2_CONSTANTS[5];
extern const precomp_basis_change_matrix_t SPLITTING_TRANSFORMS[10];
extern const precomp_basis_change_matrix_t NORMALIZATION_TRANSFORMS[6] ;
extern const int CHI_EVAL[4][4];
extern const int EVEN_INDEX[10][2];



#endif //STRUCTURES_H
