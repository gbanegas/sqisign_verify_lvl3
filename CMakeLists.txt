cmake_minimum_required(VERSION 4.0)
project(verification_skii C)

set(CMAKE_C_STANDARD 99)

add_executable(verification_skii main.c
        sqisign.c
        sqisign.h
        structures.h
        fp.h
        fp_constants.h
        tutil.h
        encoded_sizes.h
        encoded_verification.c
        verification.h
        mp.c
        mp.h
        ec_params.h
        e0_basis.c
        e0_basis.h
        fp2.c
        fp.c
        fips202.c
        fp_p65376_32.c
        constants.c
        hash_cha.c
)
