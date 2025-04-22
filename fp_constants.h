#ifndef FP_CONSTANTS_H
#define FP_CONSTANTS_H
#include <stddef.h>
#include <stdint.h>

#define RADIX 32

#define NWORDS_FIELD 14
#define NWORDS_ORDER 12
#define BITS 384
#define LOG2P 9
#define digit_t uint32_t
#define sdigit_t int32_t
#define LOG2RADIX 5

typedef digit_t fp_t[NWORDS_FIELD];



// Datatype for representing field elements

#endif

