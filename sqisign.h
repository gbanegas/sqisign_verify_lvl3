//
// Created by gustavo on 18/04/25.
//

#ifndef SQISIGN_H
#define SQISIGN_H

#include "encoded_sizes.h"
#include "verification.h"
#include "fp_constants.h"
#include "fp2.h"
#include "fp.h"
#include "tutil.h"
#include "structures.h"
#include "mp.h"


#ifndef NDEBUG
#define DEBUG_PRINT 1
#else
#define DEBUG_PRINT 0
#endif

#define HD_extra_torsion 2

#define debug_print(fmt)                                                                           \
    do {                                                                                           \
        if (DEBUG_PRINT)                                                                           \
            printf("warning: %s, file %s, line %d, function %s().\n",                              \
                   fmt,                                                                            \
                   __FILE_NAME__,                                                                  \
                   __LINE__,                                                                       \
                   __func__);                                                                      \
    } while (0)

int
sqisign_verify(const unsigned char *m,
               unsigned long long mlen,
               const unsigned char *sig,
               unsigned long long siglen,
               const unsigned char *pk);

#endif //SQISIGN_H
