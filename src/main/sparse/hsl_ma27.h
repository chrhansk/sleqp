#ifndef HSL_MA27_H
#define HSL_MA27_H

/**
 * @file  hsl_ma27.h
 * @brief A sensible interface to the Fortran HSL_MA27 solver
 * @note  Fortran 4-byte integers expected
 **/

#include <stdint.h>

typedef enum
{
        MA27_NSTEPS_OUT_OF_RANGE = -7,
        MA27_PIVOT_SIGN_CHANGE   = -6,
        MA27_SINGULAR_MATRIX     = -5,
        MA27_A_MEM_TOO_SMALL     = -4,
        MA27_IW_MEM_TOO_SMALL    = -3,
        MA27_NZ_OUT_OF_RANGE     = -2,
        MA27_N_OUT_OF_RANGE      = -1,

        MA27_SUCCESS = 0,

        MA27_WARN_IRN_ICN_OUT_OF_RANGE = +1,
        MA27_WARN_INDEFINITE           = +2,
        MA27_WARN_RANK_DEFICIENT       = +3
} MA27_ERROR;

//
// public functions
//

void ma27id_(int32_t *icntl,
             double *cntl);

void ma27ad_(const int32_t *n,
             const int32_t *nz,
             const int32_t *irn,
             const int32_t *icn,
             int32_t       *iw,
             const int32_t *liw,
             int32_t       *ikeep,
             int32_t       *iw1,
             int32_t       *nsteps,
             const int32_t *iflag,
             int32_t       *icntl,
             double        *cntl,
             int32_t       *info,
             double        *ops);

void ma27bd_(const int32_t *n,
             const int32_t *nz,
             const int32_t *irn,
             const int32_t *icn,
             double        *a,
             const int32_t *la,
             int32_t       *iw,
             const int32_t *liw,
             const int32_t *ikeep,
             const int32_t *nsteps,
             int32_t       *maxfrt,
             int32_t       *iw1,
             int32_t       *icntl,
             double        *cntl,
             int32_t       *info);

void ma27cd_(const int32_t *n,
             const double  *a,
             const int32_t *la,
             const int32_t *iw,
             const int32_t *liw,
             double        *w,
             const int32_t *maxfrt,
             double        *x_rhs,
             int32_t       *iw1,
             const int32_t *nsteps,
             int32_t       *icntl,
             int32_t       *info);


#endif // HSL_MA27_H
