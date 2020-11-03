#ifndef HSL_MA57_H
#define HSL_MA57_H

/**
 * @file  hsl_ma57.h
 * @brief A sensible interface to the Fortran HSL_MA57 solver
 * @note  Fortran 4-byte integers expected
 **/

#include <stdint.h>

/** @brief MA57 error, warning, and success codes
 */
typedef enum {
  MA57_ERROR_INSUFFICIENT_INTEGER_SPACE = +11,
  MA57_ERROR_INSUFFICIENT_REAL_SPACE    = +10,
  MA57_SUCCESS_ZERO_SOLUTION              = +8,
  MA57_SUCCESS_INDEFINITE                 = +5,
  MA57_SUCCESS_RANK_DEFICIENT             = +4,
  MA57_SUCCESS_ENTRIES_IGNORED            = +3,
  MA57_SUCCESS_DUPLICATE_ENTRIES_IGNORED  = +2,
  MA57_SUCCESS_INVALID_ENTRIES_IGNORED    = +1,

  MA57_SUCCESS                    = 0,

  MA57_N_OUT_OF_RANGE             = -1,
  MA57_NE_OUT_OF_RANGE            = -2,
  MA57_INSUFFICIENT_REAL_SPACE    = -3,
  MA57_INSUFFICIENT_INTEGER_SPACE = -4,
  MA57_TINY_PIVOT                 = -5,
  MA57_PIVOT_SIGN_CHANGE          = -6,
  MA57_COPY_TO_SMALLER_ARRAY      = -7,
  MA57_NO_REFINEMENT_CONVERGENCE  = -8,
  MA57_PERMUTATION_ERROR          = -9,
  MA57_PIVOTING_CONTROL_ERROR     = -10,
  MA57_LRHS_OUT_OF_RANGE          = -11,
  MA57_JOB_OUT_OF_RANGE           = -12,
  MA57_REFINEMENT_CONTROL_ERROR   = -13,
  MA57_CONDITION_ESTIMATE_ERROR   = -14,
  MA57_LKEEP_OUT_OF_RANGE         = -15,
  MA57_NRHS_OUT_OF_RANGE          = -16,
  MA57_LWORK_OUT_OF_RANGE         = -17,
  MA57_METIS_PACKAGE_MISSING      = -18
} MA57_ERROR;


void ma57id_(double *cntl,
             int32_t *icntl);

void ma57ad_(const int32_t* n,
             const int32_t* ne,
             const int32_t* irn,
             const int32_t* jcn,
             const int32_t* lkeep,
             int32_t* keep,
             int32_t* iwork,
             int32_t* icntl,
             int32_t* info,
             double* rinfo);

void ma57bd_(const int32_t* n,
             const int32_t* ne,
             const double* a,
             double* fact,
             const int32_t* lfact,
             int32_t* ifact,
             const int32_t* lifact,
             const int32_t* lkeep,
             int32_t* keep,
             int32_t* iwork,
             int32_t* icntl,
             double* cntl,
             int32_t* info,
             double* rinfo);

void ma57cd_(const int32_t* job,
             const int32_t* n,
             const double* fact,
             const int32_t* lfact,
             const int32_t* ifact,
             const int32_t* lifact,
             const int32_t* nrhs,
             double* rhs,
             const int32_t* lrhs,
             double* work,
             const int32_t* lwork,
             int32_t* iwork,
             int32_t* icntl,
             int32_t* info);

void ma57ed_(const int32_t* n,
             const int32_t* ic,
             int32_t* keep,
             const double* fact,
             const int32_t* lfact,
             double* newfact,
             const int32_t* lnew,
             const int32_t* ifact,
             const int32_t* lifact,
             int32_t* newifc,
             const int32_t* linew,
             int32_t* info);

#endif // HSL_MA57_H
