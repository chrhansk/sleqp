#ifndef SLEQP_PREPROCESSOR_H
#define SLEQP_PREPROCESSOR_H

#include "sleqp_iterate.h"
#include "sleqp_types.h"
#include "sleqp_problem.h"

/**
 * @file sleqp_pre.h
 * @brief Definition of the preprocessor.
 **/

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpPreprocessor SleqpPreprocessor;

  SLEQP_RETCODE sleqp_preprocessor_create(SleqpPreprocessor** star,
                                          SleqpProblem* problem,
                                          SleqpParams* params);

  SLEQP_PREPROCESSING_RESULT sleqp_preprocessor_result(SleqpPreprocessor* preprocessor);

  SleqpProblem* sleqp_preprocessor_transformed_problem(SleqpPreprocessor* preprocessor);

  SLEQP_RETCODE sleqp_preprocessor_transform_primal(const SleqpSparseVec* source,
                                                    SleqpSparseVec* target);

  SLEQP_RETCODE sleqp_preprocessor_restore_iterate(SleqpPreprocessor* preprocessor,
                                                   const SleqpIterate* transformed_iterate,
                                                   SleqpIterate* original_iterate);

  SLEQP_RETCODE sleqp_preprocessor_capture(SleqpPreprocessor* preprocessor);

  SLEQP_RETCODE sleqp_preprocessor_release(SleqpPreprocessor** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PREPROCESSOR_H */
