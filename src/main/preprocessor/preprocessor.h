#ifndef SLEQP_PREPROCESSOR_H
#define SLEQP_PREPROCESSOR_H

#include "iterate.h"
#include "problem.h"
#include "types.h"

/**
 * @file preprocessor.h
 * @brief Definition of the preprocessor.
 **/

typedef struct SleqpPreprocessor SleqpPreprocessor;

SLEQP_RETCODE
sleqp_preprocessor_create(SleqpPreprocessor** star,
                          SleqpProblem* problem,
                          SleqpParams* params);

SLEQP_PREPROCESSING_RESULT
sleqp_preprocessor_result(SleqpPreprocessor* preprocessor);

SleqpProblem*
sleqp_preprocessor_transformed_problem(SleqpPreprocessor* preprocessor);

SLEQP_RETCODE
sleqp_preprocessor_transform_primal(SleqpPreprocessor* preprocessor,
                                    const SleqpVec* source,
                                    SleqpVec* target);

SleqpTimer*
sleqp_preprocessor_get_timer(SleqpPreprocessor* preprocessor);

SLEQP_RETCODE
sleqp_preprocessor_restore_iterate(SleqpPreprocessor* preprocessor,
                                   const SleqpIterate* transformed_iterate,
                                   SleqpIterate* original_iterate);

SLEQP_RETCODE
sleqp_preprocessor_capture(SleqpPreprocessor* preprocessor);

SLEQP_RETCODE
sleqp_preprocessor_release(SleqpPreprocessor** star);

#endif /* SLEQP_PREPROCESSOR_H */
