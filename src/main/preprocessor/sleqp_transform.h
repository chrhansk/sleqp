#ifndef SLEQP_TRANSFORM_H
#define SLEQP_TRANSFORM_H

#include "sleqp_problem.h"

#include "preprocessor/sleqp_preprocessing_state.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpTransformation SleqpTransformation;

  SLEQP_RETCODE sleqp_transformation_create(SleqpTransformation** star,
                                            SleqpPreprocessingState* preprocessing_state,
                                            SleqpParams* params);

  SLEQP_RETCODE sleqp_transformation_convert_primal(SleqpTransformation* transformation,
                                                    const SleqpSparseVec* source,
                                                    SleqpSparseVec* target);

  SLEQP_RETCODE sleqp_transformation_create_transformed_problem(SleqpTransformation* transformation,
                                                                SleqpProblem** star);

  SLEQP_RETCODE sleqp_transformation_capture(SleqpTransformation* transformation);

  SLEQP_RETCODE sleqp_transformation_release(SleqpTransformation** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TRANSFORM_H */
