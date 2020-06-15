#ifndef SLEQP_HESS_STRUCT_H
#define SLEQP_HESS_STRUCT_H

#include "sleqp_types.h"

/**
 * @file sleqp_solver.h
 * @brief Definition of problem block structure.
 **/

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpHessianStruct SleqpHessianStruct;

  SLEQP_RETCODE sleqp_hessian_struct_create(SleqpHessianStruct** star,
                                            int dimension,
                                            bool empty);

  int sleqp_hessian_struct_get_num_blocks(SleqpHessianStruct* hessian_struct);

  SLEQP_RETCODE sleqp_hessian_struct_get_block_range(SleqpHessianStruct* hessian_struct,
                                                     int block,
                                                     int* begin,
                                                     int* end);

  SLEQP_RETCODE sleqp_hessian_struct_push_block(SleqpHessianStruct* hessian_struct,
                                                int end);

  SLEQP_RETCODE sleqp_hessian_struct_clear(SleqpHessianStruct* hessian_struct);

  SLEQP_RETCODE sleqp_hessian_struct_get_linear_range(SleqpHessianStruct* hessian_struct,
                                                      int* begin,
                                                      int* end);

  SLEQP_RETCODE sleqp_hessian_struct_fprintf(SleqpHessianStruct* hessian_struct,
                                             FILE* output);

  SLEQP_RETCODE sleqp_hessian_struct_free(SleqpHessianStruct** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_HESS_STRUCT_H */
