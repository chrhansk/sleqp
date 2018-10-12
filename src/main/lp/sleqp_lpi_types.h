#ifndef SLEQP_LPI_TYPES_H
#define SLEQP_LPI_TYPES_H

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpLPi SleqpLPi;

  enum SLEQP_BaseStat
  {
    SLEQP_BASESTAT_LOWER = 0,             /**< variable is at its lower bound */
    SLEQP_BASESTAT_BASIC = 1,             /**< variable is basic */
    SLEQP_BASESTAT_UPPER = 2,             /**< variable is at its upper bound */
    SLEQP_BASESTAT_ZERO  = 3              /**< free variable is non-basic and set to zero */
  };

  typedef enum SLEQP_BaseStat SLEQP_BASESTAT;

  typedef SLEQP_RETCODE (*SLEQP_LPI_CREATE)(void** lp_data,
                                            int num_variables,
                                            int num_constraints);

  typedef SLEQP_RETCODE (*SLEQP_LPI_SOLVE)(void* lp_data,
                                           double* objective,
                                           SleqpSparseMatrix* matrix,
                                           double* cons_lb,
                                           double* cons_ub,
                                           double* vars_lb,
                                           double* vars_ub);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_SOLUTION)(void* lp_data,
                                                  int num_variables,
                                                  double* objective_value,
                                                  double* solution_values);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_VARSTATS)(void* lp_data,
                                                  int num_variables,
                                                  SLEQP_BASESTAT* variable_stats);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_CONSSTATS)(void* lp_data,
                                                   int num_constraints,
                                                   SLEQP_BASESTAT* constraintstats);

  typedef SLEQP_RETCODE (*SLEQP_LPI_FREE)(void** lp_data);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LPI_TYPES_H */
