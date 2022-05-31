#ifndef SLEQP_LPI_TYPES_H
#define SLEQP_LPI_TYPES_H

/**
 * @file lpi_types.h
 * @brief Definition of the LP interface types.
 **/

#include "types.h"

#include "options.h"
#include "params.h"

#include "sparse/sparse_matrix.h"
#include "sparse/vec.h"

typedef enum
{
  SLEQP_BASESTAT_LOWER = 0, /**< variable is at its lower bound */
  SLEQP_BASESTAT_BASIC = 1, /**< variable is basic */
  SLEQP_BASESTAT_UPPER = 2, /**< variable is at its upper bound */
  SLEQP_BASESTAT_ZERO  = 3  /**< free variable is non-basic and set to zero */
} SLEQP_BASESTAT;

typedef enum
{
  SLEQP_LP_STATUS_UNKNOWN,
  SLEQP_LP_STATUS_OPTIMAL,
  SLEQP_LP_STATUS_INF,
  SLEQP_LP_STATUS_INF_OR_UNBOUNDED,
  SLEQP_LP_STATUS_UNBOUNDED,
} SLEQP_LP_STATUS;

typedef SLEQP_RETCODE (*SLEQP_LPI_CREATE)(void** lp_data,
                                          int num_variables,
                                          int num_constraints,
                                          SleqpParams* params,
                                          SleqpOptions* options);

typedef SLEQP_RETCODE (*SLEQP_LPI_SOLVE)(void* lp_data,
                                         int num_variables,
                                         int num_constraints,
                                         double time_limit);

typedef SLEQP_LP_STATUS (*SLEQP_LPI_STATUS)(void* lp_data);

typedef SLEQP_RETCODE (*SLEQP_LPI_SET_BOUNDS)(void* lp_data,
                                              int num_variables,
                                              int num_constraints,
                                              double* cons_lb,
                                              double* cons_ub,
                                              double* vars_lb,
                                              double* vars_ub);

typedef SLEQP_RETCODE (*SLEQP_LPI_SET_COEFFS)(void* lp_data,
                                              int num_variables,
                                              int num_constraints,
                                              SleqpSparseMatrix* cons_matrix);

typedef SLEQP_RETCODE (*SLEQP_LPI_SET_OBJ)(void* lp_data,
                                           int num_variables,
                                           int num_constraints,
                                           double* objective);

typedef SLEQP_RETCODE (*SLEQP_LPI_SET_BASIS)(void* lp_data,
                                             int index,
                                             const SLEQP_BASESTAT* var_stats,
                                             const SLEQP_BASESTAT* cons_stats);

typedef SLEQP_RETCODE (*SLEQP_LPI_SAVE_BASIS)(void* lp_data, int index);

typedef SLEQP_RETCODE (*SLEQP_LPI_RESTORE_BASIS)(void* lp_data, int index);

typedef SLEQP_RETCODE (*SLEQP_LPI_PRIMAL_SOL)(void* lp_data,
                                              int num_variables,
                                              int num_constraints,
                                              double* objective_value,
                                              double* solution_values);

typedef SLEQP_RETCODE (*SLEQP_LPI_DUAL_SOL)(void* lp_data,
                                            int num_variables,
                                            int num_constraints,
                                            double* vars_dual,
                                            double* cons_dual);

typedef SLEQP_RETCODE (*SLEQP_LPI_VARS_STATS)(void* lp_data,
                                              int num_variables,
                                              int num_constraints,
                                              SLEQP_BASESTAT* var_stats);

typedef SLEQP_RETCODE (*SLEQP_LPI_CONS_STATS)(void* lp_data,
                                              int num_variables,
                                              int num_constraints,
                                              SLEQP_BASESTAT* cons_stats);

typedef SLEQP_RETCODE (*SLEQP_LPI_BASIS_CONDITION_ESTIMATE)(void* lp_data,
                                                            bool* exact,
                                                            double* condition);

typedef SLEQP_RETCODE (*SLEQP_LPI_WRITE)(void* lp_data, const char* filename);

typedef SLEQP_RETCODE (*SLEQP_LPI_FREE)(void** lp_data);

typedef struct
{
  SLEQP_LPI_CREATE create_problem;
  SLEQP_LPI_SOLVE solve;
  SLEQP_LPI_STATUS status;
  SLEQP_LPI_SET_BOUNDS set_bounds;
  SLEQP_LPI_SET_COEFFS set_coeffs;
  SLEQP_LPI_SET_OBJ set_obj;
  SLEQP_LPI_SET_BASIS set_basis;
  SLEQP_LPI_SAVE_BASIS save_basis;
  SLEQP_LPI_RESTORE_BASIS restore_basis;
  SLEQP_LPI_PRIMAL_SOL primal_sol;
  SLEQP_LPI_DUAL_SOL dual_sol;
  SLEQP_LPI_VARS_STATS vars_stats;
  SLEQP_LPI_CONS_STATS cons_stats;
  SLEQP_LPI_BASIS_CONDITION_ESTIMATE basis_condition_estimate;
  SLEQP_LPI_WRITE write;
  SLEQP_LPI_FREE free_problem;
} SleqpLPiCallbacks;

#endif /* SLEQP_LPI_TYPES_H */
