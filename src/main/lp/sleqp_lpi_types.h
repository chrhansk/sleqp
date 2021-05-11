#ifndef SLEQP_LPI_TYPES_H
#define SLEQP_LPI_TYPES_H

/**
 * @file sleqp_lpi_types.h
 * @brief Definition of the LP interface types.
 **/

#include "sleqp_types.h"

#include "sleqp_params.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpLPi SleqpLPi;

  typedef enum
  {
    SLEQP_BASESTAT_LOWER = 0,             /**< variable is at its lower bound */
    SLEQP_BASESTAT_BASIC = 1,             /**< variable is basic */
    SLEQP_BASESTAT_UPPER = 2,             /**< variable is at its upper bound */
    SLEQP_BASESTAT_ZERO  = 3              /**< free variable is non-basic and set to zero */
  } SLEQP_BASESTAT;

  typedef SLEQP_RETCODE (*SLEQP_LPI_CREATE)(void** lp_data,
                                            int num_variables,
                                            int num_constraints,
                                            SleqpParams* params);

  typedef SLEQP_RETCODE (*SLEQP_LPI_SOLVE)(void* lp_data,
                                           int num_variables,
                                           int num_constraints,
                                           double time_limit);

  typedef SLEQP_RETCODE (*SLEQP_LPI_SET_BOUNDS)(void* lp_data,
                                                int num_variables,
                                                int num_constraints,
                                                double* cons_lb,
                                                double* cons_ub,
                                                double* vars_lb,
                                                double* vars_ub);

  typedef SLEQP_RETCODE (*SLEQP_LPI_SET_COEFFICIENTS)(void* lp_data,
                                                      int num_variables,
                                                      int num_constraints,
                                                      SleqpSparseMatrix* cons_matrix);

  typedef SLEQP_RETCODE (*SLEQP_LPI_SET_OBJECTIVE)(void* lp_data,
                                                   int num_variables,
                                                   int num_constraints,
                                                   double* objective);

  typedef SLEQP_RETCODE (*SLEQP_LPI_SAVE_BASIS)(void* lp_data,
                                                int index);

  typedef SLEQP_RETCODE (*SLEQP_LPI_RESTORE_BASIS)(void* lp_data,
                                                   int index);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_PRIMAL_SOL)(void* lp_data,
                                                    int num_variables,
                                                    int num_constraints,
                                                    double* objective_value,
                                                    double* solution_values);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_DUAL_SOL)(void* lp_data,
                                                  int num_variables,
                                                  int num_constraints,
                                                  double* vars_dual,
                                                  double* cons_dual);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_VARSTATS)(void* lp_data,
                                                  int num_variables,
                                                  int num_constraints,
                                                  SLEQP_BASESTAT* variable_stats);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_CONSSTATS)(void* lp_data,
                                                   int num_variables,
                                                   int num_constraints,
                                                   SLEQP_BASESTAT* constraintstats);

  typedef SLEQP_RETCODE (*SLEQP_LPI_GET_BASIS_CONDITION)(void *lp_data,
                                                         bool* exact,
                                                         double* condition);

  typedef SLEQP_RETCODE (*SLEQP_LPI_FREE)(void** lp_data);

  typedef struct {
    SLEQP_LPI_CREATE create_problem;
    SLEQP_LPI_SOLVE solve;
    SLEQP_LPI_SET_BOUNDS set_bounds;
    SLEQP_LPI_SET_COEFFICIENTS set_coefficients;
    SLEQP_LPI_SET_OBJECTIVE set_objective;
    SLEQP_LPI_SAVE_BASIS save_basis;
    SLEQP_LPI_RESTORE_BASIS restore_basis;
    SLEQP_LPI_GET_PRIMAL_SOL get_primal_sol;
    SLEQP_LPI_GET_DUAL_SOL get_dual_sol;
    SLEQP_LPI_GET_VARSTATS get_varstats;
    SLEQP_LPI_GET_CONSSTATS get_consstats;
    SLEQP_LPI_GET_BASIS_CONDITION get_basis_condition;
    SLEQP_LPI_FREE free_problem;
  } SleqpLPiCallbacks;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LPI_TYPES_H */
