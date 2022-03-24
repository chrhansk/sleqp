#include "ampl_problem.h"

#include <assert.h>

#include "ampl_func.h"
#include "ampl_util.h"

static SLEQP_RETCODE
compute_linear_coeffs(SleqpSparseMatrix* linear_coeffs, SleqpAmplData* data)
{
  ASL* asl = data->asl;

  jacval(data->x, data->jac_vals, NULL);

  int next_col = 0;

  for (int i = 0; i < nzc; ++i)
  {
    int row    = data->jac_rows[i];
    int col    = data->jac_cols[i];
    double val = data->jac_vals[i];

    while (col >= next_col)
    {
      SLEQP_CALL(sleqp_sparse_matrix_push_column(linear_coeffs, next_col++));
    }

    if (row < nlc)
    {
      continue;
    }

    row -= nlc;

    SLEQP_CALL(sleqp_sparse_matrix_push(linear_coeffs, row, col, val));
  }

  while (n_var > next_col)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push_column(linear_coeffs, next_col++));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
apply_linear_offset(const SleqpSparseMatrix* linear_coeffs,
                    SleqpParams* params,
                    SleqpAmplData* data)
{
  ASL* asl = data->asl;

  conval(data->x, data->cons_val, NULL);

  const int num_linear = n_con - nlc;

  double* linear_val = data->cons_val + nlc;

  const int linear_nnz = sleqp_sparse_matrix_nnz(linear_coeffs);

  const double* linear_data = sleqp_sparse_matrix_data(linear_coeffs);
  const int* linear_cols    = sleqp_sparse_matrix_cols(linear_coeffs);
  const int* linear_rows    = sleqp_sparse_matrix_rows(linear_coeffs);

  for (int col = 0; col < n_var; ++col)
  {
    for (int entry = linear_cols[col]; entry < linear_cols[col + 1]; ++entry)
    {
      assert(entry < linear_nnz);

      const int row      = linear_rows[entry];
      const double value = linear_data[entry];

      linear_val[row] -= value * data->x[col];
    }
  }

  double* linear_lb = data->cons_lb + nlc;
  double* linear_ub = data->cons_ub + nlc;

  const double inf = sleqp_infinity();

  for (int i = 0; i < num_linear; ++i)
  {
    if (linear_lb[i] != -inf)
    {
      linear_lb[i] -= linear_val[i];
    }

    if (linear_ub[i] != -inf)
    {
      linear_ub[i] -= linear_val[i];
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_problem_create(SleqpProblem** star,
                          SleqpAmplData* data,
                          SleqpParams* params,
                          bool halt_on_error)
{
  const int num_variables   = data->num_variables;
  const int num_constraints = data->num_constraints;

  ASL* asl = data->asl;

  SLEQP_CALL(map_ampl_inf(data->var_lb, num_variables));
  SLEQP_CALL(map_ampl_inf(data->var_ub, num_variables));

  SLEQP_CALL(map_ampl_inf(data->cons_lb, num_constraints));
  SLEQP_CALL(map_ampl_inf(data->cons_ub, num_constraints));

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SleqpVec* var_lb;
  SleqpVec* var_ub;

  SleqpVec* cons_lb;
  SleqpVec* cons_ub;

  SleqpFunc* func;

  SLEQP_CALL(sleqp_vec_create_empty(&var_lb, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&var_ub, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&cons_lb, num_constraints));
  SLEQP_CALL(sleqp_vec_create_empty(&cons_ub, num_constraints));

  SLEQP_CALL(sleqp_vec_from_raw(var_lb, data->var_lb, num_variables, zero_eps));
  SLEQP_CALL(sleqp_vec_from_raw(var_ub, data->var_ub, num_variables, zero_eps));

  const int num_linear = n_con - nlc;

  SLEQP_CALL(sleqp_ampl_func_create(&func, data, params, halt_on_error));

  if (num_linear > 0)
  {
    SleqpSparseMatrix* linear_coeffs;
    SleqpVec* linear_lb;
    SleqpVec* linear_ub;

    SleqpVec* general_lb = cons_lb;
    SleqpVec* general_ub = cons_ub;

    SLEQP_CALL(
      sleqp_sparse_matrix_create(&linear_coeffs, num_linear, n_var, nzc));

    SLEQP_CALL(sleqp_vec_create_empty(&linear_lb, num_linear));
    SLEQP_CALL(sleqp_vec_create_empty(&linear_ub, num_linear));

    SLEQP_CALL(compute_linear_coeffs(linear_coeffs, data));

    SLEQP_CALL(apply_linear_offset(linear_coeffs, params, data));

    general_lb->dim = nlc;
    general_ub->dim = nlc;

    SLEQP_CALL(sleqp_vec_from_raw(general_lb, data->cons_lb, nlc, zero_eps));
    SLEQP_CALL(sleqp_vec_from_raw(general_ub, data->cons_ub, nlc, zero_eps));

    SLEQP_CALL(
      sleqp_vec_from_raw(linear_lb, data->cons_lb + nlc, num_linear, zero_eps));
    SLEQP_CALL(
      sleqp_vec_from_raw(linear_ub, data->cons_ub + nlc, num_linear, zero_eps));

    SLEQP_CALL(sleqp_problem_create(star,
                                    func,
                                    params,
                                    var_lb,
                                    var_ub,
                                    general_lb,
                                    general_ub,
                                    linear_coeffs,
                                    linear_lb,
                                    linear_ub));

    SLEQP_CALL(sleqp_vec_free(&linear_ub));
    SLEQP_CALL(sleqp_vec_free(&linear_lb));

    SLEQP_CALL(sleqp_sparse_matrix_release(&linear_coeffs));
  }
  else
  {
    SLEQP_CALL(
      sleqp_vec_from_raw(cons_lb, data->cons_lb, num_constraints, zero_eps));
    SLEQP_CALL(
      sleqp_vec_from_raw(cons_ub, data->cons_ub, num_constraints, zero_eps));

    SLEQP_CALL(sleqp_problem_create_simple(star,
                                           func,
                                           params,
                                           var_lb,
                                           var_ub,
                                           cons_lb,
                                           cons_ub));
  }

  SLEQP_CALL(sleqp_func_release(&func));

  SLEQP_CALL(sleqp_vec_free(&cons_ub));
  SLEQP_CALL(sleqp_vec_free(&cons_lb));

  SLEQP_CALL(sleqp_vec_free(&var_ub));
  SLEQP_CALL(sleqp_vec_free(&var_lb));

  return SLEQP_OKAY;
}
