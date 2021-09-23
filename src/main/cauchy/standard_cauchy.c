#include "standard_cauchy.h"

#include <math.h>

#include "fail.h"
#include "cmp.h"
#include "mem.h"
#include "problem.h"
#include "merit.h"

typedef struct
{
  SleqpProblem* problem;
  SleqpParams* params;
  SleqpOptions* options;

  int num_lp_variables;
  int num_lp_constraints;

  double trust_radius;
  SleqpIterate* iterate;

  SLEQP_BASESTAT* var_stats;
  SLEQP_BASESTAT* cons_stats;

  bool has_basis[SLEQP_NUM_CAUCHY_OBJECTIVES];
  SLEQP_CAUCHY_OBJECTIVE_TYPE current_objective;

  bool has_coefficients;

  SleqpLPi* lp_interface;

  double* objective;
  double* cons_lb;
  double* cons_ub;
  double* vars_lb;
  double* vars_ub;

  double* solution_values;
  double* dual_values;

  SleqpSparseVec* quadratic_gradient;
} CauchyData;

static SLEQP_RETCODE
cauchy_data_create(CauchyData** star,
                   SleqpProblem* problem,
                   SleqpParams* params,
                   SleqpOptions* options,
                   SleqpLPi* lp_interface)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables = sleqp_problem_num_variables(problem);

  const int num_constraints = sleqp_problem_num_constraints(problem);

  CauchyData* data = *star;

  *data = (CauchyData) {0};

  data->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(data->problem));

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  data->options = options;

  data->num_lp_variables = num_variables + 2 * num_constraints;
  data->num_lp_constraints = num_constraints;

  data->trust_radius = SLEQP_NONE;

  assert(data->num_lp_variables == sleqp_lpi_get_num_variables(lp_interface));
  assert(data->num_lp_constraints == sleqp_lpi_get_num_constraints(lp_interface));

  SLEQP_CALL(sleqp_alloc_array(&data->var_stats, data->num_lp_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->cons_stats, data->num_lp_constraints));

  data->has_basis[SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT] = false;
  data->has_basis[SLEQP_CAUCHY_OBJECTIVE_TYPE_FEASIBILITY] = false;
  data->has_basis[SLEQP_CAUCHY_OBJECTIVE_TYPE_MIXED] = false;
  data->current_objective = SLEQP_NONE;

  data->has_coefficients = false;

  SLEQP_CALL(sleqp_lpi_capture(lp_interface));
  data->lp_interface = lp_interface;

  SLEQP_CALL(sleqp_alloc_array(&data->objective, data->num_lp_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->cons_lb, data->num_lp_constraints));
  SLEQP_CALL(sleqp_alloc_array(&data->cons_ub, data->num_lp_constraints));

  SLEQP_CALL(sleqp_alloc_array(&data->vars_lb, data->num_lp_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->vars_ub, data->num_lp_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->solution_values, data->num_lp_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->dual_values,
                               SLEQP_MAX(data->num_lp_constraints,
                                         data->num_lp_variables)));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->quadratic_gradient,
                                              num_variables));

  const double inf = sleqp_infinity();

  for(int j = num_variables; j < data->num_lp_variables; ++j)
  {
    data->vars_lb[j] = 0;
    data->vars_ub[j] = inf;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
append_identities(SleqpSparseMatrix* cons_jac,
                  int num_variables,
                  int num_constraints)
{
  assert(num_constraints == sleqp_sparse_matrix_get_num_rows(cons_jac));
  assert(num_variables == sleqp_sparse_matrix_get_num_cols(cons_jac));

  const int nnz = sleqp_sparse_matrix_get_nnz(cons_jac);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(cons_jac);
  const int num_cols = sleqp_sparse_matrix_get_num_cols(cons_jac);

  /*
   * Reserve a litte more so we can add the two
   * identity matrices afterwards
   */
  SLEQP_CALL(sleqp_sparse_matrix_reserve(cons_jac,
                                         nnz + 2*num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_resize(cons_jac,
                                        num_rows,
                                        num_cols + 2*num_constraints));

  // append a + id
  for(int i = 0; i < num_constraints; ++i)
  {
    const int con = num_variables + i;

    SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac,
                                               con));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        i,
                                        con,
                                        1.));
  }

  // append a - id
  for(int i = 0; i < num_constraints; ++i)
  {
    const int con = num_variables + num_constraints + i;

    SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac,
                                               con));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        i,
                                        con,
                                        -1.));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
remove_identities(SleqpSparseMatrix* cons_jac,
                  int num_variables,
                  int num_constraints)
{
  const int num_rows = sleqp_sparse_matrix_get_num_rows(cons_jac);
  const int num_cols = sleqp_sparse_matrix_get_num_cols(cons_jac);

  assert(num_constraints == num_rows);
  assert(num_variables + 2*num_constraints == num_cols);

  SLEQP_CALL(sleqp_sparse_matrix_resize(cons_jac,
                                        num_constraints,
                                        num_variables));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_cons_bounds(CauchyData* cauchy_data,
                   SleqpIterate* iterate,
                   int num_variables,
                   int num_constraints)
{
  int k_c = 0, k_lb = 0, k_ub = 0;

  SleqpSparseVec* lb = sleqp_problem_cons_lb(cauchy_data->problem);
  SleqpSparseVec* ub = sleqp_problem_cons_ub(cauchy_data->problem);

  SleqpSparseVec* val = sleqp_iterate_get_cons_val(iterate);

  const double inf = sleqp_infinity();

  for(int i = 0; i < num_constraints; ++i)
  {
    while(k_c < val->nnz && val->indices[k_c] < i)
    {
      ++k_c;
    }

    while(k_lb < lb->nnz && lb->indices[k_lb] < i)
    {
      ++k_lb;
    }

    while(k_ub < ub->nnz && ub->indices[k_ub] < i)
    {
      ++k_ub;
    }

    const bool valid_ub = (k_ub < ub->nnz && ub->indices[k_ub] == i);
    const bool valid_lb = (k_lb < lb->nnz && lb->indices[k_lb] == i);
    const bool valid_c = (k_c < val->nnz && val->indices[k_c] == i);

    const double ubval = valid_ub ? ub->data[k_ub] : 0;
    const double lbval = valid_lb ? lb->data[k_lb] : 0;
    const double cval = valid_c ? val->data[k_c] : 0;

    assert(!sleqp_is_infinite(lbval));
    assert(!sleqp_is_infinite(-ubval));
    assert(sleqp_is_finite(cval));

    if(sleqp_is_infinite(ubval))
    {
      cauchy_data->cons_ub[i] = inf;
    }
    else
    {
      cauchy_data->cons_ub[i] = ubval - cval;
    }

    if(sleqp_is_infinite(-lbval))
    {
      cauchy_data->cons_lb[i] = -inf;
    }
    else
    {
      cauchy_data->cons_lb[i] = lbval - cval;
    }

    assert(cauchy_data->cons_lb[i] <= cauchy_data->cons_ub[i]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_var_bounds(CauchyData* cauchy_data,
                  SleqpIterate* iterate,
                  int num_variables,
                  int num_constraints)
{
  SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
  SleqpProblem* problem = cauchy_data->problem;

  SleqpSparseVec* lb = sleqp_problem_var_lb(problem);
  SleqpSparseVec* ub = sleqp_problem_var_ub(problem);

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  int k_x = 0, k_lb = 0, k_ub = 0;

  for(int j = 0; j < num_variables; ++j)
  {
    while(k_x < x->nnz && x->indices[k_x] < j)
    {
      ++k_x;
    }

    while(k_lb < lb->nnz && lb->indices[k_lb] < j)
    {
      ++k_lb;
    }

    while(k_ub < ub->nnz && ub->indices[k_ub] < j)
    {
      ++k_ub;
    }

    const bool valid_x = (k_x < x->nnz) && (j == x->indices[k_x]);
    const bool valid_ub = (k_ub < ub->nnz) && (j == ub->indices[k_ub]);
    const bool valid_lb = (k_lb < lb->nnz) && (j == lb->indices[k_lb]);

    const double ubval = valid_ub ? ub->data[k_ub] : 0.;
    const double lbval = valid_lb ? lb->data[k_lb] : 0.;
    const double xval = valid_x ? x->data[k_x] : 0.;

    assert(!sleqp_is_infinite(lbval));
    assert(!sleqp_is_infinite(-ubval));
    assert(sleqp_is_finite(xval));

    if(sleqp_is_infinite(ubval))
    {
      cauchy_data->vars_ub[j] = trust_radius;
    }
    else
    {
      cauchy_data->vars_ub[j] = SLEQP_MIN(ubval - xval, trust_radius);
    }

    if(sleqp_is_infinite(-lbval))
    {
      cauchy_data->vars_lb[j] = -trust_radius;
    }
    else
    {
      cauchy_data->vars_lb[j] = SLEQP_MAX(lbval - xval, -trust_radius);
    }

    assert(cauchy_data->vars_lb[j] <= cauchy_data->vars_ub[j]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_objective(CauchyData* cauchy_data,
                 SleqpSparseVec* gradient,
                 double penalty)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  const int num_lp_variables = sleqp_lpi_get_num_variables(cauchy_data->lp_interface);

  assert(num_lp_variables == num_variables + 2*num_constraints);

  for(int j = 0; j < num_variables; ++j)
  {
    cauchy_data->objective[j] = 0.;
  }

  for(int j = num_variables; j < num_lp_variables; ++j)
  {
    cauchy_data->objective[j] = penalty;
  }

  if(gradient)
  {
    assert(gradient->dim == num_variables);

    for(int k = 0; k < gradient->nnz; ++k)
    {
      cauchy_data->objective[gradient->indices[k]] = gradient->data[k];
    }
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE cauchy_set_coefficients(CauchyData* cauchy_data,
                                      SleqpSparseMatrix* cons_jac)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_sparse_matrix_get_num_cols(cons_jac);
  const int num_constraints = sleqp_sparse_matrix_get_num_rows(cons_jac);

  bool fixed_jacobian = !(sleqp_problem_has_nonlinear_cons(problem));

  if(fixed_jacobian && cauchy_data->has_coefficients)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(append_identities(cons_jac,
                               num_variables,
                               num_constraints));

  assert(sleqp_sparse_matrix_is_valid(cons_jac));

  SLEQP_CALL(sleqp_lpi_set_coefficients(cauchy_data->lp_interface,
                                        cons_jac));

  SLEQP_CALL(remove_identities(cons_jac,
                               num_variables,
                               num_constraints));

  assert(sleqp_sparse_matrix_is_valid(cons_jac));

  cauchy_data->has_coefficients = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_set_iterate(SleqpIterate* iterate,
                            double trust_radius,
                            void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;
  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);

  const int num_variables = sleqp_sparse_matrix_get_num_cols(cons_jac);
  const int num_constraints = sleqp_sparse_matrix_get_num_rows(cons_jac);

  assert(trust_radius > 0.);

  cauchy_data->trust_radius = trust_radius;

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  cauchy_data->iterate = iterate;

  SLEQP_CALL(sleqp_iterate_capture(cauchy_data->iterate));

  SLEQP_CALL(create_var_bounds(cauchy_data,
                               iterate,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(create_cons_bounds(cauchy_data,
                                iterate,
                                num_variables,
                                num_constraints));

  SLEQP_CALL(sleqp_lpi_set_bounds(cauchy_data->lp_interface,
                                  cauchy_data->cons_lb,
                                  cauchy_data->cons_ub,
                                  cauchy_data->vars_lb,
                                  cauchy_data->vars_ub));

  SLEQP_CALL(cauchy_set_coefficients(cauchy_data, cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_set_trust_radius(double trust_radius,
                                 void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(cauchy_data->iterate);

  const int num_variables = sleqp_sparse_matrix_get_num_cols(cons_jac);
  const int num_constraints = sleqp_sparse_matrix_get_num_rows(cons_jac);

  cauchy_data->trust_radius = trust_radius;

  SLEQP_CALL(create_var_bounds(cauchy_data,
                               cauchy_data->iterate,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(sleqp_lpi_set_bounds(cauchy_data->lp_interface,
                                  cauchy_data->cons_lb,
                                  cauchy_data->cons_ub,
                                  cauchy_data->vars_lb,
                                  cauchy_data->vars_ub));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_basis(CauchyData* cauchy_data,
            bool* valid_basis)
{
  (*valid_basis) = true;

  int basis_size = 0;

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->var_stats));

  for(int j = 0; j < cauchy_data->num_lp_variables; ++j)
  {
    SLEQP_BASESTAT var_stat = cauchy_data->var_stats[j];

    if(var_stat != SLEQP_BASESTAT_BASIC)
    {
      ++basis_size;
    }
  }

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->cons_stats));

  for(int i = 0; i < cauchy_data->num_lp_constraints; ++i)
  {
    SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

    if(cons_stat != SLEQP_BASESTAT_BASIC)
    {
      ++basis_size;
    }
  }

  (*valid_basis) = (basis_size == cauchy_data->num_lp_variables);

  return SLEQP_OKAY;
}

/*
  static SLEQP_RETCODE
  check_direction_bounds(SleqpCauchyData* cauchy_data,
  bool* valid_direction)
  {
  SleqpProblem* problem = cauchy_data->problem;

  SleqpIterate* iterate = cauchy_data->iterate;

  assert(iterate);

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  const double eps = sleqp_params_get_eps(cauchy_data->params);

  const int num_variables = problem->num_variables;

  SLEQP_CALL(sleqp_lpi_get_primal_sol(cauchy_data->lp_interface,
  NULL,
  cauchy_data->solution_values));

  (*valid_direction) = true;

  for(int j = 0; j < num_variables; ++j)
  {
  const double dval = cauchy_data->solution_values[j];

  if(sleqp_is_gt(dval, trust_radius, eps) || sleqp_is_lt(dval, -trust_radius, eps))
  {
  (*valid_direction) = false;
  }
  }

  {
  SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
  SleqpSparseVec* lb = problem->var_lb;
  SleqpSparseVec* ub = problem->var_ub;

  int k_x = 0, k_lb = 0, k_ub = 0;

  for(int j = 0; j < num_variables; ++j)
  {
  while(k_x < x->nnz && x->indices[k_x] < j)
  {
  ++k_x;
  }

  while(k_lb < lb->nnz && lb->indices[k_lb] < j)
  {
  ++k_lb;
  }

  while(k_ub < ub->nnz && ub->indices[k_ub] < j)
  {
  ++k_ub;
  }

  const bool valid_x = (k_x < x->nnz) && (j == x->indices[k_x]);
  const bool valid_ub = (k_ub < ub->nnz) && (j == ub->indices[k_ub]);
  const bool valid_lb = (k_lb < lb->nnz) && (j == lb->indices[k_lb]);

  const double ubval = valid_ub ? ub->data[k_ub] : 0.;
  const double lbval = valid_lb ? lb->data[k_lb] : 0.;
  const double xval = valid_x ? x->data[k_x] : 0.;
  const double dval = cauchy_data->solution_values[j];

  if(sleqp_is_gt(xval + dval, ubval, eps) || sleqp_is_lt(xval + dval, lbval, eps))
  {
  (*valid_direction) = false;
  }
  }
  }

  return SLEQP_OKAY;
  }
*/

static
SLEQP_RETCODE restore_basis(CauchyData* cauchy_data,
                            SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type)
{
  if(cauchy_data->current_objective != objective_type &&
     cauchy_data->has_basis[objective_type])
  {
    SLEQP_CALL(sleqp_lpi_restore_basis(cauchy_data->lp_interface,
                                       objective_type));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_solve(SleqpSparseVec* gradient,
                      double penalty,
                      SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type,
                      void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SLEQP_CALL(create_objective(cauchy_data,
                              gradient,
                              penalty));

  SLEQP_CALL(sleqp_lpi_set_objective(cauchy_data->lp_interface,
                                     cauchy_data->objective));

  bool warm_start = sleqp_options_get_bool(cauchy_data->options,
                                           SLEQP_OPTION_BOOL_ALWAYS_WARM_START_LP);

  if(warm_start)
  {
    switch(objective_type)
    {
    case SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT:
      // fallthrough
    case SLEQP_CAUCHY_OBJECTIVE_TYPE_FEASIBILITY:
      SLEQP_CALL(restore_basis(cauchy_data, objective_type));
      break;
    case SLEQP_CAUCHY_OBJECTIVE_TYPE_MIXED:
      if(cauchy_data->current_objective != SLEQP_CAUCHY_OBJECTIVE_TYPE_MIXED)
      {
        // restart from the default, this should be closer
        // to the initial mixed one
        SLEQP_CALL(restore_basis(cauchy_data, SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));
      }
      break;
    default:
      break;
    }
  }

  cauchy_data->current_objective = objective_type;

  SLEQP_CALL(sleqp_lpi_solve(cauchy_data->lp_interface));

  SLEQP_LPI_STATUS status = sleqp_get_status(cauchy_data->lp_interface);

  if(status != SLEQP_LPI_STATUS_OPTIMAL)
  {
    sleqp_log_error("Invalid LP status: %d", status);
    return SLEQP_INTERNAL_ERROR;
  }

#if !defined(NDEBUG)
  {
    bool valid_basis = false;

    SLEQP_CALL(check_basis(cauchy_data, &valid_basis));

    assert(valid_basis);
  }

  // TODO: Find out about tolerance guarantees given by the LP solvers
  /*
    {
    bool valid_direction = false;

    SLEQP_CALL(check_direction_bounds(cauchy_data, &valid_direction));

    assert(valid_direction);
    }
  */
#endif

  if(warm_start)
  {
    SLEQP_CALL(sleqp_lpi_save_basis(cauchy_data->lp_interface, objective_type));
  }

  cauchy_data->has_basis[objective_type] = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_get_objective_value(double* objective_value,
                                    void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SLEQP_CALL(sleqp_lpi_get_primal_sol(cauchy_data->lp_interface,
                                      objective_value,
                                      NULL));

  (*objective_value) += sleqp_iterate_get_func_val(cauchy_data->iterate);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_get_working_set(SleqpIterate* iterate,
                                void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  assert(sleqp_working_set_valid(working_set));

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->var_stats));

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->cons_stats));

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  const double eps = sleqp_params_get(cauchy_data->params,
                                      SLEQP_PARAM_EPS);

  {
    SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
    SleqpSparseVec* lb = sleqp_problem_var_lb(problem);
    SleqpSparseVec* ub = sleqp_problem_var_ub(problem);

    int k_x = 0, k_lb = 0, k_ub = 0;

    for(int i = 0; i < num_variables; ++i)
    {
      while(k_x < x->nnz && x->indices[k_x] < i)
      {
        ++k_x;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_x = (k_x < x->nnz) && (i == x->indices[k_x]);
      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;
      const double xval = valid_x ? x->data[k_x] : 0.;

      assert(cauchy_data->var_stats[i] != SLEQP_BASESTAT_ZERO);

      sleqp_assert_is_leq(lbval, xval, eps);
      sleqp_assert_is_leq(xval, ubval, eps);

      if(sleqp_is_eq(lbval, ubval, eps))
      {
        SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                  i,
                                                  SLEQP_ACTIVE_BOTH));
      }
      else if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER) &&
              ((xval - lbval) < trust_radius))
      {
        SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                  i,
                                                  SLEQP_ACTIVE_LOWER));
      }
      else if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER) &&
              ((ubval - xval) < trust_radius))
      {
        SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                  i,
                                                  SLEQP_ACTIVE_UPPER));
      }
    }
  }

  {
    int k_lb = 0, k_ub = 0;

    SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
    SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_constraints;

    SleqpSparseVec* lb = sleqp_problem_cons_lb(problem);
    SleqpSparseVec* ub = sleqp_problem_cons_ub(problem);

    for(int i = 0; i < num_constraints; ++i)
    {
      const SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

      if(cons_stat == SLEQP_BASESTAT_BASIC)
      {
        continue;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;

      assert(lower_slack_stats[i] != SLEQP_BASESTAT_BASIC ||
             upper_slack_stats[i] != SLEQP_BASESTAT_BASIC);

      const bool zero_slack = lower_slack_stats[i] == SLEQP_BASESTAT_LOWER &&
        upper_slack_stats[i] == SLEQP_BASESTAT_LOWER;

      if(cons_stat == SLEQP_BASESTAT_ZERO)
      {
        assert(sleqp_is_infinite(ubval));
        assert(sleqp_is_infinite(-lbval));

        continue;
      }


      if(sleqp_is_eq(lbval, ubval, eps))
      {
        if(zero_slack)
        {
          SLEQP_CALL(sleqp_working_set_add_constraint(working_set,
                                                      i,
                                                      SLEQP_ACTIVE_BOTH));
        }
      }
      else
      {
        if(cauchy_data->cons_stats[i] == SLEQP_BASESTAT_UPPER)
        {
          if(zero_slack)
          {

            // the constraint c(x) + grad(c(x))*d +I s_l -I s_u  <= u
            // is tight at i

            SLEQP_CALL(sleqp_working_set_add_constraint(working_set,
                                                        i,
                                                        SLEQP_ACTIVE_UPPER));
          }
        }
        else if(cauchy_data->cons_stats[i] == SLEQP_BASESTAT_LOWER)
        {
          if(zero_slack)
          {

            // the constraint l <= c(x) + grad(c(x))*d +I s_l -I s_u
            // is tight at i

            SLEQP_CALL(sleqp_working_set_add_constraint(working_set,
                                                        i,
                                                        SLEQP_ACTIVE_LOWER));
          }
        }
      }
    }
  }

  const int num_active_vars = sleqp_working_set_num_active_vars(working_set);
  const int num_active_cons = sleqp_working_set_num_active_cons(working_set);

  sleqp_log_debug("Created an active set with %d variables, %d constraints",
                  num_active_vars,
                  num_active_cons);

  assert(sleqp_working_set_valid(working_set));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_get_direction(SleqpSparseVec* direction,
                              void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  const double zero_eps = sleqp_params_get(cauchy_data->params, SLEQP_PARAM_ZERO_EPS);

  const int num_variables = sleqp_problem_num_variables(cauchy_data->problem);

  SLEQP_CALL(sleqp_lpi_get_primal_sol(cauchy_data->lp_interface,
                                      NULL,
                                      cauchy_data->solution_values));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(direction,
                                          cauchy_data->solution_values,
                                          num_variables,
                                          zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_locally_infeasible(bool* locally_infeasible,
                                   void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  *locally_infeasible = false;

  SleqpIterate* iterate = cauchy_data->iterate;

  assert(iterate);

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->var_stats));

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->cons_stats));

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  const double eps = sleqp_params_get(cauchy_data->params,
                                      SLEQP_PARAM_EPS);

  bool active_trust_region = false;

  // Check if trust region is active
  {
    SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
    SleqpSparseVec* lb = sleqp_problem_var_lb(problem);
    SleqpSparseVec* ub = sleqp_problem_var_ub(problem);

    int k_x = 0, k_lb = 0, k_ub = 0;

    for(int i = 0; i < num_variables; ++i)
    {
      while(k_x < x->nnz && x->indices[k_x] < i)
      {
        ++k_x;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_x = (k_x < x->nnz) && (i == x->indices[k_x]);
      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;
      const double xval = valid_x ? x->data[k_x] : 0.;

      assert(cauchy_data->var_stats[i] != SLEQP_BASESTAT_ZERO);

      sleqp_assert_is_leq(lbval, xval, eps);
      sleqp_assert_is_leq(xval, ubval, eps);

      if(sleqp_is_eq(lbval, ubval, eps))
      {
        continue;
      }
      else if(cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER &&
              (xval - lbval) >= trust_radius)
      {
        active_trust_region = true;
        break;
      }
      else if(cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER &&
              ((ubval - xval) >= trust_radius))
      {
        active_trust_region = true;
        break;
      }
    }
  }

  bool feasible_direction = true;

  {
    int k_lb = 0, k_ub = 0;

    SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
    SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_constraints;

    SleqpSparseVec* lb = sleqp_problem_cons_lb(problem);
    SleqpSparseVec* ub = sleqp_problem_cons_ub(problem);

    for(int i = 0; i < num_constraints; ++i)
    {
      const SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

      if(cons_stat == SLEQP_BASESTAT_BASIC)
      {
        continue;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      assert(lower_slack_stats[i] != SLEQP_BASESTAT_BASIC ||
             upper_slack_stats[i] != SLEQP_BASESTAT_BASIC);

      bool zero_slack = lower_slack_stats[i] == SLEQP_BASESTAT_LOWER &&
        upper_slack_stats[i] == SLEQP_BASESTAT_LOWER;

      if(!zero_slack)
      {
        feasible_direction = false;
      }
    }
  }

  sleqp_log_debug("Trust region active: %s, feasible direction: %s",
                  sleqp_bool_string(active_trust_region),
                  sleqp_bool_string(feasible_direction));

  *locally_infeasible = !(feasible_direction || active_trust_region);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_get_dual_estimation(SleqpIterate* iterate,
                                    void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SleqpSparseVec* vars_dual = sleqp_iterate_get_vars_dual(iterate);
  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(iterate);

  const double zero_eps = sleqp_params_get(cauchy_data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  if(vars_dual)
  {
    assert(vars_dual->dim == num_variables);

    SLEQP_CALL(sleqp_lpi_get_dual_sol(cauchy_data->lp_interface,
                                      cauchy_data->dual_values,
                                      NULL));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(vars_dual,
                                            cauchy_data->dual_values,
                                            vars_dual->dim,
                                            zero_eps));

    // Note: We rescale here since sign conventions vary...
    SLEQP_CALL(sleqp_sparse_vector_scale(vars_dual, -1.));

    for(int k = 0; k < vars_dual->nnz; ++k)
    {
      const int i = vars_dual->indices[k];

      const SLEQP_ACTIVE_STATE var_state = sleqp_working_set_get_variable_state(working_set,
                                                                                i);

      if(var_state == SLEQP_INACTIVE)
      {
        vars_dual->data[k] = 0.;
      }
      else
      {
        if(var_state == SLEQP_ACTIVE_UPPER)
        {
          sleqp_assert_is_geq(vars_dual->data[k], 0., zero_eps);
        }
        else if(var_state == SLEQP_ACTIVE_LOWER)
        {
          sleqp_assert_is_leq(vars_dual->data[k], 0., zero_eps);
        }
      }
    }
  }

  if(cons_dual)
  {
    assert(cons_dual->dim == num_constraints);

    SLEQP_CALL(sleqp_lpi_get_dual_sol(cauchy_data->lp_interface,
                                      NULL,
                                      cauchy_data->dual_values));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_dual,
                                            cauchy_data->dual_values,
                                            cons_dual->dim,
                                            zero_eps));

    // Note: We rescale here since sign conventions vary...
    SLEQP_CALL(sleqp_sparse_vector_scale(cons_dual, -1.));

    for(int k = 0; k < cons_dual->nnz; ++k)
    {
      const int i = cons_dual->indices[k];

      const SLEQP_ACTIVE_STATE var_state = sleqp_working_set_get_constraint_state(working_set,
                                                                                  i);

      if(var_state == SLEQP_INACTIVE)
      {
        cons_dual->data[k] = 0.;
      }
      else
      {
        if(var_state == SLEQP_ACTIVE_UPPER)
        {
          sleqp_assert_is_geq(cons_dual->data[k], 0., zero_eps);
        }
        else if(var_state == SLEQP_ACTIVE_LOWER)
        {
          sleqp_assert_is_leq(cons_dual->data[k], 0., zero_eps);
        }
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_get_violation(double* violation,
                              void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  SLEQP_CALL(sleqp_lpi_get_primal_sol(cauchy_data->lp_interface,
                                      NULL,
                                      cauchy_data->solution_values));

  (*violation) = 0.;

  for(int j = num_variables; j < cauchy_data->num_lp_variables; ++j)
  {
    (*violation) += cauchy_data->solution_values[j];
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_get_basis_condition(bool* exact,
                                    double* condition,
                                    void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SLEQP_CALL(sleqp_lpi_get_basis_condition(cauchy_data->lp_interface,
                                           exact,
                                           condition));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_free(void* star)
{
  CauchyData* cauchy_data = (CauchyData*) star;

  SLEQP_CALL(sleqp_sparse_vector_free(&cauchy_data->quadratic_gradient));

  sleqp_free(&cauchy_data->dual_values);
  sleqp_free(&cauchy_data->solution_values);

  sleqp_free(&cauchy_data->vars_ub);
  sleqp_free(&cauchy_data->vars_lb);

  sleqp_free(&cauchy_data->cons_ub);
  sleqp_free(&cauchy_data->cons_lb);

  sleqp_free(&cauchy_data->objective);

  SLEQP_CALL(sleqp_lpi_release(&cauchy_data->lp_interface));

  sleqp_free(&cauchy_data->cons_stats);
  sleqp_free(&cauchy_data->var_stats);

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  SLEQP_CALL(sleqp_options_release(&cauchy_data->options));
  SLEQP_CALL(sleqp_params_release(&cauchy_data->params));

  SLEQP_CALL(sleqp_problem_release(&cauchy_data->problem));

  sleqp_free(&star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_standard_cauchy_create(SleqpCauchy** star,
                                           SleqpProblem* problem,
                                           SleqpParams* params,
                                           SleqpOptions* options,
                                           SleqpLPi* lp_interface)
{
  CauchyData* cauchy_data;

  SLEQP_CALL(cauchy_data_create(&cauchy_data,
                                problem,
                                params,
                                options,
                                lp_interface));

  SleqpCauchyCallbacks callbacks = {
    .set_iterate         = standard_cauchy_set_iterate,
    .set_trust_radius    = standard_cauchy_set_trust_radius,
    .solve               = standard_cauchy_solve,
    .get_objective_value = standard_cauchy_get_objective_value,
    .get_working_set     = standard_cauchy_get_working_set,
    .get_direction       = standard_cauchy_get_direction,
    .locally_infeasible  = standard_cauchy_locally_infeasible,
    .get_dual_estimation = standard_cauchy_get_dual_estimation,
    .get_violation       = standard_cauchy_get_violation,
    .get_basis_condition = standard_cauchy_get_basis_condition,
    .free                = standard_cauchy_free
  };

  SLEQP_CALL(sleqp_cauchy_create(star,
                                 &callbacks,
                                 (void*) cauchy_data));

  return SLEQP_OKAY;
}
