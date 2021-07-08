#include "lpi_gurobi.h"

#include <assert.h>
#include <gurobi_c.h>

#include "cmp.h"
#include "log.h"
#include "mem.h"

/*
 * Note: We handle the constraints \f$ l \leq Ax \leq u \f$
 * via \f$ A x - y = 0, l \leq y \leq u \f$
 */

#define GRB_BASIC 0

typedef struct SleqpLpiGRB
{
  GRBenv* env;
  GRBmodel* model;

  int num_cols;
  int num_rows;

  int num_lp_cols;

  int num_bases;
  int** vbases;
  int** cbases;

  int* slack_basis;
  int* col_basis;
  int* row_basis;

} SleqpLpiGRB;

#define SLEQP_GRB_CALL(x, env)                          \
  do                                                    \
  {                                                     \
    int grb_ret_status = (x);                           \
                                                        \
    if(grb_ret_status)                                  \
    {                                                   \
      const char* error_string = GRBgeterrormsg(env);   \
                                                        \
      sleqp_log_error("Caught Gurobi error <%d> (%s)",  \
                      grb_ret_status,                   \
                      error_string);                    \
                                                        \
      switch(grb_ret_status)                            \
      {                                                 \
      case GRB_ERROR_OUT_OF_MEMORY:                     \
        return SLEQP_NOMEM;                             \
      default:                                          \
        return SLEQP_INTERNAL_ERROR;                    \
      }                                                 \
    }                                                   \
  } while(0)

static SLEQP_RETCODE gurobi_create_problem(void** star,
                                           int num_cols,
                                           int num_rows,
                                           SleqpParams* params)
{
  SleqpLpiGRB* lp_interface = NULL;

  SLEQP_CALL(sleqp_malloc(&lp_interface));

  *lp_interface = (SleqpLpiGRB){0};

  *star = lp_interface;

  lp_interface->num_cols = num_cols;
  lp_interface->num_rows = num_rows;
  lp_interface->num_lp_cols = num_rows + num_cols;

  SLEQP_CALL(sleqp_alloc_array(&lp_interface->col_basis, num_cols));
  SLEQP_CALL(sleqp_alloc_array(&lp_interface->slack_basis, num_rows));

  SLEQP_CALL(sleqp_alloc_array(&lp_interface->row_basis, num_rows));

  int err = GRBloadenv(&lp_interface->env, NULL);

  GRBenv* env = lp_interface->env;

  if(err || env == NULL)
  {
    sleqp_log_error("Failed to create Gurobi environment");
    return SLEQP_INTERNAL_ERROR;
  }

  if(sleqp_log_level() < SLEQP_LOG_DEBUG)
  {
    SLEQP_GRB_CALL(GRBsetintparam(env,
                                  GRB_INT_PAR_OUTPUTFLAG,
                                  0), env);
  }

  SLEQP_GRB_CALL(GRBnewmodel(env,
                             &lp_interface->model,
                             "",
                             lp_interface->num_lp_cols,
                             NULL,
                             NULL,
                             NULL,
                             NULL,
                             NULL), env);

  GRBmodel* model = lp_interface->model;

  SLEQP_GRB_CALL(GRBsetintattr(model,
                               GRB_INT_ATTR_MODELSENSE,
                               GRB_MINIMIZE), env);

  for(int i = 0; i < num_rows; ++i)
  {

    SLEQP_GRB_CALL(GRBaddconstr(model,
                                0,
                                NULL,
                                NULL,
                                GRB_EQUAL,
                                0.,
                                NULL), env);
  }

  double neg_one = -1.;

  for(int i = 0; i < num_rows; ++i)
  {
    int c = num_cols + i;

    SLEQP_GRB_CALL(GRBchgcoeffs(model,
                                1,
                                &i,
                                &c,
                                &neg_one), env);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_solve(void* lp_data,
                                  int num_cols,
                                  int num_rows,
                                  double time_limit)
{
  SleqpLpiGRB* lp_interface = lp_data;
  int sol_stat;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  // SLEQP_GRB_CALL(GRBwrite(model, "file.lp"), env);

  if(time_limit != SLEQP_NONE)
  {
    GRBenv* model_env = GRBgetenv(model);

    SLEQP_GRB_CALL(GRBsetdblparam(model_env,
                                  GRB_DBL_PAR_TIMELIMIT,
                                  time_limit),
                   model_env);
  }

  SLEQP_GRB_CALL(GRBoptimize(model), env);

  SLEQP_GRB_CALL(GRBgetintattr(model, GRB_INT_ATTR_STATUS, &sol_stat), env);

  assert(sol_stat == GRB_OPTIMAL);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_set_bounds(void* lp_data,
                                       int num_cols,
                                       int num_rows,
                                       double* cons_lb,
                                       double* cons_ub,
                                       double* vars_lb,
                                       double* vars_ub)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  SLEQP_GRB_CALL(GRBsetdblattrarray(model,
                                    GRB_DBL_ATTR_LB,
                                    0,
                                    num_cols,
                                    vars_lb), env);

  SLEQP_GRB_CALL(GRBsetdblattrarray(model,
                                    GRB_DBL_ATTR_UB,
                                    0,
                                    num_cols,
                                    vars_ub), env);

  SLEQP_GRB_CALL(GRBsetdblattrarray(model,
                                    GRB_DBL_ATTR_LB,
                                    num_cols,
                                    num_rows,
                                    cons_lb), env);

  SLEQP_GRB_CALL(GRBsetdblattrarray(model,
                                    GRB_DBL_ATTR_UB,
                                    num_cols,
                                    num_rows,
                                    cons_ub), env);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_set_coefficients(void* lp_data,
                                             int num_cols,
                                             int num_rows,
                                             SleqpSparseMatrix* coeff_matrix)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  assert(sleqp_sparse_matrix_get_num_rows(coeff_matrix) == num_rows);
  assert(sleqp_sparse_matrix_get_num_cols(coeff_matrix) == num_cols);

  const int* coeff_matrix_cols = sleqp_sparse_matrix_get_cols(coeff_matrix);
  const int* coeff_matrix_rows = sleqp_sparse_matrix_get_rows(coeff_matrix);
  double* coeff_matrix_data = sleqp_sparse_matrix_get_data(coeff_matrix);

  for(int col = 0; col < num_cols; ++col)
  {
    for(int k = coeff_matrix_cols[col]; k < coeff_matrix_cols[col +1]; ++k)
    {
      int row = coeff_matrix_rows[k];
      double entry = coeff_matrix_data[k];

      SLEQP_GRB_CALL(GRBchgcoeffs(model,
                                  1,
                                  &row,
                                  &col,
                                  &entry), env);
    }
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE gurobi_set_objective(void* lp_data,
                                          int num_cols,
                                          int num_rows,
                                          double* objective)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  SLEQP_GRB_CALL(GRBsetdblattrarray(model,
                                    GRB_DBL_ATTR_OBJ,
                                    0,
                                    num_cols,
                                    objective), env);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
gurobi_reserve_bases(SleqpLpiGRB* lp_interface,
                     int size)
{
  if(size <= lp_interface->num_bases)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_realloc(&lp_interface->vbases, size));
  SLEQP_CALL(sleqp_realloc(&lp_interface->cbases, size));

  for(int j = lp_interface->num_bases; j < size; ++j)
  {
    SLEQP_CALL(sleqp_alloc_array(&lp_interface->vbases[j], lp_interface->num_lp_cols));
    SLEQP_CALL(sleqp_alloc_array(&lp_interface->cbases[j], lp_interface->num_rows));
  }

  lp_interface->num_bases = size;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_save_basis(void* lp_data,
                                       int index)
{
  SleqpLpiGRB* lp_interface = lp_data;

  assert(index >= 0);

  SLEQP_CALL(gurobi_reserve_bases(lp_interface, index + 1));

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  SLEQP_GRB_CALL(GRBgetintattrarray(model,
                                    GRB_INT_ATTR_VBASIS,
                                    0,
                                    lp_interface->num_lp_cols,
                                    lp_interface->vbases[index]), env);

  SLEQP_GRB_CALL(GRBgetintattrarray(model,
                                    GRB_INT_ATTR_CBASIS,
                                    0,
                                    lp_interface->num_rows,
                                    lp_interface->cbases[index]), env);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_restore_basis(void* lp_data,
                                          int index)
{
  SleqpLpiGRB* lp_interface = lp_data;

  assert(index >= 0);
  assert(index < lp_interface->num_bases);

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  SLEQP_GRB_CALL(GRBsetintattrarray(model,
                                    GRB_INT_ATTR_VBASIS,
                                    0,
                                    lp_interface->num_lp_cols,
                                    lp_interface->vbases[index]), env);

  SLEQP_GRB_CALL(GRBsetintattrarray(model,
                                    GRB_INT_ATTR_CBASIS,
                                    0,
                                    lp_interface->num_rows,
                                    lp_interface->cbases[index]), env);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_get_primal_sol(void* lp_data,
                                           int num_cols,
                                           int num_rows,
                                           double* objective_value,
                                           double* solution_values)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  if(objective_value)
  {
    SLEQP_GRB_CALL(GRBgetdblattr(model,
                                 GRB_DBL_ATTR_OBJVAL,
                                 objective_value), env);
  }

  if(solution_values)
  {
    SLEQP_GRB_CALL(GRBgetdblattrarray(model,
                                      GRB_DBL_ATTR_X,
                                      0,
                                      num_cols,
                                      solution_values), env);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_get_dual_sol(void* lp_data,
                                         int num_cols,
                                         int num_rows,
                                         double* vars_dual,
                                         double* cons_dual)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  assert(lp_interface->num_lp_cols == num_rows + num_cols);

  if(cons_dual)
  {
    SLEQP_GRB_CALL(GRBgetdblattrarray(model,
                                      GRB_DBL_ATTR_PI,
                                      0,
                                      num_rows,
                                      cons_dual), env);
  }

  if(vars_dual)
  {
    SLEQP_GRB_CALL(GRBgetdblattrarray(model,
                                      GRB_DBL_ATTR_RC,
                                      0,
                                      num_cols,
                                      vars_dual), env);
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE gurobi_get_varstats(void* lp_data,
                                         int num_cols,
                                         int num_rows,
                                         SLEQP_BASESTAT* variable_stats)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  SLEQP_GRB_CALL(GRBgetintattrarray(model,
                                    GRB_INT_ATTR_VBASIS,
                                    0,
                                    num_cols,
                                    lp_interface->col_basis), env);

  for(int j = 0; j < num_cols; ++j)
  {
    switch(lp_interface->col_basis[j])
    {
    case GRB_BASIC:
      variable_stats[j] = SLEQP_BASESTAT_BASIC;
      break;
    case GRB_NONBASIC_LOWER:
      variable_stats[j] = SLEQP_BASESTAT_LOWER;
      break;
    case GRB_NONBASIC_UPPER:
      variable_stats[j] = SLEQP_BASESTAT_UPPER;
      break;
    case GRB_SUPERBASIC:
      sleqp_log_error("Encountered a super-basic variable");
    default:
      assert(false);
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_get_consstats(void* lp_data,
                                          int num_cols,
                                          int num_rows,
                                          SLEQP_BASESTAT* constraint_stats)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  SLEQP_GRB_CALL(GRBgetintattrarray(model,
                                    GRB_INT_ATTR_VBASIS,
                                    num_cols,
                                    num_rows,
                                    lp_interface->slack_basis), env);

  SLEQP_GRB_CALL(GRBgetintattrarray(model,
                                    GRB_INT_ATTR_CBASIS,
                                    0,
                                    num_rows,
                                    lp_interface->row_basis), env);

  for(int i = 0; i < num_rows; ++i)
  {
    if(lp_interface->row_basis[i] == GRB_BASIC)
    {
      constraint_stats[i] = SLEQP_BASESTAT_BASIC;
      continue;
    }

    switch(lp_interface->slack_basis[i])
    {
    case GRB_BASIC:
      constraint_stats[i] = SLEQP_BASESTAT_BASIC;
      break;
    case GRB_NONBASIC_LOWER:
      constraint_stats[i] = SLEQP_BASESTAT_LOWER;
      break;
    case GRB_NONBASIC_UPPER:
      constraint_stats[i] = SLEQP_BASESTAT_UPPER;
      break;
    case GRB_SUPERBASIC:
      sleqp_log_error("Encountered a super-basic constraint");
      return SLEQP_INTERNAL_ERROR;
    default:
      assert(false);
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_get_basis_condition(void *lp_data,
                                                bool* exact,
                                                double* condition)
{
  SleqpLpiGRB* lp_interface = lp_data;

  GRBenv* env = lp_interface->env;
  GRBmodel* model = lp_interface->model;

  if(*exact)
  {
    SLEQP_GRB_CALL(GRBgetdblattr(model, GRB_DBL_ATTR_KAPPA_EXACT, condition), env);
  }
  else
  {
    SLEQP_GRB_CALL(GRBgetdblattr(model, GRB_DBL_ATTR_KAPPA, condition), env);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE gurobi_free(void** star)
{
  SleqpLpiGRB* lp_interface = *star;

  if(!lp_interface)
  {
    return SLEQP_OKAY;
  }

  if(lp_interface->model)
  {
    GRBfreemodel(lp_interface->model);
  }
  for(int i = 0; i < lp_interface->num_bases; ++i)
  {
    sleqp_free(&lp_interface->vbases[i]);
    sleqp_free(&lp_interface->cbases[i]);
  }

  sleqp_free(&lp_interface->vbases);
  sleqp_free(&lp_interface->cbases);

  /* Free environment */
  if(lp_interface->env)
  {
    GRBfreeenv(lp_interface->env);
  }

  sleqp_free(&lp_interface->row_basis);

  sleqp_free(&lp_interface->slack_basis);
  sleqp_free(&lp_interface->col_basis);


  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lpi_gurobi_create_interface(SleqpLPi** lp_star,
                                                int num_cols,
                                                int num_rows,
                                                SleqpParams* params)
{
  SleqpLPiCallbacks callbacks = {
    .create_problem = gurobi_create_problem,
    .solve = gurobi_solve,
    .set_bounds = gurobi_set_bounds,
    .set_coefficients = gurobi_set_coefficients,
    .set_objective = gurobi_set_objective,
    .save_basis = gurobi_save_basis,
    .restore_basis = gurobi_restore_basis,
    .get_primal_sol = gurobi_get_primal_sol,
    .get_dual_sol = gurobi_get_dual_sol,
    .get_varstats = gurobi_get_varstats,
    .get_consstats = gurobi_get_consstats,
    .get_basis_condition = gurobi_get_basis_condition,
    .free_problem = gurobi_free
  };

  return sleqp_lpi_create_interface(lp_star,
                                    num_cols,
                                    num_rows,
                                    params,
                                    &callbacks);
}


SLEQP_RETCODE sleqp_lpi_create_default_interface(SleqpLPi** lp_interface,
                                                 int num_variables,
                                                 int num_constraints,
                                                 SleqpParams* params)
{
  SLEQP_CALL(sleqp_lpi_gurobi_create_interface(lp_interface,
                                               num_variables,
                                               num_constraints,
                                               params));

  return SLEQP_OKAY;
}
