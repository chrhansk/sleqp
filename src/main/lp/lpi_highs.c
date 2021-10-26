#include "lpi_highs.h"

#include <assert.h>
#include <math.h>

#include <interfaces/highs_c_api.h>

#include "cmp.h"
#include "defs.h"
#include "log.h"
#include "mem.h"

#define HIGHS_LOWER    0
#define HIGHS_BASIC    1
#define HIGHS_UPPER    2
#define HIGHS_ZERO     3
#define HIGHS_NONBASIC 4

#define HIGHS_NOTSET                  0
#define HIGHS_ERROR                   1
#define HIGHS_MODEL_ERROR             2
#define HIGHS_PRESOLVE_ERROR          3
#define HIGHS_SOLVE_ERROR             4
#define HIGHS_POSTSOLVE_ERROR         5
#define HIGHS_MODEL_EMPTY             6
#define HIGHS_OPTIMAL                 7
#define HIGHS_INFEASIBLE              8
#define HIGHS_UNBOUNDED_OR_INFEASIBLE 9
#define HIGHS_UNBOUNDED               10
#define HIGHS_OBJECTIVE_BOUND         11
#define HIGHS_OBJECTIVE_TARGET        12
#define HIGHS_TIME_LIMIT              13
#define HIGHS_ITERATION_LIMIT         14
#define HIGHS_UNKNOWN                 15

typedef struct SleqpLpiHIGHS
{
  void* highs;

  SLEQP_LPI_STATUS status;

  int num_cols;
  int num_rows;

  int num_bases;
  int** vbases;
  int** cbases;

  int* col_basis;
  int* row_basis;

  double* cols_primal_dummysol;
  double* rows_primal_dummysol;
  double* cols_dual_dummysol;
  double* rows_dual_dummysol;

} SleqpLpiHIGHS;

#define SLEQP_HIGHS_CALL(x, highs)                \
  do                                              \
  {                                               \
    int highs_ret_status = (x);                   \
                                                  \
    if(highs_ret_status != HighsStatuskOk)        \
    {                                             \
      sleqp_log_error("Caught HiGHS error <%d>",  \
                      highs_ret_status);          \
                                                  \
      return SLEQP_INTERNAL_ERROR;                \
    }                                             \
  } while(0)

static SLEQP_RETCODE highs_create_problem(void** star,
                                          int num_cols,
                                          int num_rows,
                                          SleqpParams* params,
                                          SleqpOptions* options)
{
  SleqpLpiHIGHS* lp_interface = NULL;

  const int sense_min = 1;
  const int format_sparse_columns = 1;
  const int dummy_nnz = 0;
  const double objective_offset = 0.;

  double *zero_vector_columns;
  double *zero_vector_rows;
  int *zero_matrix_colptr;

  SLEQP_CALL(sleqp_malloc(&lp_interface));

  *lp_interface = (SleqpLpiHIGHS){0};

  *star = lp_interface;

  lp_interface->num_cols = num_cols;
  lp_interface->num_rows = num_rows;
  lp_interface->status = SLEQP_LPI_STATUS_UNKNOWN;

  SLEQP_CALL(sleqp_alloc_array(&lp_interface->col_basis, num_cols));
  SLEQP_CALL(sleqp_alloc_array(&lp_interface->row_basis, num_rows));
  SLEQP_CALL(sleqp_alloc_array(&lp_interface->cols_primal_dummysol,
                               lp_interface->num_cols));
  SLEQP_CALL(sleqp_alloc_array(&lp_interface->rows_primal_dummysol,
                               lp_interface->num_rows));
  SLEQP_CALL(sleqp_alloc_array(&lp_interface->cols_dual_dummysol,
                               lp_interface->num_cols));
  SLEQP_CALL(sleqp_alloc_array(&lp_interface->rows_dual_dummysol,
                               lp_interface->num_rows));


  lp_interface->highs = Highs_create();
  void* highs = lp_interface->highs;

  if(sleqp_log_level() < SLEQP_LOG_DEBUG)
  {
    SLEQP_HIGHS_CALL(Highs_setBoolOptionValue(highs,
                                              "output_flag",
                                              false), highs);
  }

  {
    const int num_threads = sleqp_options_get_int(options,
                                                  SLEQP_OPTION_INT_NUM_THREADS);

    if(num_threads != SLEQP_NONE)
    {
      SLEQP_HIGHS_CALL(Highs_setIntOptionValue(highs,
                                               "highs_max_threads",
                                               num_threads), highs);
    }
  }

  SLEQP_HIGHS_CALL(Highs_setDoubleOptionValue(highs,
                                              "infinite_cost",
                                              sleqp_infinity()), highs);
  SLEQP_HIGHS_CALL(Highs_setDoubleOptionValue(highs,
                                              "infinite_bound",
                                              sleqp_infinity()), highs);

  // setting these option causes HiGHS to fail for the
  // constrained test example without quadratic model while the defaults work

  //const double feas_eps = sleqp_params_get(params,
  //                                         SLEQP_PARAM_FEASIBILITY_TOL);

  //const double stat_eps = sleqp_params_get(params,
  //                                         SLEQP_PARAM_STATIONARITY_TOL);

  //SLEQP_HIGHS_CALL(Highs_setDoubleOptionValue(highs,
  //                                            "primal_feasibility_tolerance",
  //                                            feas_eps), highs);

  //SLEQP_HIGHS_CALL(Highs_setDoubleOptionValue(highs,
  //                                            "dual_feasibility_tolerance",
  //                                            stat_eps), highs);

  // allocate dummy vectors to pass an empty linear program of appropriate size

  SLEQP_CALL(sleqp_alloc_array(&zero_vector_columns, lp_interface->num_cols));
  SLEQP_CALL(sleqp_alloc_array(&zero_vector_rows, lp_interface->num_rows));
  SLEQP_CALL(sleqp_alloc_array(&zero_matrix_colptr, lp_interface->num_cols+1));

  for(int i = 0; i < num_cols; ++i)
  {
    zero_vector_columns[i] = 0.;
  }
  for(int i = 0; i < num_rows; ++i)
  {
    zero_vector_rows[i] = 0.;
  }
  for(int i = 0; i < num_cols + 1; ++i)
  {
    zero_matrix_colptr[i] = 0;
  }

  SLEQP_HIGHS_CALL(Highs_passLp(highs,
                                lp_interface->num_cols,
                                lp_interface->num_rows,
                                dummy_nnz,
                                format_sparse_columns,
                                sense_min,
                                objective_offset,
                                zero_vector_columns,
                                zero_vector_columns,
                                zero_vector_columns,
                                zero_vector_rows,
                                zero_vector_rows,
                                zero_matrix_colptr,
                                NULL,
                                NULL), highs);

  sleqp_free(&zero_vector_columns);
  sleqp_free(&zero_vector_rows);
  sleqp_free(&zero_matrix_colptr);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_write(void* lp_data, const char* filename)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  SLEQP_HIGHS_CALL(Highs_writeModel(highs, filename), highs);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_solve(void* lp_data,
                                 int num_cols,
                                 int num_rows,
                                 double time_limit)
{
  int model_status;
  int scaled_model_status;

  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  if(time_limit != SLEQP_NONE)
  {
    SLEQP_HIGHS_CALL(Highs_setDoubleOptionValue(highs,
                                                "time_limit",
                                                time_limit), highs);
  }

  SLEQP_HIGHS_CALL(Highs_run(highs), highs);

  model_status = Highs_getModelStatus(highs);
  scaled_model_status = Highs_getScaledModelStatus(highs);

  switch(model_status)
  {
  case HIGHS_NOTSET:
    // This might be related to https://github.com/ERGO-Code/HiGHS/issues/494!
    sleqp_log_warn("Received HiGHS status NOTSET for original model after run");
    switch(scaled_model_status)
    {
    case HIGHS_OPTIMAL:
      lp_interface->status = SLEQP_LPI_STATUS_OPTIMAL;
      break;
    case HIGHS_INFEASIBLE:
      lp_interface->status = SLEQP_LPI_STATUS_INF;
      break;
    case HIGHS_UNBOUNDED_OR_INFEASIBLE:
      lp_interface->status = SLEQP_LPI_STATUS_INF_OR_UNBOUNDED;
      break;
    case HIGHS_UNBOUNDED:
      lp_interface->status = SLEQP_LPI_STATUS_UNBOUNDED;
      break;
    case HIGHS_TIME_LIMIT:
      lp_interface->status = SLEQP_LPI_STATUS_UNKNOWN;
      return SLEQP_ABORT_TIME;
      break;
    default:
      sleqp_log_error("Invalid HiGHS status for scaled model: %d", model_status);
      lp_interface->status = SLEQP_LPI_STATUS_UNKNOWN;
      return SLEQP_INTERNAL_ERROR;
    }
    break;
    //case HIGHS_ERROR:
    //case HIGHS_MODEL_ERROR:
    //case HIGHS_PRESOLVE_ERROR:
    //case HIGHS_SOLVE_ERROR:
    //case HIGHS_POSTSOLVE_ERROR:
    //case HIGHS_MODEL_ERROR:
  case HIGHS_OPTIMAL:
    lp_interface->status = SLEQP_LPI_STATUS_OPTIMAL;
    break;
  case HIGHS_INFEASIBLE:
    lp_interface->status = SLEQP_LPI_STATUS_INF;
    break;
  case HIGHS_UNBOUNDED_OR_INFEASIBLE:
    lp_interface->status = SLEQP_LPI_STATUS_INF_OR_UNBOUNDED;
    break;
  case HIGHS_UNBOUNDED:
    lp_interface->status = SLEQP_LPI_STATUS_UNBOUNDED;
    break;
    //case HIGHS_OBJECTIVE_BOUND:
    //case HIGHS_OBJECTIVE_TARGET:
  case HIGHS_TIME_LIMIT:
    lp_interface->status = SLEQP_LPI_STATUS_UNKNOWN;
    return SLEQP_ABORT_TIME;
    break;
  case HIGHS_ITERATION_LIMIT:  // fallthrough
  case HIGHS_UNKNOWN:          // fallthrough
  default:
    sleqp_log_error("Invalid HiGHS status: %d", model_status);
    lp_interface->status = SLEQP_LPI_STATUS_UNKNOWN;
    return SLEQP_INTERNAL_ERROR;
  }

  return SLEQP_OKAY;
}

static SLEQP_LPI_STATUS highs_get_status(void* lp_data)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;

  return lp_interface->status;
}

static double adjust_neg_inf(double value)
{
  return sleqp_is_infinite(-value) ? -INFINITY : value;
}

static double adjust_pos_inf(double value)
{
  return sleqp_is_infinite(value) ? INFINITY : value;
}

static SLEQP_RETCODE highs_set_bounds(void* lp_data,
                                      int num_cols,
                                      int num_rows,
                                      double* cons_lb,
                                      double* cons_ub,
                                      double* vars_lb,
                                      double* vars_ub)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  for(int i = 0; i < num_cols; ++i)
  {
    SLEQP_HIGHS_CALL(Highs_changeColBounds(highs,
                                           i,
                                           adjust_neg_inf(vars_lb[i]),
                                           adjust_pos_inf(vars_ub[i])), highs);
  }

  for(int i = 0; i < num_rows; ++i)
  {
    SLEQP_HIGHS_CALL(Highs_changeRowBounds(highs,
                                           i,
                                           adjust_neg_inf(cons_lb[i]),
                                           adjust_pos_inf(cons_ub[i])), highs);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_set_coefficients(void* lp_data,
                                            int num_cols,
                                            int num_rows,
                                            SleqpSparseMatrix* coeff_matrix)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

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

      SLEQP_HIGHS_CALL(Highs_changeCoeff(highs,
                                         row,
                                         col,
                                         entry), highs);
    }
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE highs_set_objective(void* lp_data,
                                         int num_cols,
                                         int num_rows,
                                         double* objective)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  for(int i = 0; i < num_cols; ++i)
  {

    SLEQP_HIGHS_CALL(Highs_changeColCost(highs,
                                         i,
                                         objective[i]), highs);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
highs_reserve_bases(SleqpLpiHIGHS* lp_interface,
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
    SLEQP_CALL(sleqp_alloc_array(&lp_interface->vbases[j], lp_interface->num_cols));
    SLEQP_CALL(sleqp_alloc_array(&lp_interface->cbases[j], lp_interface->num_rows));
  }

  lp_interface->num_bases = size;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_save_basis(void* lp_data,
                                      int index)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  assert(index >= 0);

  SLEQP_CALL(highs_reserve_bases(lp_interface, index + 1));

  SLEQP_HIGHS_CALL(Highs_getBasis(highs,
                                  lp_interface->vbases[index],
                                  lp_interface->cbases[index]), highs);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_restore_basis(void* lp_data,
                                         int index)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  assert(index >= 0);
  assert(index < lp_interface->num_bases);

  SLEQP_HIGHS_CALL(Highs_setBasis(highs,
                                  lp_interface->vbases[index],
                                  lp_interface->cbases[index]), highs);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_get_primal_sol(void* lp_data,
                                          int num_cols,
                                          int num_rows,
                                          double* objective_value,
                                          double* solution_values)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  if(objective_value)
  {
    *objective_value = Highs_getObjectiveValue(highs);
  }

  if(solution_values)
  {
    SLEQP_HIGHS_CALL(Highs_getSolution(highs,
                                       solution_values,
                                       lp_interface->cols_dual_dummysol,
                                       lp_interface->rows_primal_dummysol,
                                       lp_interface->rows_dual_dummysol), highs);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_get_dual_sol(void* lp_data,
                                        int num_cols,
                                        int num_rows,
                                        double* vars_dual,
                                        double* cons_dual)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  assert(lp_interface->num_cols == num_cols);

  if(!vars_dual)
  {
    vars_dual = lp_interface->cols_dual_dummysol;
  }

  if(!cons_dual)
  {
    cons_dual = lp_interface->rows_dual_dummysol;
  }

  SLEQP_HIGHS_CALL(Highs_getSolution(highs,
                                     lp_interface->cols_primal_dummysol,
                                     vars_dual,
                                     cons_dual,
                                     cons_dual), highs);

  return SLEQP_OKAY;
}

static SLEQP_BASESTAT basestat_for(int status)
{
  switch (status)
  {
  case HIGHS_BASIC:
    return SLEQP_BASESTAT_BASIC;
  case HIGHS_LOWER:
    return SLEQP_BASESTAT_LOWER;
  case HIGHS_UPPER:
    return SLEQP_BASESTAT_UPPER;
  case HIGHS_ZERO:
    return SLEQP_BASESTAT_ZERO;
  case HIGHS_NONBASIC:
    sleqp_log_error("Encountered an unspecific non-basic variable");
    break;
  default:
    break;
  }

  // Invalid basis status reported
  assert(false);

  return SLEQP_BASESTAT_ZERO;
}

static SLEQP_RETCODE highs_get_varstats(void* lp_data,
                                        int num_cols,
                                        int num_rows,
                                        SLEQP_BASESTAT* variable_stats)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  SLEQP_HIGHS_CALL(Highs_getBasis(highs,
                                  lp_interface->col_basis,
                                  lp_interface->row_basis), highs);

  for(int j = 0; j < num_cols; ++j)
  {
    variable_stats[j] = basestat_for(lp_interface->col_basis[j]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_get_consstats(void* lp_data,
                                         int num_cols,
                                         int num_rows,
                                         SLEQP_BASESTAT* constraint_stats)
{
  SleqpLpiHIGHS* lp_interface = (SleqpLpiHIGHS*) lp_data;
  void* highs = lp_interface->highs;

  SLEQP_HIGHS_CALL(Highs_getBasis(highs,
                                  lp_interface->col_basis,
                                  lp_interface->row_basis), highs);

  for(int j = 0; j < num_rows; ++j)
  {
    constraint_stats[j] = basestat_for(lp_interface->row_basis[j]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_get_basis_condition(void *lp_data,
                                               bool* exact,
                                               double* condition)
{
  *exact = false;
  *condition = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE highs_free(void** star)
{
  SleqpLpiHIGHS* lp_interface = *star;
  void* highs = lp_interface->highs;

  if(!lp_interface)
  {
    return SLEQP_OKAY;
  }

  if(lp_interface->highs)
  {
    Highs_destroy(highs);
  }

  for(int i = 0; i < lp_interface->num_bases; ++i)
  {
    sleqp_free(&lp_interface->vbases[i]);
    sleqp_free(&lp_interface->cbases[i]);
  }

  sleqp_free(&lp_interface->vbases);
  sleqp_free(&lp_interface->cbases);

  sleqp_free(&lp_interface->row_basis);
  sleqp_free(&lp_interface->col_basis);
  sleqp_free(&lp_interface->cols_primal_dummysol);
  sleqp_free(&lp_interface->rows_primal_dummysol);
  sleqp_free(&lp_interface->cols_dual_dummysol);
  sleqp_free(&lp_interface->rows_dual_dummysol);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lpi_highs_create_interface(SleqpLPi** lp_star,
                                               int num_cols,
                                               int num_rows,
                                               SleqpParams* params,
                                               SleqpOptions* options)
{
  SleqpLPiCallbacks callbacks = {
    .create_problem      = highs_create_problem,
    .solve               = highs_solve,
    .get_status          = highs_get_status,
    .set_bounds          = highs_set_bounds,
    .set_coefficients    = highs_set_coefficients,
    .set_objective       = highs_set_objective,
    .save_basis          = highs_save_basis,
    .restore_basis       = highs_restore_basis,
    .get_primal_sol      = highs_get_primal_sol,
    .get_dual_sol        = highs_get_dual_sol,
    .get_varstats        = highs_get_varstats,
    .get_consstats       = highs_get_consstats,
    .get_basis_condition = highs_get_basis_condition,
    .free_problem        = highs_free
  };

  return sleqp_lpi_create_interface(lp_star,
                                    SLEQP_LP_SOLVER_HIGHS_NAME,
                                    SLEQP_LP_SOLVER_HIGHS_VERSION,
                                    num_cols,
                                    num_rows,
                                    params,
                                    options,
                                    &callbacks);
}


SLEQP_RETCODE sleqp_lpi_create_default_interface(SleqpLPi** lp_interface,
                                                 int num_variables,
                                                 int num_constraints,
                                                 SleqpParams* params,
                                                 SleqpOptions* options)
{
  SLEQP_CALL(sleqp_lpi_highs_create_interface(lp_interface,
                                              num_variables,
                                              num_constraints,
                                              params,
                                              options));

  return SLEQP_OKAY;
}
