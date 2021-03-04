#include "sleqp_lpi_soplex.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include <iostream>
#include <soplex.h>

struct SleqpLpiSoplex
{
  soplex::SoPlex* soplex;

  soplex::SPxSolver::VarStatus* basis_rows;
  soplex::SPxSolver::VarStatus* basis_cols;

  int num_cols;
  int num_rows;
};

static SLEQP_RETCODE soplex_create_problem(void** lp_data,
                                           int num_cols,
                                           int num_rows,
                                           SleqpParams* params)
{
  SleqpLpiSoplex* spx = new SleqpLpiSoplex;

  *spx = (SleqpLpiSoplex){0};

  spx->soplex = new soplex::SoPlex();
  soplex::SoPlex& soplex = *(spx->soplex);

  spx->basis_rows = new soplex::SPxSolver::VarStatus[num_rows];
  spx->basis_cols = new soplex::SPxSolver::VarStatus[num_cols];

  soplex::SPxOut spxout;

  spxout.setStream(soplex::SPxOut::INFO1, std::cerr);
  spxout.setStream(soplex::SPxOut::INFO2, std::cerr);
  spxout.setStream(soplex::SPxOut::INFO3, std::cerr);
  spxout.setStream(soplex::SPxOut::DEBUG, std::cerr);
  spxout.setStream(soplex::SPxOut::ERROR, std::cerr);
  spxout.setStream(soplex::SPxOut::WARNING, std::cerr);

  soplex.spxout = spxout;

  const double zero_eps = sleqp_params_get(params,
                                           SLEQP_PARAM_ZERO_EPS);

  soplex.setRealParam(soplex::SoPlex::EPSILON_ZERO, zero_eps);

  if(sleqp_log_level() >= SLEQP_LOG_DEBUG)
  {
    soplex.setIntParam(soplex::SoPlex::VERBOSITY,
                       soplex::SoPlex::VERBOSITY_HIGH);
  }
  else if(sleqp_log_level() >= SLEQP_LOG_INFO)
  {
    soplex.setIntParam(soplex::SoPlex::VERBOSITY,
                       soplex::SoPlex::VERBOSITY_NORMAL);
  }
  else if(sleqp_log_level() >= SLEQP_LOG_WARN)
  {
    soplex.setIntParam(soplex::SoPlex::VERBOSITY,
                       soplex::SoPlex::VERBOSITY_WARNING);
  }
  else if(sleqp_log_level() >= SLEQP_LOG_ERROR)
  {
    soplex.setIntParam(soplex::SoPlex::VERBOSITY,
                       soplex::SoPlex::VERBOSITY_ERROR);
  }

  spx->num_cols = num_cols;
  spx->num_rows = num_rows;

  soplex.setIntParam(soplex::SoPlex::OBJSENSE,
                     soplex::SoPlex::OBJSENSE_MINIMIZE);

  // add dummy (empty) rows / cols
  soplex::DSVector vec(0);

  double inf = soplex::infinity;

  {
    soplex::LPColSetReal cols(num_cols, 0);

    for(int j = 0; j < num_cols; ++j)
    {
      cols.add(soplex::LPCol(0., vec, inf, -inf));
    }

    soplex.addColsReal(cols);
  }

  {
    soplex::LPRowSetReal rows(num_rows ,0);

    for(int i = 0; i < num_rows; ++i)
    {
      rows.add(soplex::LPRow(-inf, vec, inf));
    }

    soplex.addRowsReal(rows);
  }

  *lp_data = spx;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_solve(void* lp_data,
                                  int num_cols,
                                  int num_rows,
                                  double time_limit)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  if(time_limit != SLEQP_NONE)
  {
    soplex.setRealParam(soplex::SoPlex::TIMELIMIT, time_limit);
  }

  soplex::SPxSolver::Status status = soplex.optimize();

  // retry the solve from scratch
  if(status == soplex::SPxSolver::SINGULAR)
  {
    sleqp_log_warn("Initial basis is singular, resolving from scratch");
    soplex.clearBasis();
    status = soplex.optimize();
  }

  //soplex.writeFileReal("test.lp");

  /*
  if(status != soplex::SPxSolver::OPTIMAL)
  {
    bool do_resolve = true;

    bool is_primal_feasible = soplex.isPrimalFeasible();
    bool is_dual_feasible = soplex.isDualFeasible();

    sleqp_log_debug("Solution of LP not optimal (pfeas=%d, dfeas=%d)",
                    is_primal_feasible,
                    is_dual_feasible);


    while (do_resolve)
    {
      do_resolve = false;

      is_primal_feasible = soplex.isPrimalFeasible();
      is_dual_feasible = soplex.isDualFeasible();

      double feas_tol = soplex.realParam(soplex::SoPlex::FEASTOL);
      double opt_tol = soplex.realParam(soplex::SoPlex::OPTTOL);

      if(!is_primal_feasible && !sleqp_is_zero(feas_tol, spx->eps))
      {
        sleqp_log_debug("Solving again with higher feasibility tolerance");

        feas_tol *= 1e-3;

        feas_tol = SLEQP_MAX(feas_tol, spx->eps);

        soplex.setRealParam(soplex::SoPlex::FEASTOL, feas_tol);

        do_resolve = true;
      }
      else if(!is_dual_feasible && !sleqp_is_zero(opt_tol, spx->eps))
      {
        sleqp_log_debug("Solving again with higher optimality tolerance");

        opt_tol *= 1e-3;

        opt_tol = SLEQP_MAX(opt_tol, spx->eps);

        soplex.setRealParam(soplex::SoPlex::OPTTOL, opt_tol);

        do_resolve = true;
      }

      if(do_resolve)
      {
        status = soplex.optimize();
      }
    }

    assert(status == soplex::SPxSolver::OPTIMAL);
  }
  */

  if(!(status == soplex::SPxSolver::OPTIMAL &&
       soplex.hasBasis()))
  {
    sleqp_log_error("Failed to solve LP");
    return SLEQP_INTERNAL_ERROR;
  }

  return SLEQP_OKAY;
}

static double adjust_inf(double value)
{
  if(sleqp_is_inf(value))
  {
    return soplex::infinity;
  }
  else if(sleqp_is_inf(-value))
  {
    return -(soplex::infinity);
  }

  return value;
}

static SLEQP_RETCODE soplex_set_bounds(void* lp_data,
                                       int num_cols,
                                       int num_rows,
                                       double* cons_lb,
                                       double* cons_ub,
                                       double* vars_lb,
                                       double* vars_ub)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int i = 0; i < num_rows; ++i)
  {
    assert(cons_lb[i] <= cons_ub[i]);
    soplex.changeRangeReal(i, adjust_inf(cons_lb[i]), adjust_inf(cons_ub[i]));
  }

  for(int j = 0; j < num_cols; ++j)
  {
    assert(vars_lb[j] <= vars_ub[j]);
    soplex.changeBoundsReal(j, adjust_inf(vars_lb[j]), adjust_inf(vars_ub[j]));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_set_coefficients(void* lp_data,
                                             int num_cols,
                                             int num_rows,
                                             SleqpSparseMatrix* coeff_matrix)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  assert(sleqp_sparse_matrix_get_num_rows(coeff_matrix) == num_rows);
  assert(sleqp_sparse_matrix_get_num_cols(coeff_matrix) == num_cols);

  const int* coeff_matrix_cols = sleqp_sparse_matrix_get_cols(coeff_matrix);
  const int* coeff_matrix_rows = sleqp_sparse_matrix_get_rows(coeff_matrix);
  double* coeff_matrix_data = sleqp_sparse_matrix_get_data(coeff_matrix);

  // Note: We save / restore the basis in order to
  //       warm-start the iteration.
  soplex.getBasis(spx->basis_rows, spx->basis_cols);

  soplex.clearBasis();
  assert(soplex.status() == soplex::SPxSolver::NO_PROBLEM);

  for(int j = 0; j < num_cols; ++j)
  {
    int num_entries = coeff_matrix_cols[j + 1] - coeff_matrix_cols[j];

    int offset = coeff_matrix_cols[j];

    soplex::DSVectorReal soplex_col(num_entries);

    soplex_col.add(num_entries,
                   coeff_matrix_rows + offset,
                   coeff_matrix_data + offset);

    double objective = soplex.objReal(j);

    double lb = soplex.lowerReal(j);
    double ub = soplex.upperReal(j);

    soplex.changeColReal(j, soplex::LPCol(objective, soplex_col, ub, lb));
  }

  soplex.setBasis(spx->basis_rows, spx->basis_cols);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_set_objective(void* lp_data,
                                          int num_cols,
                                          int num_rows,
                                          double* objective)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int j = 0; j < num_cols; ++j)
  {
    soplex.changeObjReal(j, adjust_inf(objective[j]));
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE soplex_get_primal_sol(void* lp_data,
                                           int num_cols,
                                           int num_rows,
                                           double* objective_value,
                                           double* solution_values)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  if(objective_value)
  {
    *objective_value = soplex.objValueReal();
  }

  bool found_solution = soplex.getPrimalReal(solution_values,
                                             num_cols);

  assert(found_solution);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_get_dual_sol(void* lp_data,
                                         int num_cols,
                                         int num_rows,
                                         double* vars_dual,
                                         double* cons_dual)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  if(vars_dual)
  {
    bool found_solution = soplex.getRedCostReal(vars_dual,
                                                num_cols);

    assert(found_solution);
  }

  if(cons_dual)
  {
    bool found_solution = soplex.getDualReal(cons_dual,
                                             num_rows);

    assert(found_solution);
  }

  return SLEQP_OKAY;
}

static SLEQP_BASESTAT basestat_for(soplex::SPxSolver::VarStatus status)
{
  switch (status)
  {
  case soplex::SPxSolver::ON_LOWER:
    return SLEQP_BASESTAT_LOWER;
  case soplex::SPxSolver::ON_UPPER:
    return SLEQP_BASESTAT_UPPER;
  case soplex::SPxSolver::ZERO:
    return SLEQP_BASESTAT_ZERO;
  case soplex::SPxSolver::FIXED:
    return SLEQP_BASESTAT_UPPER;
  case soplex::SPxSolver::BASIC:
    return SLEQP_BASESTAT_BASIC;
  default:
    break;
  }

  // Invalid basis status reported by Soplex
  assert(false);

  return SLEQP_BASESTAT_ZERO;
}

static SLEQP_RETCODE soplex_get_varstats(void* lp_data,
                                         int num_cols,
                                         int num_rows,
                                         SLEQP_BASESTAT* variable_stats)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int i = 0; i < num_cols; ++i)
  {
    variable_stats[i] = basestat_for(soplex.basisColStatus(i));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_get_consstats(void* lp_data,
                                          int num_cols,
                                          int num_rows,
                                          SLEQP_BASESTAT* constraint_stats)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int i = 0; i < num_rows; ++i)
  {
    constraint_stats[i] = basestat_for(soplex.basisRowStatus(i));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_get_basis_condition(void* lp_data,
                                                bool* exact,
                                                double* condition)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;

  soplex::SoPlex& soplex = *(spx->soplex);

  if(*exact)
  {
    bool success = soplex.getExactCondition(*condition);

    if(!success)
    {
      sleqp_log_error("Failed to get basis condition");
      return SLEQP_INTERNAL_ERROR;
    }
  }
  else
  {
    bool success = soplex.getEstimatedCondition(*condition);

    if(!success)
    {
      sleqp_log_error("Failed to get basis condition");
      return SLEQP_INTERNAL_ERROR;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_free(void** lp_data)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) *lp_data;

  delete[] spx->basis_cols;
  delete[] spx->basis_rows;
  spx->basis_cols = NULL;
  spx->basis_rows = NULL;

  delete spx->soplex;
  spx->soplex = NULL;

  delete spx;

  *lp_data = NULL;

  return SLEQP_OKAY;
}

extern "C"
{
  SLEQP_RETCODE sleqp_lpi_soplex_create_interface(SleqpLPi** lp_star,
                                                  int num_cols,
                                                  int num_rows,
                                                  SleqpParams* params)
  {
    SleqpLPiCallbacks callbacks = {
      .create_problem = soplex_create_problem,
      .solve = soplex_solve,
      .set_bounds = soplex_set_bounds,
      .set_coefficients = soplex_set_coefficients,
      .set_objective = soplex_set_objective,
      .get_primal_sol = soplex_get_primal_sol,
      .get_dual_sol = soplex_get_dual_sol,
      .get_varstats = soplex_get_varstats,
      .get_consstats = soplex_get_consstats,
      .get_basis_condition = soplex_get_basis_condition,
      .free_problem = soplex_free
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
    SLEQP_CALL(sleqp_lpi_soplex_create_interface(lp_interface,
                                                 num_variables,
                                                 num_constraints,
                                                 params));

    return SLEQP_OKAY;
  }

}
