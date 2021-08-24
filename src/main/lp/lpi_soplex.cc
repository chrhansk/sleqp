#include "lpi_soplex.h"

#include <iostream>
#include <vector>

#include <soplex.h>

#include "cmp.h"
#include "defs.h"
#include "log.h"
#include "mem.h"

static const double tolerance_factor = 1e-1;

struct SoPlexBasis
{
  SoPlexBasis(int num_cols, int num_rows)
    : status(num_cols + num_rows, soplex::SPxSolver::ZERO),
      basis_rows(status.data()),
      basis_cols(status.data() + num_rows)
  {
  }

  std::vector<soplex::SPxSolver::VarStatus> status;

  soplex::SPxSolver::VarStatus* const basis_rows;
  soplex::SPxSolver::VarStatus* const basis_cols;
};

struct SleqpLpiSoplex
{
  SleqpLpiSoplex(int num_cols, int num_rows)
    : current_basis(num_cols, num_rows),
      num_cols(num_cols),
      num_rows(num_rows)
  {}

  soplex::SoPlex* soplex;

  std::vector<SoPlexBasis> bases;

  SoPlexBasis current_basis;

  int num_cols;
  int num_rows;

  SLEQP_LPI_STATUS status;
};

static SLEQP_RETCODE soplex_create_problem(void** lp_data,
                                           int num_cols,
                                           int num_rows,
                                           SleqpParams* params,
                                           SleqpOptions* options)
{
  SleqpLpiSoplex* spx = new SleqpLpiSoplex(num_cols, num_rows);

  spx->soplex = new soplex::SoPlex();
  soplex::SoPlex& soplex = *(spx->soplex);

  assert(soplex.setRealParam(soplex::SoPlexBase<double>::INFTY,
                             sleqp_infinity()));

  const double feas_eps = sleqp_params_get(params,
                                           SLEQP_PARAM_FEASIBILITY_TOL);

  const double stat_eps = sleqp_params_get(params,
                                           SLEQP_PARAM_STATIONARITY_TOL);

  assert(soplex.setRealParam(soplex::SoPlexBase<double>::FEASTOL,
                             feas_eps*tolerance_factor));

  assert(soplex.setRealParam(soplex::SoPlexBase<double>::OPTTOL,
                             stat_eps*tolerance_factor));

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
  spx->status   = SLEQP_LPI_STATUS_UNKNOWN;

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

static SLEQP_RETCODE soplex_write(void* lp_data, const char* filename)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  soplex.writeFileReal(filename);

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

  switch(status)
  {
  case soplex::SPxSolver::ABORT_TIME:
    spx->status = SLEQP_LPI_STATUS_UNKNOWN;
    return SLEQP_ABORT_TIME;
    break;
  case soplex::SPxSolver::OPTIMAL:
    spx->status = SLEQP_LPI_STATUS_OPTIMAL;
    break;
  case soplex::SPxSolver::UNBOUNDED:
    spx->status = SLEQP_LPI_STATUS_UNBOUNDED;
    break;
  case soplex::SPxSolver::INFEASIBLE:
    spx->status = SLEQP_LPI_STATUS_INF;
    break;
  case soplex::SPxSolver::INForUNBD:
    spx->status = SLEQP_LPI_STATUS_INF_OR_UNBOUNDED;
    break;
  default:
    sleqp_log_error("Invalid SoPlex status: %d", status);
    spx->status = SLEQP_LPI_STATUS_UNKNOWN;
    return SLEQP_INTERNAL_ERROR;
  }

  assert(soplex.hasBasis());

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

  return SLEQP_OKAY;
}

static SLEQP_LPI_STATUS soplex_get_status(void* lp_data)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;

  return spx->status;
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
    soplex.changeRangeReal(i, cons_lb[i], cons_ub[i]);
  }

  for(int j = 0; j < num_cols; ++j)
  {
    assert(vars_lb[j] <= vars_ub[j]);
    soplex.changeBoundsReal(j, vars_lb[j], vars_ub[j]);
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
  soplex.getBasis(spx->current_basis.basis_rows, spx->current_basis.basis_cols);

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

  soplex.setBasis(spx->current_basis.basis_rows, spx->current_basis.basis_cols);

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
    soplex.changeObjReal(j, objective[j]);
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE soplex_save_basis(void* lp_data,
                                       int index)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  assert(index >= 0);

  unsigned int uindex = index;

  while(uindex >= spx->bases.size())
  {
    spx->bases.push_back(SoPlexBasis(spx->num_cols, spx->num_rows));
  }

  SoPlexBasis& basis = spx->bases[index];

  soplex.getBasis(basis.basis_rows, basis.basis_cols);

  return SLEQP_OKAY;
}


static SLEQP_RETCODE soplex_restore_basis(void* lp_data,
                                          int index)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  assert(index >= 0);

  unsigned int uindex = index;

  assert(uindex < spx->bases.size());

  const SoPlexBasis& basis = spx->bases[uindex];

  soplex.setBasis(basis.basis_rows, basis.basis_cols);

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

  if(solution_values)
  {
    const bool found_solution = soplex.getPrimalReal(solution_values,
                                                     num_cols);

    assert(found_solution);
  }

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
                                                  SleqpParams* params,
                                                  SleqpOptions* options)
  {
    SleqpLPiCallbacks callbacks = {
      .create_problem      = soplex_create_problem,
      .solve               = soplex_solve,
      .get_status          = soplex_get_status,
      .set_bounds          = soplex_set_bounds,
      .set_coefficients    = soplex_set_coefficients,
      .set_objective       = soplex_set_objective,
      .save_basis          = soplex_save_basis,
      .restore_basis       = soplex_restore_basis,
      .get_primal_sol      = soplex_get_primal_sol,
      .get_dual_sol        = soplex_get_dual_sol,
      .get_varstats        = soplex_get_varstats,
      .get_consstats       = soplex_get_consstats,
      .get_basis_condition = soplex_get_basis_condition,
      .free_problem        = soplex_free
    };

    return sleqp_lpi_create_interface(lp_star,
                                      SLEQP_LP_SOLVER_SOPLEX_NAME,
                                      SLEQP_LP_SOLVER_SOPLEX_VERSION,
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
    SLEQP_CALL(sleqp_lpi_soplex_create_interface(lp_interface,
                                                 num_variables,
                                                 num_constraints,
                                                 params,
                                                 options));

    return SLEQP_OKAY;
  }

}
