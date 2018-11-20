#include "sleqp_lpi_soplex.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include <soplex.h>

struct SleqpLpiSoplex
{
  soplex::SoPlex* soplex;
  int num_variables;
  int num_constraints;
};

static SLEQP_RETCODE soplex_create_problem(void** lp_data,
                                           int num_variables,
                                           int num_constraints)
{
  SleqpLpiSoplex* spx = new SleqpLpiSoplex;

  spx->soplex = new soplex::SoPlex();
  soplex::SoPlex& soplex = *(spx->soplex);

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

  spx->num_variables = num_variables;
  spx->num_constraints = num_constraints;

  soplex.setIntParam(soplex::SoPlex::OBJSENSE,
                     soplex::SoPlex::OBJSENSE_MINIMIZE);

  // add dummy (empty) rows / cols
  soplex::DSVector vec(0);

  double inf = soplex::infinity;

  for(int j = 0; j < num_variables; ++j)
  {
    soplex.addColReal(soplex::LPCol(0., vec, inf, -inf));
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    soplex.addRowReal(soplex::LPRow(-inf, vec, inf));
  }

  *lp_data = spx;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_solve(void* lp_data,
                                  int num_variables,
                                  int num_constraints)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  soplex::SPxSolver::Status status = soplex.optimize();

  //soplex.writeFileReal("test.lp");

  assert(status == soplex::SPxSolver::OPTIMAL);

  assert(soplex.hasBasis());

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
                                       int num_variables,
                                       int num_constraints,
                                       double* cons_lb,
                                       double* cons_ub,
                                       double* vars_lb,
                                       double* vars_ub)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int i = 0; i < num_constraints; ++i)
  {
    soplex.changeRangeReal(i, adjust_inf(cons_lb[i]), adjust_inf(cons_ub[i]));
  }

  for(int j = 0; j < num_variables; ++j)
  {
    soplex.changeBoundsReal(j, adjust_inf(vars_lb[j]), adjust_inf(vars_ub[j]));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_set_coefficients(void* lp_data,
                                             int num_variables,
                                             int num_constraints,
                                             SleqpSparseMatrix* coeff_matrix)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  assert(num_variables == coeff_matrix->num_cols);
  assert(num_constraints == coeff_matrix->num_rows);

  // TODO: Find out whether storing / restoring the basis
  //       is a good idea.
  soplex.clearBasis();
  assert(soplex.status() == soplex::SPxSolver::NO_PROBLEM);

  for(int j = 0; j < num_variables; ++j)
  {
    int num_entries = coeff_matrix->cols[j + 1] - coeff_matrix->cols[j];

    int offset = coeff_matrix->cols[j];

    soplex::DSVectorReal soplex_col(num_entries);

    soplex_col.add(num_entries,
                   coeff_matrix->rows + offset,
                   coeff_matrix->data + offset);

    double objective = soplex.objReal(j);

    double lb = soplex.lowerReal(j);
    double ub = soplex.upperReal(j);

    soplex.changeColReal(j, soplex::LPCol(objective, soplex_col, ub, lb));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_set_objective(void* lp_data,
                                          int num_variables,
                                          int num_constraints,
                                          double* objective)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int j = 0; j < num_variables; ++j)
  {
    soplex.changeObjReal(j, adjust_inf(objective[j]));
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE soplex_get_solution(void* lp_data,
                                         int num_variables,
                                         int num_constraints,
                                         double* objective_value,
                                         double* solution_values)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  if(objective_value)
  {
    *objective_value = soplex.objValueReal();
  }

  soplex::VectorReal solution(num_variables, solution_values);

  bool found_solution = soplex.getPrimalReal(solution);

  assert(found_solution);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_get_varstats(void* lp_data,
                                         int num_variables,
                                         int num_constraints,
                                         SLEQP_BASESTAT* variable_stats)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int i = 0; i < num_variables; ++i)
  {
    switch (soplex.basisColStatus(i))
    {
    case soplex::SPxSolver::ON_LOWER:
      variable_stats[i] = SLEQP_BASESTAT_LOWER;
      break;
    case soplex::SPxSolver::ON_UPPER:
      variable_stats[i] = SLEQP_BASESTAT_UPPER;
      break;
    case soplex::SPxSolver::ZERO:
      variable_stats[i] = SLEQP_BASESTAT_ZERO;
      break;
    case soplex::SPxSolver::FIXED:
    case soplex::SPxSolver::BASIC:
      variable_stats[i] = SLEQP_BASESTAT_BASIC;
      break;
    default:
      break;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_get_consstats(void* lp_data,
                                          int num_variables,
                                          int num_constraints,
                                          SLEQP_BASESTAT* constraint_stats)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(int i = 0; i < num_constraints; ++i)
  {
    switch (soplex.basisRowStatus(i))
    {
    case soplex::SPxSolver::ON_LOWER:
      constraint_stats[i] = SLEQP_BASESTAT_LOWER;
      break;
    case soplex::SPxSolver::ON_UPPER:
      constraint_stats[i] = SLEQP_BASESTAT_UPPER;
      break;
    case soplex::SPxSolver::ZERO:
      constraint_stats[i] = SLEQP_BASESTAT_ZERO;
      break;
    case soplex::SPxSolver::FIXED:
    case soplex::SPxSolver::BASIC:
      constraint_stats[i] = SLEQP_BASESTAT_BASIC;
      break;
    default:
      break;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_free(void** lp_data)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) *lp_data;

  delete spx->soplex;

  delete spx;

  *lp_data = NULL;

  return SLEQP_OKAY;
}

extern "C"
{
  SLEQP_RETCODE sleqp_lpi_soplex_create_interface(SleqpLPi** lp_star,
                                                  int num_variables,
                                                  int num_constraints)
  {
    return sleqp_lpi_create_interface(lp_star,
                                      num_variables,
                                      num_constraints,
                                      soplex_create_problem,
                                      soplex_solve,
                                      soplex_set_bounds,
                                      soplex_set_coefficients,
                                      soplex_set_objective,
                                      soplex_get_solution,
                                      soplex_get_varstats,
                                      soplex_get_consstats,
                                      soplex_free);
  }
}
