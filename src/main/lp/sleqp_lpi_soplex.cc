#include "sleqp_lpi_soplex.h"

#include "sleqp_mem.h"

#include <soplex.h>

struct SleqpLpiSoplex
{
  soplex::SoPlex* soplex;
  size_t num_variables;
  size_t num_constraints;
};

static SLEQP_RETCODE soplex_create_problem(void** lp_data,
                                           size_t num_variables,
                                           size_t num_constraints)
{
  SleqpLpiSoplex* spx = new SleqpLpiSoplex;

  spx->soplex = new soplex::SoPlex();
  soplex::SoPlex& soplex = *(spx->soplex);

  spx->num_variables = num_variables;
  spx->num_constraints = num_constraints;

  soplex.setIntParam(soplex::SoPlex::OBJSENSE,
                     soplex::SoPlex::OBJSENSE_MINIMIZE);

  // add dummy (empty) rows / cols
  soplex::DSVector vec(0);

  double inf = soplex::infinity;

  for(size_t j = 0; j < num_variables; ++j)
  {
    soplex.addColReal(soplex::LPCol(0., vec, inf, -inf));
  }

  for(size_t i = 0; i < num_constraints; ++i)
  {
    soplex.addRowReal(soplex::LPRow(-inf, vec, inf));
  }

  *lp_data = spx;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_solve(void* lp_data,
                                  double* objective,
                                  SleqpSparseMatrix* cons_matrix,
                                  double* cons_lb,
                                  double* cons_ub,
                                  double* vars_lb,
                                  double* vars_ub)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  //soplex.changeObjReal(soplex::VectorReal(spx->num_variables, objective));

  soplex.changeRangeReal(soplex::VectorReal(spx->num_variables, cons_lb),
                         soplex::VectorReal(spx->num_variables, cons_ub));

  /*
  soplex.changeLowerReal(soplex::VectorReal(spx->num_constraints, cons_lb));
  soplex.changeUpperReal(soplex::VectorReal(spx->num_constraints, cons_ub));
  */


  assert(cons_matrix->num_cols == spx->num_variables);
  assert(cons_matrix->num_rows == spx->num_constraints);

  size_t j = 0;

  for(size_t k = 0; k < cons_matrix->nnz; ++k)
  {
    while(cons_matrix->cols[j] < k)
    {
      ++j;
    }

    assert(cons_matrix->cols[j + 1] >= cons_matrix->cols[j]);

    size_t num_entries = cons_matrix->cols[j + 1] - cons_matrix->cols[j];

    soplex::Nonzero<soplex::Real>* entries = new soplex::Nonzero<soplex::Real>[num_entries];

    soplex::SVectorReal column(num_entries, entries);

    /* This does not work... size_t* vs int* for second arg.
    column.add(num_entries,
               cons_matrix->rows + k,
               cons_matrix->data + k);
    */

    for(size_t i = 0; i < num_entries; ++i)
    {
      column.add(cons_matrix->rows[k + i], cons_matrix->data[k + i]);
    }

    soplex.changeColReal(j, soplex::LPCol(objective[j],
                                          column,
                                          vars_ub[j],
                                          vars_lb[j]));

    delete[] entries;
  }

  soplex::SPxSolver::Status status = soplex.optimize();

  assert(status == soplex::SPxSolver::OPTIMAL);

  return SLEQP_OKAY;
}


static SLEQP_RETCODE soplex_get_solution(void* lp_data,
                                         size_t num_variables,
                                         double* objective_value,
                                         double* solution_values)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  *objective_value = soplex.objValueReal();

  soplex::VectorReal solution(num_variables, solution_values);

  soplex.getPrimalReal(solution);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE soplex_get_varstats(void* lp_data,
                                         size_t num_variables,
                                         SLEQP_BASESTAT* variable_stats)
{
  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  for(size_t i = 0; i < num_variables; ++i)
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
  SLEQP_RETCODE sleqp_lpi_soplex_create_interface(SleqpLPi** lp_star)
  {
    return sleqp_lpi_create_interface(lp_star,
                                      soplex_create_problem,
                                      soplex_solve,
                                      soplex_get_solution,
                                      soplex_get_varstats,
                                      soplex_free);
  }
}
