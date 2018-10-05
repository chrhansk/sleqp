#include "sleqp_lpi_soplex.h"

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
                                  double* objective,
                                  SleqpSparseMatrix* cons_matrix,
                                  double* cons_lb,
                                  double* cons_ub,
                                  double* vars_lb,
                                  double* vars_ub)
{
  // TODO: what about inf values?!

  SleqpLpiSoplex* spx = (SleqpLpiSoplex*) lp_data;
  soplex::SoPlex& soplex = *(spx->soplex);

  assert(soplex.numRowsReal() == (int) spx->num_constraints);
  assert(soplex.numColsReal() == (int) spx->num_variables);


  assert(cons_matrix->num_cols == spx->num_variables);
  assert(cons_matrix->num_rows == spx->num_constraints);

  soplex.changeRangeReal(soplex::VectorReal(spx->num_constraints, cons_lb),
                         soplex::VectorReal(spx->num_constraints, cons_ub));

  int j = 0;

  for(int k = 0; k < cons_matrix->nnz; ++k)
  {
    while(cons_matrix->cols[j] < k)
    {
      ++j;
    }

    assert(cons_matrix->cols[j + 1] >= cons_matrix->cols[j]);

    int num_entries = cons_matrix->cols[j + 1] - cons_matrix->cols[j];

    soplex::Nonzero<soplex::Real>* entries = new soplex::Nonzero<soplex::Real>[num_entries];

    soplex::SVectorReal column(num_entries, entries);

    /* This does not work... int* vs int* for second arg.
     * TODO: It DOES work, reimplement!
    column.add(num_entries,
               cons_matrix->rows + k,
               cons_matrix->data + k);
    */

    for(int i = 0; i < num_entries; ++i)
    {
      column.add(cons_matrix->rows[k + i], cons_matrix->data[k + i]);
    }

    soplex.changeColReal(j, soplex::LPCol(objective[j],
                                          column,
                                          vars_ub[j],
                                          vars_lb[j]));

    delete[] entries;
  }

  soplex.changeObjReal(soplex::VectorReal(spx->num_variables,
                                          objective));

  soplex.changeBoundsReal(soplex::VectorReal(spx->num_variables, vars_lb),
                          soplex::VectorReal(spx->num_variables, vars_ub));

  //soplex.writeFileReal("test.lp");

  soplex::SPxSolver::Status status = soplex.optimize();

  assert(status == soplex::SPxSolver::OPTIMAL);

  return SLEQP_OKAY;
}


static SLEQP_RETCODE soplex_get_solution(void* lp_data,
                                         int num_variables,
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
