#include "eqp.h"

#include "mem.h"

struct SleqpEQPSolver
{
  int refcount;
  SleqpTimer* timer;
  SleqpEQPCallbacks callbacks;
  void* eqp_data;
};

SLEQP_RETCODE sleqp_eqp_solver_create(SleqpEQPSolver** star,
                                      SleqpEQPCallbacks* callbacks,
                                      void* eqp_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpEQPSolver* solver = *star;

  *solver = (SleqpEQPSolver) {0};

  solver->refcount = 1;

  SLEQP_CALL(sleqp_timer_create(&solver->timer));

  solver->callbacks = *callbacks;
  solver->eqp_data = eqp_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_eqp_solver_set_iterate(SleqpEQPSolver* solver,
                                           SleqpIterate* iterate,
                                           SleqpAugJac* jacobian,
                                           double trust_radius,
                                           double penalty_parameter)
{
  SLEQP_CALL(solver->callbacks.set_iterate(iterate,
                                           jacobian,
                                           trust_radius,
                                           penalty_parameter,
                                           solver->eqp_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_eqp_solver_set_time_limit(SleqpEQPSolver* solver,
                                              double time_limit)
{
  SLEQP_CALL(solver->callbacks.set_time_limit(time_limit,
                                              solver->eqp_data));

  return SLEQP_OKAY;
}

SleqpTimer* sleqp_eqp_solver_get_timer(SleqpEQPSolver* solver)
{
  return solver->timer;
}

SLEQP_RETCODE sleqp_eqp_solver_add_violated_multipliers(SleqpEQPSolver* solver,
                                                        SleqpSparseVec* multipliers)
{
  SLEQP_CALL(solver->callbacks.add_violated_multipliers(multipliers,
                                                        solver->eqp_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_eqp_solver_compute_step(SleqpEQPSolver* solver,
                                            const SleqpSparseVec* multipliers,
                                            SleqpSparseVec* newton_step)
{
  SLEQP_CALL(sleqp_timer_start(solver->timer));

  SLEQP_RETCODE status = solver->callbacks.compute_step(multipliers,
                                                        newton_step,
                                                        solver->eqp_data);

  SLEQP_CALL(sleqp_timer_stop(solver->timer));

  return status;
}

SLEQP_RETCODE sleqp_eqp_solver_current_rayleigh(SleqpEQPSolver* solver,
                                                double* min_rayleigh,
                                                double* max_rayleigh)
{
  SLEQP_CALL(solver->callbacks.current_rayleigh(min_rayleigh,
                                                max_rayleigh,
                                                solver->eqp_data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
eqp_solver_free(SleqpEQPSolver** star)
{
  SleqpEQPSolver* solver = *star;

  SLEQP_CALL(solver->callbacks.free(solver->eqp_data));

  SLEQP_CALL(sleqp_timer_free(&solver->timer));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_eqp_solver_capture(SleqpEQPSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_eqp_solver_release(SleqpEQPSolver** star)
{
  SleqpEQPSolver* eqp_solver = *star;

  if(!eqp_solver)
  {
    return SLEQP_OKAY;
  }

  if(--eqp_solver->refcount == 0)
  {
    SLEQP_CALL(eqp_solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
