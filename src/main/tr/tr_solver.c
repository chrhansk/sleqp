#include "tr_solver.h"

#include "log.h"
#include "mem.h"

struct SleqpTRSolver
{
  int refcount;

  SleqpTRCallbacks callbacks;
  void* solver_data;

  double time_limit;
  SleqpTimer* timer;
};

SLEQP_RETCODE sleqp_tr_solver_set_time_limit(SleqpTRSolver* data,
                                             double time_limit)
{
  data->time_limit = time_limit;

  return SLEQP_OKAY;
}

SleqpTimer* sleqp_tr_solver_get_solve_timer(SleqpTRSolver* data)
{
  return data->timer;
}

SLEQP_RETCODE sleqp_tr_solver_create(SleqpTRSolver** star,
                                     SleqpTRCallbacks* callbacks,
                                     void* solver_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpTRSolver* solver = *star;

  *solver = (SleqpTRSolver) {0};

  solver->refcount = 1;

  solver->callbacks = (*callbacks);
  solver->solver_data = solver_data;

  solver->time_limit = SLEQP_NONE;
  SLEQP_CALL(sleqp_timer_create(&(solver->timer)));

  return SLEQP_OKAY;
}



SLEQP_RETCODE sleqp_tr_solver_solve(SleqpTRSolver* solver,
                                    SleqpAugJac* jacobian,
                                    SleqpSparseVec* multipliers,
                                    SleqpSparseVec* gradient,
                                    SleqpSparseVec* newton_step,
                                    double trust_radius,
                                    double* tr_dual)
{
  SLEQP_CALL(sleqp_timer_start(solver->timer));

  SLEQP_CALL(solver->callbacks.solve(jacobian,
                                     multipliers,
                                     gradient,
                                     newton_step,
                                     trust_radius,
                                     tr_dual,
                                     solver->time_limit,
                                     solver->solver_data));

  SLEQP_CALL(sleqp_timer_stop(solver->timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_tr_solver_current_rayleigh(SleqpTRSolver* solver,
                                               double* min_rayleigh,
                                               double* max_rayleigh)
{
  SLEQP_CALL(solver->callbacks.rayleigh(min_rayleigh,
                                        max_rayleigh,
                                        solver->solver_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_tr_solver_capture(SleqpTRSolver* solver)
{
  ++(solver->refcount);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE tr_solver_free(SleqpTRSolver** star)
{
  SleqpTRSolver* solver = *star;

  SLEQP_CALL(sleqp_timer_free(&(solver->timer)));

  SLEQP_CALL(solver->callbacks.free(&(solver->solver_data)));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_tr_solver_release(SleqpTRSolver** star)
{
  SleqpTRSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  if(--(solver->refcount) == 0)
  {
    SLEQP_CALL(tr_solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
