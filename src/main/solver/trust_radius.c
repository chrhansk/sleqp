#include "solver.h"

#include "cmp.h"

SLEQP_RETCODE sleqp_solver_update_lp_trust_radius(SleqpSolver* solver,
                                                  bool trial_step_accepted,
                                                  double trial_step_infnorm,
                                                  double cauchy_step_infnorm,
                                                  double cauchy_step_length,
                                                  double eps,
                                                  double* lp_trust_radius)
{
  if(trial_step_accepted)
  {
    double norm_increase_factor = 1.2;

    trial_step_infnorm *= norm_increase_factor;
    cauchy_step_infnorm *= norm_increase_factor;

    double scaled_trust_radius = .1 * (*lp_trust_radius);

    double update_lhs = SLEQP_MAX(trial_step_infnorm,
                                  cauchy_step_infnorm);

    update_lhs = SLEQP_MAX(update_lhs, scaled_trust_radius);

    if(sleqp_is_eq(cauchy_step_length, 1., eps))
    {
      (*lp_trust_radius) *= 7.;
    }

    *lp_trust_radius = SLEQP_MIN(update_lhs, *lp_trust_radius);

  }
  else
  {
    double half_norm = .5 * trial_step_infnorm;
    double small_radius = .1 * (*lp_trust_radius);

    double reduced_radius = SLEQP_MAX(half_norm, small_radius);

    *lp_trust_radius = SLEQP_MIN(reduced_radius, *lp_trust_radius);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_update_trust_radius(SleqpSolver* solver,
                                               double reduction_ratio,
                                               bool trial_step_accepted,
                                               double direction_norm)
{
  const double eps = sleqp_params_get(solver->params, SLEQP_PARAM_EPS);

  double* trust_radius = &(solver->trust_radius);

  if(reduction_ratio >= 0.9)
  {
    *trust_radius = SLEQP_MAX(*trust_radius,
                              7*direction_norm);
  }
  else if(reduction_ratio >= 0.3)
  {
    *trust_radius = SLEQP_MAX(*trust_radius,
                              2*direction_norm);
  }
  else if(trial_step_accepted)
  {
    // stays the same
  }
  else
  {
    // filter out very small steps
    if(sleqp_is_zero(direction_norm, eps))
    {
      *trust_radius *= .5;
    }
    else
    {
      *trust_radius = SLEQP_MIN(.5 * (*trust_radius),
                                .5 * direction_norm);
    }
  }

  return SLEQP_OKAY;
}
