#include "penalty.h"

#include "log.h"

const double penalty_increase    = 10.;
const double violation_tolerance = 1e-8;
const double min_decrease        = .1;
const int max_increases          = 100;

SLEQP_RETCODE
sleqp_update_penalty(SleqpProblem* problem,
                     SleqpIterate* iterate,
                     SleqpCauchy* cauchy_data,
                     double* penalty_parameter,
                     bool* locally_infeasible)
{
  const int num_constraints = sleqp_problem_num_cons(problem);

  (*locally_infeasible) = false;

  if (num_constraints == 0)
  {
    return SLEQP_OKAY;
  }

  double current_violation;

  SLEQP_CALL(sleqp_cauchy_get_violation(cauchy_data, &current_violation));

  current_violation /= num_constraints;

  sleqp_log_debug("Updating penalty parameter, average violation is %.10e",
                  current_violation);

  if (current_violation <= violation_tolerance)
  {
    sleqp_log_debug("Average violation is already below the tolerance of %.10e",
                    violation_tolerance);

    return SLEQP_OKAY;
  }

  sleqp_log_debug(
    "Resolving linearization to compute minimum average violation");

  SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                NULL,
                                (*penalty_parameter),
                                SLEQP_CAUCHY_OBJECTIVE_TYPE_FEASIBILITY));

  {
    SLEQP_CALL(
      sleqp_cauchy_locally_infeasible(cauchy_data, locally_infeasible));

    if (*locally_infeasible)
    {
      sleqp_log_warn("Current iterate is locally infeasible");
    }
  }

  double inf_violation;

  SLEQP_CALL(sleqp_cauchy_get_violation(cauchy_data, &inf_violation));

  inf_violation /= num_constraints;

  sleqp_log_debug("Minimum average violation: %.10e", inf_violation);

  // sleqp_assert_is_geq(current_violation, inf_violation, eps);

  if (inf_violation <= violation_tolerance)
  {
    sleqp_log_debug("Minimum average violation is below tolerance");

    for (int i = 0; i < max_increases; ++i)
    {
      (*penalty_parameter) *= penalty_increase;

      sleqp_log_debug("Resolving linearization to compute average violation "
                      "for penalty value %e",
                      (*penalty_parameter));

      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    sleqp_iterate_obj_grad(iterate),
                                    (*penalty_parameter),
                                    SLEQP_CAUCHY_OBJECTIVE_TYPE_MIXED));

      double next_violation;

      SLEQP_CALL(sleqp_cauchy_get_violation(cauchy_data, &next_violation));

      next_violation /= num_constraints;

      sleqp_log_debug("Average violation for penalty value %e is %.10e",
                      (*penalty_parameter),
                      next_violation);

      if (next_violation <= violation_tolerance)
      {
        sleqp_log_debug("Average violation is below the tolerance of %e",
                        (*penalty_parameter));

        return SLEQP_OKAY;
      }
      else
      {
        sleqp_log_debug(
          "Average violation is above the tolerance of %e, continuing",
          (*penalty_parameter));
      }
    }
  }
  else
  {
    sleqp_log_debug("Minimum average violation is above tolerance");

    if (current_violation - inf_violation <= violation_tolerance)
    {
      sleqp_log_debug("Cannot make progress towards feasibility, aborting");
      // we can't make progress in feasibility, no need for an increase
      return SLEQP_OKAY;
    }

    for (int i = 0; i < max_increases; ++i)
    {
      (*penalty_parameter) *= penalty_increase;

      sleqp_log_debug("Resolving linearization to compute average violation "
                      "for penalty value %e",
                      (*penalty_parameter));

      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    sleqp_iterate_obj_grad(iterate),
                                    (*penalty_parameter),
                                    SLEQP_CAUCHY_OBJECTIVE_TYPE_MIXED));

      double next_violation;

      SLEQP_CALL(sleqp_cauchy_get_violation(cauchy_data, &next_violation));

      next_violation /= num_constraints;

      sleqp_log_debug("Average violation for penalty value %e is %.10e",
                      (*penalty_parameter),
                      next_violation);

      if ((current_violation - next_violation)
          >= min_decrease * (current_violation - inf_violation))
      {
        sleqp_log_debug("Penalty value of %e achieves sufficiently high "
                        "reduction in average violation",
                        (*penalty_parameter));

        return SLEQP_OKAY;
      }
      else
      {
        sleqp_log_debug("Penalty value of %e does not achieve sufficiently "
                        "high reduction in average violation",
                        (*penalty_parameter));
      }
    }
  }

  return SLEQP_OKAY;
}
