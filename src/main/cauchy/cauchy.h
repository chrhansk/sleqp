#ifndef SLEQP_CAUCHY_H
#define SLEQP_CAUCHY_H

#include "cauchy_types.h"

typedef struct SleqpCauchy SleqpCauchy;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_create(SleqpCauchy** cauchy,
                    SleqpCauchyCallbacks* callbacks,
                    void* cauchy_data);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_set_iterate(SleqpCauchy* cauchy,
                         SleqpIterate* iterate,
                         double trust_radius);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_set_trust_radius(SleqpCauchy* cauchy, double trust_radius);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_solve(SleqpCauchy* cauchy,
                   SleqpVec* gradient,
                   double penalty,
                   SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_obj_val(SleqpCauchy* cauchy, double* objective_value);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_working_set(SleqpCauchy* cauchy, SleqpIterate* iterate);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_lp_step(SleqpCauchy* cauchy, SleqpVec* direction);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_locally_infeasible(SleqpCauchy* cauchy, bool* locally_infeasible);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_estimate_duals(SleqpCauchy* cauchy,
                            const SleqpWorkingSet* working_set,
                            SleqpVec* cons_dual,
                            SleqpVec* vars_dual);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_violation(SleqpCauchy* cauchy, double* violation);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_set_time_limit(SleqpCauchy* cauchy, double time_limit);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_basis_condition(SleqpCauchy* cauchy,
                             bool* exact,
                             double* condition);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_print_stats(SleqpCauchy* cauchy, double total_elapsed);

// Bound on the criticality measure used in
// "On the Convergence of Successive Linear Programming Algorithms"
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_compute_criticality_bound(SleqpCauchy* cauchy,
                                       double merit_value,
                                       double* criticality_bound);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_capture(SleqpCauchy* cauchy);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cauchy_release(SleqpCauchy** star);

#endif /* SLEQP_CAUCHY_H */
