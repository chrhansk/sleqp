#ifndef SLEQP_PENALTY_H
#define SLEQP_PENALTY_H

#include "cauchy/cauchy.h"
#include "problem.h"

SLEQP_RETCODE
sleqp_update_penalty(SleqpProblem* problem,
                     SleqpIterate* iterate,
                     SleqpCauchy* cauchy_data,
                     double* penalty_parameter,
                     bool* locally_infeasible);

#endif /* SLEQP_PENALTY_H */
