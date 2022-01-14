#ifndef SLEQ_AMPL_OUTPUT_H
#define SLEQ_AMPL_OUTPUT_H

#include <asl.h>

#include "sleqp.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_ampl_report(SleqpProblem* problem,
                  SleqpSolver* solver,
                  ASL* asl,
                  Option_Info* option_info,
                  bool error_occurred);

#endif /* SLEQ_AMPL_OUTPUT_H */
