#ifndef SLEQP_TYPES_H
#define SLEQP_TYPES_H

#include "pub_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SLEQP_SOLVER_PHASE_OPTIMIZATION = 0,
    SLEQP_SOLVER_PHASE_RESTORATION,
    SLEQP_SOLVER_NUM_PHASES
  } SLEQP_SOLVER_PHASE;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TYPES_H */
