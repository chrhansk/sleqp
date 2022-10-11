#ifndef DYN_CONSTRAINED_FIXTURE_H
#define DYN_CONSTRAINED_FIXTURE_H

#include "func.h"

#include "constrained_fixture.h"

extern SleqpFunc* dyn_constrained_func;

void
dyn_constrained_setup();

void
dyn_constrained_teardown();

#endif /* DYN_CONSTRAINED_FIXTURE_H */
