#ifndef DYN_ROSENBROCK_FIXTURE_H
#define DYN_ROSENBROCK_FIXTURE_H

#include "rosenbrock_fixture.h"

#ifdef __cplusplus
extern "C" {
#endif

  extern SleqpFunc* dyn_rosenbrock_func;

  void dyn_rosenbrock_setup();

  void dyn_rosenbrock_teardown();

#ifdef __cplusplus
}
#endif

#endif /* DYN_ROSENBROCK_FIXTURE_H */
