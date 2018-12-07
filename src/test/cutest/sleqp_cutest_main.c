#include "sleqp.h"

#include "sleqp_cutest_driver.h"

int main(int argc, char *argv[])
{
  const char *fname = PROBLEM_OUTSDIF_D; /* CUTEst data file */

  return sleqp_cutest_run(fname);
}
