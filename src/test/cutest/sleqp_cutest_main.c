#include "sleqp.h"

#include "sleqp_cutest_driver.h"

int main(int argc, char *argv[])
{
  const char* filename = PROBLEM_OUTSDIF_D; /* CUTEst data file */
  const char* probname = PROBLEM_NAME;

  return sleqp_cutest_run(filename, probname);
}
