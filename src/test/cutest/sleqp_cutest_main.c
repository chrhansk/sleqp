#include <stdlib.h>
#include <getopt.h>

#include "pub_log.h"
#include "pub_types.h"

#include "sleqp_cutest_driver.h"
#include "sleqp_cutest_options.h"

int
parse_command_line_options(int argc, char *argv[], SleqpCutestOptions* options)
{
  while(true)
  {
    int option_index = 0;

    static struct option long_options[] = {
      {"enable_logging",              no_argument,       0, 'l'},
      {"enable_preprocessing",        no_argument,       0, 'p'},
      {"force_nonlinear_constraints", no_argument,       0, 'n'},
      {"max_num_threads",             required_argument, 0, 't'},
      {0,                             0,                 0,  0}
    };

    int c = getopt_long(argc,
                        argv,
                        "pn",
                        long_options,
                        &option_index);
    if (c == -1)
      break;

    switch(c) {
    case 'l':
      break;
      sleqp_log_debug("Enabling logging");
      options->enable_logging = true;
      break;
    case 'p':
      sleqp_log_debug("Enabling preprocessing");
      options->enable_preprocessing = true;
      break;

    case 'n':
      sleqp_log_debug("Forcing nonlinear constraints");
      options->force_nonlinear_constraints = true;
      break;

    case 't':
      options->max_num_threads = atoi(optarg);
      sleqp_log_debug("Using up to %d threads", options->max_num_threads);
      break;

    default:
      sleqp_log_error("Invalid option %o", c);
      return EXIT_FAILURE;
    }
  }

  if(optind < argc)
  {
    sleqp_log_error("Unexpected positional arguments");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  const char* filename = PROBLEM_OUTSDIF_D; /* CUTEst data file */
  const char* probname = PROBLEM_NAME;

  SleqpCutestOptions options = (SleqpCutestOptions) {0};
  options.max_num_threads = SLEQP_NONE;

  if(parse_command_line_options(argc, argv, &options) != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  return sleqp_cutest_run(filename, probname, &options);
}
