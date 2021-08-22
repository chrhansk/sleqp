#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "cutest_defs.h"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "pub_log.h"
#include "pub_types.h"

#include "sleqp_cutest_driver.h"
#include "sleqp_cutest_options.h"

const int minute = 60;
const int grace_minutes = 1;

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
      {"time_limit"     ,             required_argument, 0, 's'},
      {"output",                      required_argument, 0, 'o'},
      {0,                             0,                 0,  0}
    };

    int c = getopt_long(argc,
                        argv,
                        "lpntso",
                        long_options,
                        &option_index);
    if (c == -1)
      break;

    switch(c) {
    case 'l':
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

#ifdef WITH_OPENMP
      omp_set_num_threads(options->max_num_threads);
#endif
      break;

    case 's':
      options->time_limit = atof(optarg);
      sleqp_log_debug("Setting time limit to %s", optarg);
      break;

    case 'o':
      options->output = optarg;
      sleqp_log_debug("Setting output to %s", optarg);
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

int handle_cutest_child(const char* probname,
                        double time_limit,
                        pid_t child_pid)
{
  sigset_t child_mask, old_mask;
  sigemptyset(&child_mask);
  sigaddset(&child_mask, SIGCHLD);

  if (sigprocmask(SIG_BLOCK, &child_mask, &old_mask) == -1)
  {
    sleqp_log_error("Failed to set signal mask: %s",
                    strerror(errno));

    return EXIT_FAILURE;
  }

  struct timespec ts = {0};
  ts.tv_sec = minute;

  bool exited = false;

  const int num_minutes = (time_limit / (double) minute) + grace_minutes;

  for(int i = 0; i < num_minutes; ++i)
  {
    int ret = sigtimedwait(&child_mask, NULL, &ts);

    if(ret > 0)
    {
      exited = true;
      break;
    }
    else
    {
      assert(ret == -1);

      if(errno == EAGAIN)
      {
        continue;
      }

      sleqp_log_error("Failed to wait for child termination: %s",
                      strerror(errno));

      break;
    }
  }

  if(!exited)
  {
    sleqp_log_error("Killing solver %s after deadline expired",
                    probname);

    kill(child_pid, SIGTERM);
  }

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  const char* filename = PROBLEM_OUTSDIF_D; /* CUTEst data file */
  const char* probname = PROBLEM_NAME;

  SleqpCutestOptions options = (SleqpCutestOptions) {0};

  if(parse_command_line_options(argc, argv, &options) != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  pid_t pid = fork();

  if(pid == -1)
  {
    sleqp_log_error("Failed to fork(): %s", strerror(errno));
    return EXIT_FAILURE;
  }
  else if(pid == 0)
  {
    // child
    return sleqp_cutest_run(filename, probname, &options);
  }
  else
  {
    // parent
    assert(pid > 0);

    return handle_cutest_child(probname,
                               options.time_limit,
                               pid);
  }

  return EXIT_SUCCESS;
}
