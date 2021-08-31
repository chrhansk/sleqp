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

#include "log.h"
#include "types.h"
#include "timer.h"

#include "sleqp_cutest_driver.h"
#include "sleqp_cutest_options.h"

const char* filename = PROBLEM_OUTSDIF_D; /* CUTEst data file */
const char* probname = PROBLEM_NAME;

const int grace_seconds = 60;


static int
parse_command_line_options(int argc,
                           char *argv[],
                           SleqpCutestOptions* options)
{
  while(true)
  {
    int option_index = 0;

    static struct option long_options[] = {
      {"enable_logging",              no_argument,       0, 'l'},
      {"enable_preprocessing",        no_argument,       0, 'p'},
      {"force_nonlinear_constraints", no_argument,       0, 'n'},
      {"no_fork",                     no_argument,       0, 'f'},
      {"max_num_threads",             required_argument, 0, 't'},
      {"time_limit"     ,             required_argument, 0, 's'},
      {"output",                      required_argument, 0, 'o'},
      {0,                             0,                 0,  0}
    };

    int c = getopt_long(argc,
                        argv,
                        "lpnftso",
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

    case 'f':
      options->no_fork = true;
      break;

    case 't':
      options->max_num_threads = atoi(optarg);
      sleqp_log_debug("Using up to %d threads",
                      options->max_num_threads);

#ifdef WITH_OPENMP
      omp_set_num_threads(options->max_num_threads);
#endif
      break;

    case 's':
      options->time_limit = atof(optarg);
      sleqp_log_debug("Setting time limit to %ss",
                      optarg);
      break;

    case 'o':
      options->output = optarg;
      sleqp_log_debug("Setting output to %s",
                      optarg);
      break;
    default:
      sleqp_log_error("Invalid option");
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

static int
handle_cutest_child(const char* probname,
                    double time_limit,
                    pid_t child_pid)
{
  SleqpTimer* timer;

  sigset_t child_mask, old_mask;
  sigemptyset(&child_mask);
  sigaddset(&child_mask, SIGCHLD);

  SLEQP_CALL(sleqp_timer_create(&timer));

  if(sigprocmask(SIG_BLOCK, &child_mask, &old_mask) == -1)
  {
    sleqp_log_error("Failed to set signal mask: %s",
                    strerror(errno));

    return EXIT_FAILURE;
  }

  bool child_exited = false;

  const double total_time = time_limit + grace_seconds;
  int remaining_time = total_time;

  while(true)
  {
    const double elapsed = sleqp_timer_get_ttl(timer);
    remaining_time = total_time - elapsed;

    if(remaining_time <= 0)
    {
      break;
    }

    struct timespec ts = {0};
    ts.tv_sec = remaining_time;

    SLEQP_CALL(sleqp_timer_start(timer));

    int ret = sigtimedwait(&child_mask, NULL, &ts);

    SLEQP_CALL(sleqp_timer_stop(timer));

    if(ret > 0)
    {
      child_exited = true;
      break;
    }
    else
    {
      assert(ret == -1);

      // sigtimedwait() reached limit
      if(errno == EAGAIN)
      {
        break;
      }
      // sigtimedwait() was interrupted by signal
      else if(errno == EINTR)
      {
        continue;
      }
      else
      {
        sleqp_log_error("Failed to wait for child termination: %s",
                        strerror(errno));
      }

      break;
    }
  }

  if(!child_exited)
  {
    sleqp_log_error("Killing solver %s after deadline expired",
                    probname);

    kill(child_pid, SIGTERM);
  }

  SLEQP_CALL(sleqp_timer_free(&timer));

  return EXIT_SUCCESS;
}

static int
run_cutest_forking(SleqpCutestOptions *options)
{
  pid_t pid = fork();

  if(pid == -1)
  {
    sleqp_log_error("Failed to fork(): %s", strerror(errno));
    return EXIT_FAILURE;
  }
  else if(pid == 0)
  {
    // child
    return sleqp_cutest_run(filename, probname, options);
  }
  else
  {
    // parent
    assert(pid > 0);

    return handle_cutest_child(probname,
                               options->time_limit,
                               pid);
  }

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  SleqpCutestOptions options = (SleqpCutestOptions) {0};

  sleqp_cutest_options_default(&options);

  if(parse_command_line_options(argc, argv, &options) != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  if(options.no_fork)
  {
    return sleqp_cutest_run(filename, probname, &options);
  }
  else
  {
    return run_cutest_forking(&options);
  }

  return EXIT_SUCCESS;
}
