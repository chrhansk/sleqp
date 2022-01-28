#include <mex.h>

#include <assert.h>
#include <string.h>

#include "mex_fields.h"
#include "mex_info.h"
#include "mex_solve.h"
#include "mex_solve_lsq.h"

#include "sleqp.h"

#define COMMAND_BUFSIZE 256

typedef SLEQP_RETCODE (*MEX_FUNCTION)(int nlhs,
                                      mxArray* plhs[],
                                      int nrhs,
                                      const mxArray* prhs[]);

typedef struct
{
  const char* name;
  MEX_FUNCTION command;
  int nlhs;
  int nrhs;
} MexCommand;

static void
mex_log_handler(SLEQP_LOG_LEVEL level, time_t time, const char* message)
{
  const char* prefix = "";

  switch (level)
  {
  case SLEQP_LOG_ERROR:
    prefix = "error";
    break;
  case SLEQP_LOG_WARN:
    prefix = "warn";
    break;
  case SLEQP_LOG_INFO:
    prefix = "info";
    break;
  case SLEQP_LOG_DEBUG:
    prefix = "debug";
    break;
  default:
    break;
  }

  mexPrintf("%5s: %s\n", prefix, message);
}

const char*
error_prefix(SLEQP_ERROR_TYPE error_type)
{
  switch (error_type)
  {
  case SLEQP_NOMEM:
    return MEX_IDENTIFIER_NOMEM;
  case SLEQP_INTERNAL_ERROR:
    return MEX_IDENTIFIER_INTERNAL_ERROR;
  case SLEQP_FUNC_EVAL_ERROR:
    return MEX_IDENTIFIER_FUNC_EVAL_ERROR;
  case SLEQP_MATH_ERROR:
    return MEX_IDENTIFIER_MATH_ERROR;
  case SLEQP_INVALID_DERIV:
    return MEX_IDENTIFIER_INVALID_DERIV;
  case SLEQP_ILLEGAL_ARGUMENT:
    return MEX_IDENTIFIER_ILLEGAL_ARGUMENT;
  }

  return MEX_IDENTIFIER_DEFAULT;
}

SLEQP_EXPORT void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  char command_name[COMMAND_BUFSIZE];

  MexCommand commands[]
    = {{MEX_COMMAND_INFO, mex_command_info, 1, 0},
       {MEX_COMMAND_SOLVE, mex_command_solve, 2, 3},
       {MEX_COMMAND_SOLVE_LSQ, mex_command_solve_lsq, 2, 3}};

  const int num_commands = sizeof(commands) / sizeof(MexCommand);

  if (nrhs < 1)
  {
    mexErrMsgIdAndTxt(MEX_IDENTIFIER_DEFAULT, "Need at least one argument");

    return;
  }

  if (!mxIsChar(prhs[0]))
  {
    mexErrMsgIdAndTxt(MEX_IDENTIFIER_DEFAULT,
                      "First argument needs to be a string");

    return;
  }

  sleqp_log_set_handler(mex_log_handler);

  if (mxGetString(prhs[0], command_name, COMMAND_BUFSIZE) == 1)
  {
    mexErrMsgIdAndTxt(MEX_IDENTIFIER_DEFAULT, "Failed to extract command name");

    return;
  }

  for (int i = 0; i < num_commands; ++i)
  {
    if (!strcmp(command_name, commands[i].name))
    {
      --nrhs;
      ++prhs;

      if (commands[i].nlhs != nlhs)
      {
        mexErrMsgIdAndTxt(MEX_IDENTIFIER_DEFAULT,
                          "Invalid number of return values for function '%s', "
                          "expected %d, received %d",
                          command_name,
                          commands[i].nlhs,
                          nlhs);
        return;
      }

      if (commands[i].nrhs != nrhs)
      {
        mexErrMsgIdAndTxt(MEX_IDENTIFIER_DEFAULT,
                          "Invalid number of parameters for function '%s', "
                          "expected %d, received %d",
                          command_name,
                          commands[i].nrhs,
                          nrhs);
        return;
      }

      SLEQP_RETCODE retcode = commands[i].command(nlhs, plhs, nrhs, prhs);

      if (retcode != SLEQP_OKAY)
      {
        mexErrMsgIdAndTxt(error_prefix(sleqp_error_type()), sleqp_error_msg());
      }

      return;
    }
  }

  mexErrMsgIdAndTxt(MEX_IDENTIFIER_DEFAULT, "Invalid command %s", command_name);

  return;
}
