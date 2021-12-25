#include <mex.h>

#include <assert.h>
#include <string.h>

#include "mex_fields.h"
#include "mex_solve.h"

#include "sleqp.h"

#define COMMAND_BUFSIZE 256

// TODO: Better logging
// TODO: Pass more options / params

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
  mexPrintf("%s\n", message);
}

void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  char command_name[COMMAND_BUFSIZE];

  MexCommand commands[] = {{MEX_COMMAND_SOLVE, mex_command_solve, 2, 3}};

  const int num_commands = 1;

  if (nrhs < 1)
  {
    mexErrMsgIdAndTxt(MEX_MSG_IDENTIFIER, "Need at least one argument");

    return;
  }

  if (!mxIsChar(prhs[0]))
  {
    mexErrMsgIdAndTxt(MEX_MSG_IDENTIFIER,
                      "First argument needs to be a string");

    return;
  }

  sleqp_log_set_handler(mex_log_handler);

  if (mxGetString(prhs[0], command_name, COMMAND_BUFSIZE) == 1)
  {
    mexErrMsgIdAndTxt(MEX_MSG_IDENTIFIER, "Failed to extract command name");

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
        mexErrMsgIdAndTxt(MEX_MSG_IDENTIFIER,
                          "Invalid number of return values for function '%s', "
                          "expected %d, received %d",
                          command_name,
                          commands[i].nlhs,
                          nlhs);
        return;
      }

      if (commands[i].nrhs != nrhs)
      {
        mexErrMsgIdAndTxt(MEX_MSG_IDENTIFIER,
                          "Invalid number of parameters for function '%s', "
                          "expected %d, received %d",
                          command_name,
                          commands[i].nrhs,
                          nrhs);
        return;
      }

      SLEQP_STATUS status = commands[i].command(nlhs, plhs, nrhs, prhs);

      if (status != SLEQP_OKAY)
      {
        mexErrMsgIdAndTxt(MEX_MSG_IDENTIFIER,
                          "Error executing command %s",
                          command_name);
      }

      return;
    }
  }

  mexErrMsgIdAndTxt(MEX_MSG_IDENTIFIER, "Invalid command %s", command_name);

  return;
}
