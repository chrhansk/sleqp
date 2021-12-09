#include "mpi_utils.h"

#include <mpi.h>

#include "log.h"

#define SLEQP_MPI_CALL(x)                                                      \
  do                                                                           \
  {                                                                            \
    const int _status = (x);                                                   \
    if (_status != MPI_SUCCESS)                                                \
    {                                                                          \
      sleqp_log_error("MPI error in function %s", __func__);                   \
      return SLEQP_INTERNAL_ERROR;                                             \
    }                                                                          \
  } while (0)

SLEQP_RETCODE
sleqp_mpi_initialize()
{
  int flag;

  SLEQP_MPI_CALL(MPI_Initialized(&flag));

  if (!flag)
  {
    sleqp_log_debug("Manually initializing MPI");

    SLEQP_MPI_CALL(MPI_Init(NULL, NULL));
  }

  return SLEQP_OKAY;
}
