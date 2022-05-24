#include "mpi_utils.h"

#include <mpi.h>

#include "error.h"
#include "log.h"
#include "pub_types.h"

#define SLEQP_MPI_CALL(x)                                                      \
  do                                                                           \
  {                                                                            \
    const int _status = (x);                                                   \
    if (_status != MPI_SUCCESS)                                                \
    {                                                                          \
      sleqp_raise(SLEQP_INTERNAL_ERROR, "MPI error in function %s", __func__); \
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
