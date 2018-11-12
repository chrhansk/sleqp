#include <cutest.h>

#include "sleqp.h"

int main(int argc, char *argv[])
{
  //char *fname = "OUTSDIF.d"; /* CUTEst data file */

  char *fname = PROBLEM_OUTSDIF_D; /* CUTEst data file */

  integer funit = 42;   /* FORTRAN unit number for OUTSDIF.d */
  integer ierr;         /* Exit flag from OPEN and CLOSE */
  integer status;       /* Exit flag from CUTEst tools */

  integer CUTEst_nvar;        /* number of variables */
  integer CUTEst_ncon;        /* number of constraints */
  integer CUTEst_lcjac;       /* length of Jacobian arrays */
  integer CUTEst_nnzj;        /* number of nonzeros in Jacobian */
  integer CUTEst_nnzh;        /* number of nonzeros in upper triangular
                                 part of the Hessian of the Lagrangian */
  ierr = 0;
  FORTRAN_open(&funit, fname, &ierr);

  if(ierr)
  {
    sleqp_log_error("Failed to open %s, aborting.", fname);
    return 1;
  }

  CUTEST_cdimen( &status, &funit, &CUTEst_nvar, &CUTEst_ncon);

  sleqp_log_info("Problem has %d variables, %d constraints",
                 CUTEst_nvar,
                 CUTEst_ncon);

  FORTRAN_close(&funit, &ierr);

  if(ierr)
  {
    sleqp_log_error("Error closing %s on unit %d.\n",
                    fname,
                    funit);
  }

  return 0;
}
