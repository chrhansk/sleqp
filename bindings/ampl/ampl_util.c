#include "ampl_util.h"

bool
sleqp_ampl_max_problem(ASL* asl)
{
  if (n_obj > 0 && objtype[obj_no] != 0)
  {
    return true;
  }

  return false;
}
