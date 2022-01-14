#include "mex_info.h"

#include "mex_fields.h"

#include <assert.h>

static SLEQP_RETCODE
set_struct_field_to_string(mxArray* info, const char* name, const char* value)
{
  mxSetField(info, 0, name, mxCreateString(value));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_struct_field_to_real(mxArray* info, const char* name, double value)
{
  mxArray* array;

  array = mxCreateDoubleScalar(value);

  mxSetField(info, 0, name, array);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_command_info(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  assert(nlhs == 1);
  assert(nrhs == 0);

  const SLEQP_LOG_LEVEL level = sleqp_log_level();

  sleqp_log_set_level(SLEQP_LOG_SILENT);

  // plhs[0] = mxCreateString(str);

  sleqp_log_set_level(level);

  const char* fieldnames[] = {MEX_INFO_VERSION,
                              MEX_INFO_VERSION_MAJOR,
                              MEX_INFO_VERSION_MINOR,
                              MEX_INFO_VERSION_PATCH,
                              MEX_INFO_FACT_NAME,
                              MEX_INFO_FACT_VERSION,
                              MEX_INFO_LPS_NAME,
                              MEX_INFO_LPS_VERSION};

  const int num_fields = sizeof(fieldnames) / sizeof(const char*);

  plhs[0] = mxCreateStructMatrix(1, 1, num_fields, fieldnames);

  mxArray* info_array = plhs[0];

  SLEQP_CALL(
    set_struct_field_to_string(info_array, MEX_INFO_VERSION, SLEQP_VERSION));

  SLEQP_CALL(set_struct_field_to_real(info_array,
                                      MEX_INFO_VERSION_MAJOR,
                                      SLEQP_VERSION_MAJOR));

  SLEQP_CALL(set_struct_field_to_real(info_array,
                                      MEX_INFO_VERSION_MINOR,
                                      SLEQP_VERSION_MINOR));

  SLEQP_CALL(set_struct_field_to_real(info_array,
                                      MEX_INFO_VERSION_PATCH,
                                      SLEQP_VERSION_PATCH));

  SLEQP_CALL(set_struct_field_to_string(info_array,
                                        MEX_INFO_FACT_NAME,
                                        SLEQP_FACT_NAME));

  SLEQP_CALL(set_struct_field_to_string(info_array,
                                        MEX_INFO_FACT_VERSION,
                                        SLEQP_FACT_VERSION));

  SLEQP_CALL(set_struct_field_to_string(info_array,
                                        MEX_INFO_LPS_NAME,
                                        SLEQP_LP_SOLVER_NAME));

  SLEQP_CALL(set_struct_field_to_string(info_array,
                                        MEX_INFO_LPS_VERSION,
                                        SLEQP_LP_SOLVER_VERSION));

  return SLEQP_OKAY;
}
