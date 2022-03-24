#include "mex_solve_common.h"

#include <assert.h>
#include <mex.h>

#include "mex_fields.h"
#include "mex_func_common.h"
#include "mex_output.h"
#include "mex_problem.h"

typedef struct
{
  const char* name;
  int value;
} Name;

static const Name param_names[] = {
  {MEX_PARAM_ZERO_EPS, SLEQP_PARAM_ZERO_EPS},
  {MEX_PARAM_EPS, SLEQP_PARAM_EPS},
  {MEX_PARAM_OBJ_LOWER, SLEQP_PARAM_OBJ_LOWER},
  {MEX_PARAM_DERIV_PERTURBATION, SLEQP_PARAM_DERIV_PERTURBATION},
  {MEX_PARAM_DERIV_TOL, SLEQP_PARAM_DERIV_TOL},
  {MEX_PARAM_CAUCHY_TAU, SLEQP_PARAM_CAUCHY_TAU},
  {MEX_PARAM_CAUCHY_ETA, SLEQP_PARAM_CAUCHY_ETA},
  {MEX_PARAM_LINESEARCH_TAU, SLEQP_PARAM_LINESEARCH_TAU},
  {MEX_PARAM_LINESEARCH_ETA, SLEQP_PARAM_LINESEARCH_ETA},
  {MEX_PARAM_LINESEARCH_CUTOFF, SLEQP_PARAM_LINESEARCH_CUTOFF},
  {MEX_PARAM_FEASIBILITY_TOL, SLEQP_PARAM_FEAS_TOL},
  {MEX_PARAM_SLACKNESS_TOL, SLEQP_PARAM_SLACK_TOL},
  {MEX_PARAM_STATIONARITY_TOL, SLEQP_PARAM_STAT_TOL},
  {MEX_PARAM_ACCEPTED_REDUCTION, SLEQP_PARAM_ACCEPTED_REDUCTION},
  {MEX_PARAM_DEADPOINT_BOUND, SLEQP_PARAM_DEADPOINT_BOUND},
};

static const Name enum_option_names[]
  = {{MEX_DERIV_CHECK, SLEQP_OPTION_ENUM_DERIV_CHECK},
     {MEX_HESS_EVAL, SLEQP_OPTION_ENUM_HESS_EVAL},
     {MEX_DUAL_ESTIMATION_TYPE, SLEQP_OPTION_ENUM_DUAL_ESTIMATION_TYPE},
     {MEX_BFGS_SIZING, SLEQP_OPTION_ENUM_BFGS_SIZING},
     {MEX_TR_SOLVER, SLEQP_OPTION_ENUM_TR_SOLVER},
     {MEX_POLISHING_TYPE, SLEQP_OPTION_ENUM_POLISHING_TYPE},
     {MEX_STEP_RULE, SLEQP_OPTION_ENUM_STEP_RULE},
     {MEX_LINESEARCH, SLEQP_OPTION_ENUM_LINESEARCH},
     {MEX_PARAMETRIC_CAUCHY, SLEQP_OPTION_ENUM_PARAMETRIC_CAUCHY},
     {MEX_INITIAL_TR_CHOICE, SLEQP_OPTION_ENUM_INITIAL_TR_CHOICE}};

static const Name int_option_names[] = {
  {MEX_NUM_QUASI_NEWTON_ITERATES, SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES},
  {MEX_MAX_NEWTON_ITERATIONS, SLEQP_OPTION_INT_MAX_NEWTON_ITERATIONS},
  {MEX_NUM_THREADS, SLEQP_OPTION_INT_NUM_THREADS}};

static const Name bool_option_names[]
  = {{MEX_PERFORM_NEWTON_STEP, SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP},
     {MEX_PERFORM_SOC, SLEQP_OPTION_BOOL_PERFORM_SOC},
     {MEX_USE_QUADRATIC_MODEL, SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL},
     {MEX_ALWAYS_WARM_START_LP, SLEQP_OPTION_BOOL_ALWAYS_WARM_START_LP},
     {MEX_ENABLE_RESTORATION_PHASE, SLEQP_OPTION_BOOL_ENABLE_RESTORATION_PHASE},
     {MEX_ENABLE_PREPROCESSOR, SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR}};

static SLEQP_RETCODE
read_option_entry(const mxArray* mex_options,
                  const char* name,
                  double* scalar_value,
                  bool* present)
{
  const mxArray* value = mxGetField(mex_options, 0, name);

  if (!value)
  {
    *present = false;
    return SLEQP_OKAY;
  }

  *present = true;

  if (!(mxIsScalar(value) && mxIsDouble(value)))
  {
    return SLEQP_ERROR;
  }

  const double* ptr = mxGetPr(value);

  assert(ptr);

  *scalar_value = *ptr;

  return SLEQP_OKAY;
}

typedef SLEQP_RETCODE (*SET_VALUE)(int name, double value, void* data);

static SLEQP_RETCODE
read_values(const mxArray* mex_options,
            const Name* names,
            int num_names,
            SET_VALUE set_value,
            void* data)
{
  for (int i = 0; i < num_names; ++i)
  {
    const Name* name = names + i;
    double param_value;
    bool present;

    SLEQP_CALL(
      read_option_entry(mex_options, name->name, &param_value, &present));

    if (present)
    {
      SLEQP_CALL(set_value(name->value, param_value, data));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_param_value(int name, double value, void* data)
{
  return sleqp_params_set_value((SleqpParams*)data, name, value);
}

static SLEQP_RETCODE
read_params(SleqpParams* params, const mxArray* mex_options)
{
  if (!mxIsStruct(mex_options))
  {
    return SLEQP_ERROR;
  }

  const int num_params = sizeof(param_names) / sizeof(param_names[0]);

  SLEQP_CALL(
    read_values(mex_options, param_names, num_params, set_param_value, params));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_option_bool_value(int name, double value, void* data)
{
  return sleqp_options_set_bool_value((SleqpOptions*)data, name, !!(value));
}

static SLEQP_RETCODE
set_option_int_value(int name, double value, void* data)
{
  return sleqp_options_set_int_value((SleqpOptions*)data, name, (int)value);
}

static SLEQP_RETCODE
set_option_enum_value(int name, double value, void* data)
{
  return sleqp_options_set_enum_value((SleqpOptions*)data, name, (int)value);
}

static SLEQP_RETCODE
read_options(SleqpOptions* options, const mxArray* mex_options)
{
  if (!mxIsStruct(mex_options))
  {
    return SLEQP_ERROR;
  }

  const int num_bool_options
    = sizeof(bool_option_names) / sizeof(bool_option_names[0]);

  SLEQP_CALL(read_values(mex_options,
                         bool_option_names,
                         num_bool_options,
                         set_option_bool_value,
                         options));

  const int num_int_options
    = sizeof(int_option_names) / sizeof(int_option_names[0]);

  SLEQP_CALL(read_values(mex_options,
                         int_option_names,
                         num_int_options,
                         set_option_int_value,
                         options));

  const int num_enum_options
    = sizeof(enum_option_names) / sizeof(enum_option_names[0]);

  SLEQP_CALL(read_values(mex_options,
                         enum_option_names,
                         num_enum_options,
                         set_option_enum_value,
                         options));

  return SLEQP_OKAY;
}

typedef struct
{
  mxArray* primal;

  // Callbacks
  struct
  {
    mxArray* accepted_iterate;
  } callbacks;

} CallbackData;

static SLEQP_RETCODE
accepted_iterate(SleqpSolver* solver, SleqpIterate* iterate, void* data)
{
  CallbackData* callback_data = (CallbackData*)data;

  SleqpVec* primal = sleqp_iterate_primal(iterate);

  SLEQP_CALL(sleqp_vec_to_raw(primal, mxGetPr(callback_data->primal)));

  mxArray* rhs[]
    = {callback_data->callbacks.accepted_iterate, callback_data->primal};

  const int nrhs = sizeof(rhs) / sizeof(rhs[0]);

  bool value = false;

  SLEQP_CALL(mex_eval_into_bool(nrhs, rhs, &value));

  if (value)
  {
    SLEQP_CALL(sleqp_solver_abort(solver));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_callback_data(SleqpProblem* problem, CallbackData* callback_data)
{
  const int num_vars = sleqp_problem_num_vars(problem);

  callback_data->primal = mxCreateDoubleMatrix(num_vars, 1, mxREAL);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
register_callbacks(SleqpProblem* problem,
                   SleqpSolver* solver,
                   const mxArray* mex_funcs,
                   CallbackData* callback_data)
{
  bool has_callbacks = false;

  bool has_accepted_iterate;

  SLEQP_CALL(mex_callback_has_field(mex_funcs,
                                    MEX_CALLBACK_ACCEPTED_ITERATE,
                                    &has_accepted_iterate));

  if (has_accepted_iterate)
  {
    SLEQP_CALL(
      mex_callback_from_struct(mex_funcs,
                               MEX_CALLBACK_ACCEPTED_ITERATE,
                               &callback_data->callbacks.accepted_iterate));

    SLEQP_CALL(sleqp_solver_add_callback(solver,
                                         SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                                         accepted_iterate,
                                         (void*)callback_data));

    has_callbacks = true;
  }

  if (has_callbacks)
  {
    SLEQP_CALL(create_callback_data(problem, callback_data));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
destroy_callback_data(CallbackData* callback_data)
{
  mxDestroyArray(callback_data->primal);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_solve(mxArray** sol_star,
          mxArray** info_star,
          bool lsq,
          const mxArray* mex_x0,
          const mxArray* mex_funcs,
          const mxArray* mex_options)
{
  SleqpOptions* options;
  SleqpParams* params;
  SleqpProblem* problem;
  SleqpSolver* solver;
  SleqpVec* initial;

  CallbackData callback_data = (CallbackData){0};

  SLEQP_CALL(sleqp_options_create(&options));
  SLEQP_CALL(sleqp_params_create(&params));

  SLEQP_CALL(read_params(params, mex_options));

  SLEQP_CALL(read_options(options, mex_options));

  SLEQP_CALL(mex_problem_create(&problem,
                                params,
                                options,
                                lsq,
                                mex_x0,
                                mex_funcs,
                                mex_options));

  SLEQP_CALL(mex_create_vec_from_array(&initial, mex_x0));

  SLEQP_CALL(
    sleqp_solver_create(&solver, problem, params, options, initial, NULL));

  SLEQP_CALL(register_callbacks(problem, solver, mex_funcs, &callback_data));

  SLEQP_CALL(sleqp_solver_solve(solver, SLEQP_NONE, SLEQP_NONE));

  SLEQP_CALL(mex_create_solver_output(problem, solver, sol_star, info_star));

  SLEQP_CALL(destroy_callback_data(&callback_data));

  SLEQP_CALL(sleqp_vec_free(&initial));
  SLEQP_CALL(sleqp_solver_release(&solver));
  SLEQP_CALL(sleqp_problem_release(&problem));

  SLEQP_CALL(sleqp_params_release(&params));
  SLEQP_CALL(sleqp_options_release(&options));

  return SLEQP_OKAY;
}
