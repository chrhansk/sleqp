#include "types.h"

static const SleqpEnum active_enum = {.name    = "ActiveState",
                                      .flags   = true,
                                      .entries = {{"Inactive", SLEQP_INACTIVE},
                                                  {"Lower", SLEQP_ACTIVE_LOWER},
                                                  {"Upper", SLEQP_ACTIVE_UPPER},
                                                  {NULL, 0}}};

static const SleqpEnum status_enum
  = {.name    = "Status",
     .flags   = false,
     .entries = {{"Unknown", SLEQP_STATUS_UNKNOWN},
                 {"Running", SLEQP_STATUS_RUNNING},
                 {"Optimal", SLEQP_STATUS_OPTIMAL},
                 {"Infeasible", SLEQP_STATUS_INFEASIBLE},
                 {"Unbounded", SLEQP_STATUS_UNBOUNDED},
                 {"AbortDeadpoint", SLEQP_STATUS_ABORT_DEADPOINT},
                 {"AbortIter", SLEQP_STATUS_ABORT_ITER},
                 {"AbortManual", SLEQP_STATUS_ABORT_MANUAL},
                 {"AbortTime", SLEQP_STATUS_ABORT_TIME},
                 {NULL, 0}}};

static const SleqpEnum deriv_check_enum
  = {.name    = "DerivCheck",
     .flags   = true,
     .entries = {{"Skip", SLEQP_DERIV_CHECK_SKIP},
                 {"FirstObj", SLEQP_DERIV_CHECK_FIRST_OBJ},
                 {"FirstCons", SLEQP_DERIV_CHECK_FIRST_CONS},
                 {"SecondSimple", SLEQP_DERIV_CHECK_SECOND_SIMPLE},
                 {"SecondObj", SLEQP_DERIV_CHECK_SECOND_OBJ},
                 {"SecondCons", SLEQP_DERIV_CHECK_SECOND_CONS},
                 {NULL, 0}}};

static const SleqpEnum hess_eval_enum
  = {.name    = "HessEval",
     .flags   = false,
     .entries = {{"Exact", SLEQP_HESS_EVAL_EXACT},
                 {"SR1", SLEQP_HESS_EVAL_SR1},
                 {"SimpleBFGS", SLEQP_HESS_EVAL_SIMPLE_BFGS},
                 {"DampedBFGS", SLEQP_HESS_EVAL_DAMPED_BFGS},
                 {NULL, 0}}};

static const SleqpEnum bfgs_sizing_enum
  = {.name    = "BFGSSizing",
     .flags   = false,
     .entries = {{"None", SLEQP_BFGS_SIZING_NONE},
                 {"CenteredOL", SLEQP_BFGS_SIZING_CENTERED_OL},
                 {NULL, 0}}};

static const SleqpEnum steptype_enum
  = {.name    = "StepType",
     .flags   = false,
     .entries = {{"None", SLEQP_STEPTYPE_NONE},
                 {"Accepted", SLEQP_STEPTYPE_ACCEPTED},
                 {"AcceptedFull", SLEQP_STEPTYPE_ACCEPTED_FULL},
                 {"AcceptedSOC", SLEQP_STEPTYPE_ACCEPTED_SOC},
                 {"Rejected", SLEQP_STEPTYPE_REJECTED},
                 {NULL, 0}}};

static const SleqpEnum dual_estimation_enum
  = {.name    = "DualEstimation",
     .flags   = false,
     .entries = {{"LP", SLEQP_DUAL_ESTIMATION_TYPE_LP},
                 {"LSQ", SLEQP_DUAL_ESTIMATION_TYPE_LSQ},
                 {"Mixed", SLEQP_DUAL_ESTIMATION_TYPE_MIXED},
                 {NULL, 0}}};

static const SleqpEnum tr_solver_enum
  = {.name    = "TRSolver",
     .flags   = false,
     .entries = {{"Trlib", SLEQP_TR_SOLVER_TRLIB},
                 {"CG", SLEQP_TR_SOLVER_CG},
                 {"LSQR", SLEQP_TR_SOLVER_LSQR},
                 {"Auto", SLEQP_TR_SOLVER_AUTO},
                 {NULL, 0}}};

static const SleqpEnum polishing_enum
  = {.name    = "Polishing",
     .flags   = false,
     .entries = {{"None", SLEQP_POLISHING_NONE},
                 {"ZeroDual", SLEQP_POLISHING_ZERO_DUAL},
                 {"Inactive", SLEQP_POLISHING_INACTIVE},
                 {NULL, 0}}};

static const SleqpEnum parametric_cauchy_enum
  = {.name    = "ParametricCauchy",
     .flags   = false,
     .entries = {{"Disabled", SLEQP_PARAMETRIC_CAUCHY_DISABLED},
                 {"Coarse", SLEQP_PARAMETRIC_CAUCHY_COARSE},
                 {"Fine", SLEQP_PARAMETRIC_CAUCHY_FINE},
                 {NULL, 0}}};

static const SleqpEnum initial_tr_enum
  = {.name    = "InitialTRChoice",
     .flags   = false,
     .entries = {{"Narrow", SLEQP_INITIAL_TR_CHOICE_NARROW},
                 {"Wide", SLEQP_INITIAL_TR_CHOICE_WIDE},
                 {NULL, 0}}};

static const SleqpEnum linesearch_enum
  = {.name    = "Linesearch",
     .flags   = false,
     .entries = {{"Exact", SLEQP_LINESEARCH_EXACT},
                 {"Approx", SLEQP_LINESEARCH_APPROX},
                 {NULL, 0}}};

static const SleqpEnum step_rule_enum
  = {.name    = "StepRule",
     .flags   = false,
     .entries = {{"Direct", SLEQP_STEP_RULE_DIRECT},
                 {"Window", SLEQP_STEP_RULE_WINDOW},
                 {"Minstep", SLEQP_STEP_RULE_MINSTEP},
                 {NULL, 0}}};

const SleqpEnum*
sleqp_enum_active_state()
{
  return &active_enum;
}

const SleqpEnum*
sleqp_enum_status()
{
  return &status_enum;
}

const SleqpEnum*
sleqp_enum_deriv_check()
{
  return &deriv_check_enum;
}

const SleqpEnum*
sleqp_enum_hess_eval()
{
  return &hess_eval_enum;
}

const SleqpEnum*
sleqp_enum_bfgs_sizing()
{
  return &bfgs_sizing_enum;
}

const SleqpEnum*
sleqp_enum_steptype()
{
  return &steptype_enum;
}

const SleqpEnum*
sleqp_enum_dual_estimation()
{
  return &dual_estimation_enum;
}

const SleqpEnum*
sleqp_enum_tr_solver()
{
  return &tr_solver_enum;
}

const SleqpEnum*
sleqp_enum_polishing_type()
{
  return &polishing_enum;
}

const SleqpEnum*
sleqp_enum_parametric_cauchy()
{
  return &parametric_cauchy_enum;
}

const SleqpEnum*
sleqp_enum_initial_tr()
{
  return &initial_tr_enum;
}

const SleqpEnum*
sleqp_enum_linesearch()
{
  return &linesearch_enum;
}

const SleqpEnum*
sleqp_enum_step_rule()
{
  return &step_rule_enum;
}
