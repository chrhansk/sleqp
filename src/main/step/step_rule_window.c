#include "step_rule_window.h"

#include <assert.h>

#include "cmp.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

typedef struct
{
  // max size
  int num;
  // curr size
  int len;
  // curr index
  int curr;

  double* exact_merit_values;
  double* model_reductions;

  SleqpParams* params;

} StepRule;

static int
data_index(const StepRule* step_rule, int index)
{
  if (step_rule->len == 0)
  {
    return 0;
  }

  int data_index = index % (step_rule->num);

  return (data_index < 0) ? (data_index + step_rule->num) : data_index;
}

static int
data_begin(const StepRule* step_rule)
{
  return step_rule->curr - step_rule->len + 1;
}

static int
data_end(const StepRule* step_rule)
{
  return step_rule->curr;
}

static int
compute_ref_index(const StepRule* step_rule)
{
  const int begin = data_begin(step_rule);
  const int end   = data_end(step_rule);

  assert(step_rule->len > 0);

  double max_merit_value = -sleqp_infinity();
  int max_index          = 0;

  for (int i = begin; i <= end; ++i)
  {
    int j                    = data_index(step_rule, i);
    double exact_merit_value = step_rule->exact_merit_values[j];

    if (exact_merit_value > max_merit_value)
    {
      max_index       = i;
      max_merit_value = exact_merit_value;
    }
  }

  return max_index;
}

static double
compute_historic_reduction_ratio(const StepRule* step_rule,
                                 const double iterate_merit,
                                 const double trial_exact_merit,
                                 const double model_reduction)
{
  const double current_reduction_ratio
    = (iterate_merit - trial_exact_merit) / model_reduction;

  if (step_rule->len == 0)
  {
    return current_reduction_ratio;
  }

  const int ref_index = compute_ref_index(step_rule);
  const int end       = data_end(step_rule);

  const int j            = data_index(step_rule, ref_index);
  const double ref_merit = step_rule->exact_merit_values[j];

  // Maximum is attained at current iterate,
  // historic reduction ratio == current reduction ratio
  if (ref_merit < trial_exact_merit)
  {
    return current_reduction_ratio;
  }

  assert(data_begin(step_rule) <= ref_index);
  assert(ref_index <= end);

  double model_reduction_sum = 0.;

  for (int i = ref_index; i <= end; ++i)
  {
    int j = data_index(step_rule, i);

    model_reduction_sum += step_rule->model_reductions[j];
  }

  const double total_model_reduction = model_reduction_sum + model_reduction;

  return (ref_merit - trial_exact_merit) / total_model_reduction;
}

static SLEQP_RETCODE
add_accepted_step(StepRule* step_rule,
                  double iterate_merit,
                  double model_reduction)
{
  if (step_rule->num == 0)
  {
    return SLEQP_OKAY;
  }

  const double eps = sleqp_params_value(step_rule->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  sleqp_assert_is_geq(model_reduction, 0., eps);

  int next = data_end(step_rule) + 1;

  int j = data_index(step_rule, next);

  step_rule->exact_merit_values[j] = iterate_merit;
  step_rule->model_reductions[j]   = model_reduction;

  if (step_rule->len < step_rule->num)
  {
    ++step_rule->len;
  }

  ++(step_rule->curr);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
step_rule_window_apply(double iterate_merit,
                       double trial_exact_merit,
                       double trial_model_merit,
                       bool* accept_step,
                       double* reduction_ratio,
                       void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  const double exact_reduction = iterate_merit - trial_exact_merit;
  const double model_reduction = iterate_merit - trial_model_merit;

  double current_reduction_ratio = 1.;

  if (exact_reduction != model_reduction)
  {
    current_reduction_ratio = exact_reduction / model_reduction;
  }

  const double historic_reduction_ratio
    = compute_historic_reduction_ratio(step_rule,
                                       iterate_merit,
                                       trial_exact_merit,
                                       model_reduction);

  *reduction_ratio
    = SLEQP_MAX(current_reduction_ratio, historic_reduction_ratio);

  const double accepted_reduction
    = sleqp_params_value(step_rule->params, SLEQP_PARAM_ACCEPTED_REDUCTION);

  *accept_step = (*reduction_ratio >= accepted_reduction);

  sleqp_log_debug("Step with current reduction ratio %e, historic reduction "
                  "ratio %e, accepted: %s",
                  current_reduction_ratio,
                  historic_reduction_ratio,
                  sleqp_bool_string(*accept_step));

  if (*accept_step)
  {
    SLEQP_CALL(add_accepted_step(step_rule, iterate_merit, model_reduction));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
step_rule_window_reset(void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  step_rule->curr = -1;
  step_rule->len  = 0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
step_rule_window_free(void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  SLEQP_CALL(sleqp_params_release(&step_rule->params));

  sleqp_free(&step_rule->model_reductions);
  sleqp_free(&step_rule->exact_merit_values);

  sleqp_free(&step_rule);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_step_rule_window_create(SleqpStepRule** star,
                              SleqpProblem* problem,
                              SleqpParams* params,
                              int window_size)
{
  assert(window_size >= 0);

  SleqpStepRuleCallbacks callbacks = {.rule_apply = step_rule_window_apply,
                                      .rule_reset = step_rule_window_reset,
                                      .rule_free  = step_rule_window_free};

  StepRule* step_rule = NULL;

  SLEQP_CALL(sleqp_malloc(&step_rule));

  step_rule->num = window_size;
  SLEQP_CALL(step_rule_window_reset(step_rule));

  SLEQP_CALL(sleqp_alloc_array(&step_rule->exact_merit_values, window_size));
  SLEQP_CALL(sleqp_alloc_array(&step_rule->model_reductions, window_size));

  SLEQP_CALL(sleqp_params_capture(params));
  step_rule->params = params;

  SLEQP_CALL(sleqp_step_rule_create(star,
                                    problem,
                                    params,
                                    &callbacks,
                                    (void*)step_rule));

  return SLEQP_OKAY;
}
