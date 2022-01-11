#include "ampl_keywords.h"

#include "ampl_mem.h"

typedef struct
{
  union
  {
    SleqpOptions* options;
    SleqpParams* params;
    SleqpAmplKeywords* keywords;
  } opt_params;

  int index;

} CallbackData;

enum
{
  ITER_LIMIT_MAXITER,
  ITER_LIMIT_MAXIT,
  TIME_LIMIT,
  LOG_LEVEL_PRINT_LEVEL,
  LOG_LEVEL_OUTLEV,
  WANTSOL,
  NUM_EXTRA
};

typedef enum
{
  POS_ENUM          = 0,
  POS_INT           = SLEQP_NUM_ENUM_OPTIONS,
  POS_BOOL          = POS_INT + SLEQP_NUM_INT_OPTIONS,
  POS_PAR           = POS_BOOL + SLEQP_NUM_BOOL_OPTIONS,
  POS_EXTRA         = POS_PAR + SLEQP_NUM_PARAMS,
  AMPL_NUM_KEYWORDS = POS_EXTRA + NUM_EXTRA
} AMPL_KEYWORDS;

static char*
kwdfunc_enum(Option_Info* oi, keyword* kw, char* value)
{
  sleqp_log_debug("Inside kwdfunc_enum, kw: %s, value as char: '%s'",
                  kw->name,
                  value);

  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  SLEQP_RETCODE retcode
    = sleqp_options_set_enum_value(callback_data->opt_params.options,
                                   callback_data->index,
                                   int_val);

  if (retcode != SLEQP_OKAY)
  {
    return badval_ASL(oi, kw, value, retval);
  }

  return retval;
}

static char*
kwdfunc_int(Option_Info* oi, keyword* kw, char* value)
{
  sleqp_log_debug("Inside kwdfunc_int, kw: %s, value as char: '%s'",
                  kw->name,
                  value);

  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  SLEQP_RETCODE retcode
    = sleqp_options_set_int_value(callback_data->opt_params.options,
                                  callback_data->index,
                                  int_val);

  if (retcode != SLEQP_OKAY)
  {
    return badval_ASL(oi, kw, value, retval);
  }

  return retval;
}

static char*
kwdfunc_bool(Option_Info* oi, keyword* kw, char* value)
{
  sleqp_log_debug("Inside kwdfunc_bool, kw: %s, value as char: '%s'",
                  kw->name,
                  value);

  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = IK1_val(oi, kw, value);

  SLEQP_RETCODE retcode
    = sleqp_options_set_bool_value(callback_data->opt_params.options,
                                   callback_data->index,
                                   true);

  if (retcode != SLEQP_OKAY)
  {
    return badval_ASL(oi, kw, value, retval);
  }

  return retval;
}

static char*
kwdfunc_param(Option_Info* oi, keyword* kw, char* value)
{
  sleqp_log_debug("Inside kwdfunc_param, kw: %s, value as char: '%s'",
                  kw->name,
                  value);

  CallbackData* callback_data = (CallbackData*)kw->info;

  double real_val;
  kw->info = &real_val;

  char* retval = D_val(oi, kw, value);

  SLEQP_RETCODE retcode
    = sleqp_params_set_value(callback_data->opt_params.params,
                             callback_data->index,
                             real_val);

  if (retcode != SLEQP_OKAY)
  {
    return badval_ASL(oi, kw, value, retval);
  }

  return retval;
}

struct SleqpAmplKeywords
{
  SleqpOptions* options;
  SleqpParams* params;

  CallbackData callback_data[AMPL_NUM_KEYWORDS];
  keyword keywds[AMPL_NUM_KEYWORDS];

  double time_limit;
  int iteration_limit;
};

static char*
kwdfunc_iterlimit(Option_Info* oi, keyword* kw, char* value)
{
  sleqp_log_debug("Inside kwdfunc_iterlimit, kw: %s, value as char: '%s'",
                  kw->name,
                  value);

  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  if (int_val != SLEQP_NONE)
  {
    if (int_val < 0)
    {
      return badval_ASL(oi, kw, value, retval);
    }

    callback_data->opt_params.keywords->iteration_limit = int_val;
  }

  return retval;
}

static char*
kwdfunc_timelimit(Option_Info* oi, keyword* kw, char* value)
{
  sleqp_log_debug("Inside kwdfunc_timelimit, kw: %s, value as char: '%s'",
                  kw->name,
                  value);

  CallbackData* callback_data = (CallbackData*)kw->info;

  double real_val;
  kw->info = &real_val;

  char* retval = D_val(oi, kw, value);

  if (real_val != SLEQP_NONE)
  {
    if (real_val < 0)
    {
      return badval_ASL(oi, kw, value, retval);
    }

    callback_data->opt_params.keywords->time_limit = real_val;
  }

  return retval;
}

static char*
kwdfunc_log_level(Option_Info* oi, keyword* kw, char* value)
{
  sleqp_log_debug("Inside kwdfunc_log_level, kw: %s, value as char: '%s'",
                  kw->name,
                  value);

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  if (int_val != SLEQP_NONE)
  {
    if (int_val < SLEQP_LOG_SILENT || int_val > SLEQP_LOG_DEBUG)
    {
      return badval_ASL(oi, kw, value, retval);
    }

    sleqp_log_set_level(int_val);
  }

  return retval;
}

static int
compare_kwds(const void* first, const void* second)
{
  return strcmp(((keyword*)first)->name, ((keyword*)second)->name);
}

static SLEQP_RETCODE
keywords_fill(SleqpAmplKeywords* ampl_keywords,
              SleqpOptions* options,
              SleqpParams* params,
              keyword* kwds)
{
  CallbackData* callback_data = ampl_keywords->callback_data;

  int pos = POS_ENUM;

  for (; pos < POS_INT; ++pos)
  {
    callback_data[pos]
      = (CallbackData){.opt_params.options = options, .index = pos};

    kwds[pos] = (keyword){.name = strdup(sleqp_options_enum_name(pos)),
                          .kf   = kwdfunc_enum,
                          .info = callback_data + pos,
                          .desc = "Description"};
  }

  for (; pos < POS_BOOL; ++pos)
  {
    callback_data[pos]
      = (CallbackData){.opt_params.options = options, .index = pos - POS_INT};

    kwds[pos] = (keyword){.name = strdup(sleqp_options_int_name(pos - POS_INT)),
                          .kf   = kwdfunc_int,
                          .info = callback_data + pos,
                          .desc = "Description"};
  }

  for (; pos < POS_PAR; ++pos)
  {
    callback_data[pos]
      = (CallbackData){.opt_params.options = options, .index = pos - POS_BOOL};

    kwds[pos]
      = (keyword){.name = strdup(sleqp_options_bool_name(pos - POS_BOOL)),
                  .kf   = kwdfunc_bool,
                  .info = callback_data + pos,
                  .desc = "Description"};
  }

  for (; pos < POS_EXTRA; ++pos)
  {
    SLEQP_PARAM param_value = pos - POS_PAR;

    callback_data[pos]
      = (CallbackData){.opt_params.params = params, .index = param_value};

    kwds[pos] = (keyword){.name = strdup(sleqp_params_name(param_value)),
                          .kf   = kwdfunc_param,
                          .info = callback_data + pos,
                          .desc = strdup(sleqp_params_desc(param_value))};
  }

  for (pos = POS_EXTRA; pos < AMPL_NUM_KEYWORDS; ++pos)
  {
    callback_data[pos]
      = (CallbackData){.opt_params.keywords = ampl_keywords, .index = 0};
  }

  {
    pos = POS_EXTRA + ITER_LIMIT_MAXITER;

    kwds[pos] = (keyword){.name = "max_iter",
                          .kf   = kwdfunc_iterlimit,
                          .info = callback_data + pos,
                          .desc = "Maximum number of iterations"};
  }

  {
    pos = POS_EXTRA + ITER_LIMIT_MAXIT;

    kwds[pos] = (keyword){.name = "maxit",
                          .kf   = kwdfunc_iterlimit,
                          .info = callback_data + pos,
                          .desc = "Alias for 'max_iter'"};
  }

  {
    pos = POS_EXTRA + TIME_LIMIT;

    kwds[pos] = (keyword){.name = "max_wall_time",
                          .kf   = kwdfunc_timelimit,
                          .info = callback_data + pos,
                          .desc = "Wallclock time limit"};
  }

  {
    pos = POS_EXTRA + LOG_LEVEL_PRINT_LEVEL;

    kwds[pos] = (keyword){.name = "print_level",
                          .kf   = kwdfunc_log_level,
                          .info = NULL,
                          .desc = "Verbosity level"};
  }

  {
    pos = POS_EXTRA + LOG_LEVEL_OUTLEV;

    kwds[pos] = (keyword){.name = "outlev",
                          .kf   = kwdfunc_log_level,
                          .info = NULL,
                          .desc = "Alias for 'print_level'"};
  }

  {
    pos = POS_EXTRA + WANTSOL;

    kwds[pos] = (keyword){
      .name = "wantsol",
      .kf   = WS_val,
      .info = NULL,
      .desc = "solution report without -AMPL: sum of 1 (write .sol file), 2 "
              "(print primal variable values), 4 (print dual variable values), "
              "8 (do not print solution message)"};
  }

  // Keywords must be sorted alphabetically
  qsort(kwds, AMPL_NUM_KEYWORDS, sizeof(keyword), compare_kwds);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_keywords_create(SleqpAmplKeywords** star,
                           SleqpOptions* options,
                           SleqpParams* params)
{
  SLEQP_CALL(sleqp_ampl_malloc(star));

  SleqpAmplKeywords* ampl_keywords = *star;

  (*ampl_keywords) = (SleqpAmplKeywords){0};

  ampl_keywords->time_limit      = SLEQP_NONE;
  ampl_keywords->iteration_limit = SLEQP_NONE;

  SLEQP_CALL(sleqp_options_capture(options));
  ampl_keywords->options = options;

  SLEQP_CALL(sleqp_params_capture(params));
  ampl_keywords->params = params;

  SLEQP_CALL(
    keywords_fill(ampl_keywords, options, params, ampl_keywords->keywds));

  return SLEQP_OKAY;
}

double
sleqp_ampl_keywords_iter_limit(SleqpAmplKeywords* ampl_keywords)
{
  return ampl_keywords->iteration_limit;
}

double
sleqp_ampl_keywords_time_limit(SleqpAmplKeywords* ampl_keywords)
{
  return ampl_keywords->time_limit;
}

SLEQP_RETCODE
sleqp_ampl_keywords_get(SleqpAmplKeywords* ampl_keywords,
                        keyword** star,
                        int* num_keywords)
{
  *star = ampl_keywords->keywds;

  *num_keywords = AMPL_NUM_KEYWORDS;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_keywords_free(SleqpAmplKeywords** star)
{
  SleqpAmplKeywords* sleqp_keywords = *star;

  if (!sleqp_keywords)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_params_release(&sleqp_keywords->params));

  SLEQP_CALL(sleqp_options_release(&sleqp_keywords->options));

  sleqp_free(star);

  return SLEQP_OKAY;
}
