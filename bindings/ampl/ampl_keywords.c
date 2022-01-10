#include "ampl_keywords.h"

#include "ampl_mem.h"

typedef struct
{
  union
  {
    SleqpOptions* options;
    SleqpParams* params;
  } opt_params;

  int index;

} CallbackData;

typedef enum
{
  POS_ENUM          = 0,
  POS_INT           = SLEQP_NUM_ENUM_OPTIONS,
  POS_BOOL          = POS_INT + SLEQP_NUM_INT_OPTIONS,
  POS_PAR           = POS_BOOL + SLEQP_NUM_BOOL_OPTIONS,
  AMPL_NUM_KEYWORDS = POS_PAR + SLEQP_NUM_PARAMS
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

static int
compare_kwds(const void* first, const void* second)
{
  return strcmp(((keyword*)first)->name, ((keyword*)second)->name);
}

static SLEQP_RETCODE
keywords_fill(SleqpOptions* options,
              SleqpParams* params,
              keyword* kwds,
              CallbackData* callback_data)
{
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

  for (; pos < AMPL_NUM_KEYWORDS; ++pos)
  {
    callback_data[pos]
      = (CallbackData){.opt_params.params = params, .index = pos - POS_PAR};

    kwds[pos] = (keyword){.name = strdup(sleqp_params_name(pos - POS_PAR)),
                          .kf   = kwdfunc_param,
                          .info = callback_data + pos,
                          .desc = "Description"};
  }

  // Keywords must be sorted alphabetically
  qsort(kwds, AMPL_NUM_KEYWORDS, sizeof(keyword), compare_kwds);

  return SLEQP_OKAY;
}

struct SleqpAmplKeywords
{
  SleqpOptions* options;
  SleqpParams* params;

  CallbackData callback_data[AMPL_NUM_KEYWORDS];
  keyword keywds[AMPL_NUM_KEYWORDS];
};

SLEQP_RETCODE
sleqp_ampl_keywords_create(SleqpAmplKeywords** star,
                           SleqpOptions* options,
                           SleqpParams* params)
{
  SLEQP_CALL(sleqp_ampl_malloc(star));

  SleqpAmplKeywords* sleqp_keywords = *star;

  SLEQP_CALL(sleqp_options_capture(options));
  sleqp_keywords->options = options;

  SLEQP_CALL(sleqp_params_capture(params));
  sleqp_keywords->params = params;

  SLEQP_CALL(keywords_fill(options,
                           params,
                           sleqp_keywords->keywds,
                           sleqp_keywords->callback_data));

  return SLEQP_OKAY;
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
