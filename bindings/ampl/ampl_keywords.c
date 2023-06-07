#include "ampl_keywords.h"
#include "sleqp/pub_settings.h"

#include <string.h>

typedef struct
{
  union
  {
    SleqpSettings* settings;
    SleqpAmplKeywords* keywords;
  } data;

  int index;

} CallbackData;

enum
{
  ITER_LIMIT_MAXITER = 0,
  ITER_LIMIT_MAXIT,
  TIME_LIMIT,
  LOG_LEVEL_PRINT_LEVEL,
  LOG_LEVEL_OUTLEV,
  WANTSOL,
  HALTONERROR,
  NUM_EXTRA
};

typedef enum
{
  POS_ENUM          = 0,
  POS_INT           = SLEQP_NUM_ENUM_SETTINGS,
  POS_BOOL          = POS_INT + SLEQP_NUM_INT_SETTINGS,
  POS_PAR           = POS_BOOL + SLEQP_NUM_BOOL_SETTINGS,
  POS_EXTRA         = POS_PAR + SLEQP_NUM_REAL_SETTINGS,
  AMPL_NUM_KEYWORDS = POS_EXTRA + NUM_EXTRA
} AMPL_KEYWORDS;

struct SleqpAmplKeywords
{
  SleqpSettings* settings;

  CallbackData callback_data[AMPL_NUM_KEYWORDS];
  keyword keywds[AMPL_NUM_KEYWORDS];

  double time_limit;
  int iteration_limit;
  bool halt_on_error;
};

static char*
kwdfunc_enum(Option_Info* oi, keyword* kw, char* value)
{
  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  SLEQP_RETCODE retcode
    = sleqp_settings_set_enum_value(callback_data->data.settings,
                                    callback_data->index,
                                    int_val);

  if (retcode != SLEQP_OKAY)
  {
    fprintf(stderr,
            "Invalid value \"%d\" for keyword \"%s\"\n",
            int_val,
            kw->name);
    badopt_ASL(oi);
  }

  return retval;
}

static char*
kwdfunc_int(Option_Info* oi, keyword* kw, char* value)
{
  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  SLEQP_RETCODE retcode
    = sleqp_settings_set_int_value(callback_data->data.settings,
                                  callback_data->index,
                                  int_val);

  if (retcode != SLEQP_OKAY)
  {
    fprintf(stderr,
            "Invalid value \"%d\" for keyword \"%s\"\n",
            int_val,
            kw->name);
    badopt_ASL(oi);
  }

  return retval;
}

static char*
kwdfunc_bool(Option_Info* oi, keyword* kw, char* value)
{
  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  if (!(int_val == 0 || int_val == 1))
  {
    fprintf(stderr,
            "Invalid value \"%d\" for keyword \"%s\"\n",
            int_val,
            kw->name);
    badopt_ASL(oi);
  }

  SLEQP_RETCODE retcode
    = sleqp_settings_set_bool_value(callback_data->data.settings,
                                   callback_data->index,
                                   !!(int_val));

  if (retcode != SLEQP_OKAY)
  {
    badopt_ASL(oi);
  }

  return retval;
}

static char*
kwdfunc_param(Option_Info* oi, keyword* kw, char* value)
{
  CallbackData* callback_data = (CallbackData*)kw->info;

  double real_val;
  kw->info = &real_val;

  char* retval = D_val(oi, kw, value);

  SLEQP_RETCODE retcode = sleqp_settings_set_real_value(callback_data->data.settings,
                                                        callback_data->index,
                                                        real_val);

  if (retcode != SLEQP_OKAY)
  {
    fprintf(stderr,
            "Invalid value \"%f\" for keyword \"%s\"\n",
            real_val,
            kw->name);
    badopt_ASL(oi);
  }

  return retval;
}

static char*
kwdfunc_iterlimit(Option_Info* oi, keyword* kw, char* value)
{
  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  if (int_val != SLEQP_NONE)
  {
    if (int_val < 0)
    {
      fprintf(stderr,
              "Invalid value \"%d\" for keyword \"%s\"\n",
              int_val,
              kw->name);
      badopt_ASL(oi);
    }

    callback_data->data.keywords->iteration_limit = int_val;
  }

  return retval;
}

static char*
kwdfunc_timelimit(Option_Info* oi, keyword* kw, char* value)
{
  CallbackData* callback_data = (CallbackData*)kw->info;

  double real_val;
  kw->info = &real_val;

  char* retval = D_val(oi, kw, value);

  if (real_val != SLEQP_NONE)
  {
    if (real_val < 0)
    {
      fprintf(stderr,
              "Invalid value \"%f\" for keyword \"%s\"\n",
              real_val,
              kw->name);
      badopt_ASL(oi);
    }

    callback_data->data.keywords->time_limit = real_val;
  }

  return retval;
}

static char*
kwdfunc_log_level(Option_Info* oi, keyword* kw, char* value)
{
  int int_val;
  kw->info = &int_val;

  char* retval = I_val(oi, kw, value);

  if (int_val != SLEQP_NONE)
  {
    if (int_val < SLEQP_LOG_SILENT || int_val > SLEQP_LOG_DEBUG)
    {
      fprintf(stderr,
              "Invalid value \"%d\" for keyword \"%s\"\n",
              int_val,
              kw->name);
      badopt_ASL(oi);
    }

    sleqp_log_set_level(int_val);
  }

  return retval;
}

static char*
kwdfunc_haltonerror(Option_Info* oi, keyword* kw, char* value)
{
  CallbackData* callback_data = (CallbackData*)kw->info;

  int int_val;
  kw->info = &int_val;

  char* str_val;
  kw->info     = &str_val;
  char* retval = C_val(oi, kw, value);

  if (strcmp(str_val, "yes") == 0)
  {
    callback_data->data.keywords->halt_on_error = true;
  }
  else if (strcmp(str_val, "no") == 0)
  {
    callback_data->data.keywords->halt_on_error = false;
  }
  else
  {
    fprintf(stderr,
            "Invalid value \"%s\" for keyword \"%s\"\n",
            str_val,
            kw->name);
    badopt_ASL(oi);
  }

  return retval;
}

static int
compare_kwds(const void* first, const void* second)
{
  return strcmp(((keyword*)first)->name, ((keyword*)second)->name);
}

static SLEQP_RETCODE
keywords_fill(SleqpAmplKeywords* sleqp_keywords,
              SleqpSettings* settings)
{
  CallbackData* callback_data = sleqp_keywords->callback_data;

  int pos = POS_ENUM;

  keyword* kwds = sleqp_keywords->keywds;

  for (; pos < POS_INT; ++pos)
  {
    SLEQP_SETTINGS_ENUM option_value = pos;

    callback_data[pos]
      = (CallbackData){.data.settings = settings, .index = option_value};

    kwds[pos]
      = (keyword){.name = strdup(sleqp_settings_enum_name(option_value)),
                  .kf   = kwdfunc_enum,
                  .info = callback_data + pos,
                  .desc = strdup(sleqp_settings_enum_desc(option_value))};
  }

  for (; pos < POS_BOOL; ++pos)
  {
    SLEQP_SETTINGS_INT option_value = pos - POS_INT;

    callback_data[pos]
      = (CallbackData){.data.settings = settings, .index = option_value};

    kwds[pos] = (keyword){.name = strdup(sleqp_settings_int_name(option_value)),
                          .kf   = kwdfunc_int,
                          .info = callback_data + pos,
                          .desc = strdup(sleqp_settings_int_desc(option_value))};
  }

  for (; pos < POS_PAR; ++pos)
  {
    SLEQP_SETTINGS_BOOL option_value = pos - POS_BOOL;

    callback_data[pos]
      = (CallbackData){.data.settings = settings, .index = option_value};

    kwds[pos]
      = (keyword){.name = strdup(sleqp_settings_bool_name(option_value)),
                  .kf   = kwdfunc_bool,
                  .info = callback_data + pos,
                  .desc = strdup(sleqp_settings_bool_desc(option_value))};
  }

  for (; pos < POS_EXTRA; ++pos)
  {
    SLEQP_SETTINGS_REAL real_value = pos - POS_PAR;

    callback_data[pos]
      = (CallbackData){.data.settings = settings, .index = real_value};

    kwds[pos] = (keyword){.name = strdup(sleqp_settings_real_name(real_value)),
                          .kf   = kwdfunc_param,
                          .info = callback_data + pos,
                          .desc = strdup(sleqp_settings_real_desc(real_value))};
  }

  for (pos = POS_EXTRA; pos < AMPL_NUM_KEYWORDS; ++pos)
  {
    callback_data[pos]
      = (CallbackData){.data.keywords = sleqp_keywords, .index = 0};
  }

  {
    pos = POS_EXTRA + ITER_LIMIT_MAXITER;

    kwds[pos] = (keyword){.name = strdup("max_iter"),
                          .kf   = kwdfunc_iterlimit,
                          .info = callback_data + pos,
                          .desc = strdup("Maximum number of iterations")};
  }

  {
    pos = POS_EXTRA + ITER_LIMIT_MAXIT;

    kwds[pos] = (keyword){.name = strdup("maxit"),
                          .kf   = kwdfunc_iterlimit,
                          .info = callback_data + pos,
                          .desc = strdup("Alias for 'max_iter'")};
  }

  {
    pos = POS_EXTRA + TIME_LIMIT;

    kwds[pos] = (keyword){.name = strdup("max_wall_time"),
                          .kf   = kwdfunc_timelimit,
                          .info = callback_data + pos,
                          .desc = strdup("Wallclock time limit")};
  }

  {
    pos = POS_EXTRA + LOG_LEVEL_PRINT_LEVEL;

    kwds[pos] = (keyword){.name = strdup("print_level"),
                          .kf   = kwdfunc_log_level,
                          .info = NULL,
                          .desc = strdup("Verbosity level")};
  }

  {
    pos = POS_EXTRA + LOG_LEVEL_OUTLEV;

    kwds[pos] = (keyword){.name = strdup("outlev"),
                          .kf   = kwdfunc_log_level,
                          .info = NULL,
                          .desc = strdup("Alias for 'print_level'")};
  }

  {
    pos = POS_EXTRA + WANTSOL;

    kwds[pos] = (keyword){
      .name = strdup("wantsol"),
      .kf   = WS_val,
      .info = NULL,
      .desc = strdup(
        "solution report without -AMPL: sum of 1 (write .sol file), 2 "
        "(print primal variable values), 4 (print dual variable values), "
        "8 (do not print solution message)")};
  }

  {
    pos = POS_EXTRA + HALTONERROR;

    kwds[pos]
      = (keyword){.name = strdup("halt_on_ampl_error"),
                  .kf   = kwdfunc_haltonerror,
                  .info = callback_data + pos,
                  .desc = strdup("Exit with message on evaluation error")};
  }

  // Keywords must be sorted alphabetically
  qsort(kwds, AMPL_NUM_KEYWORDS, sizeof(keyword), compare_kwds);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_keywords_create(SleqpAmplKeywords** star,
                           SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpAmplKeywords* ampl_keywords = *star;

  (*ampl_keywords) = (SleqpAmplKeywords){0};

  ampl_keywords->time_limit      = SLEQP_NONE;
  ampl_keywords->iteration_limit = SLEQP_NONE;
  ampl_keywords->halt_on_error   = false;

  SLEQP_CALL(sleqp_settings_capture(settings));
  ampl_keywords->settings = settings;

  SLEQP_CALL(keywords_fill(ampl_keywords, settings));

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

bool
sleqp_ampl_keywords_halt_on_error(SleqpAmplKeywords* ampl_keywords)
{
  return ampl_keywords->halt_on_error;
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

  SLEQP_CALL(sleqp_settings_release(&sleqp_keywords->settings));

  keyword* kwds = sleqp_keywords->keywds;

  for (int pos = POS_ENUM; pos < AMPL_NUM_KEYWORDS; ++pos)
  {
    sleqp_free(&(kwds[pos].name));
    sleqp_free(&(kwds[pos].desc));
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}
