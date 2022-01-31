#include <stdio.h>
#include <stdlib.h>

#include "enum.h"
#include "options.h"
#include "types.h"

#define BUF_SIZE 256

void
generate_enum_entries(FILE* file, const SleqpEnumEntry* entry)
{
  for (; entry->name; ++entry)
  {
    fprintf(file, "    function [out] = %s()\n", entry->name);

    fprintf(file, "      out = %d;\n", entry->value);

    fprintf(file, "    end\n");
  }
}

void
generate_mex_enum(const SleqpEnum* sleqp_enum, const char* desc)
{
  char buffer[BUF_SIZE];

  snprintf(buffer, BUF_SIZE, "%s.m", sleqp_enum->name);

  FILE* file = fopen(buffer, "w");

  fprintf(file, "%% Auto-generated file\n\n");

  fprintf(file, "classdef %s\n", sleqp_enum->name);

  fprintf(file, "  %% %s %s\n", sleqp_enum->name, desc);

  fprintf(file, "  methods (Static = true)\n");

  generate_enum_entries(file, sleqp_enum->entries);

  fprintf(file, "  end\n");

  fprintf(file, "end\n");

  fclose(file);
}

int
main()
{
  SleqpOptions* options;

  SLEQP_CALL(sleqp_options_create(&options));

  {
    const SleqpEnum* option_enums[SLEQP_NUM_ENUM_OPTIONS]
      = {sleqp_enum_deriv_check(),
         sleqp_enum_hess_eval(),
         sleqp_enum_dual_estimation(),
         NULL, // Floating point warning flags
         NULL, // Floating point error flags
         sleqp_enum_bfgs_sizing(),
         sleqp_enum_tr_solver(),
         sleqp_enum_polishing_type(),
         sleqp_enum_step_rule(),
         sleqp_enum_linesearch(),
         sleqp_enum_parametric_cauchy(),
         sleqp_enum_initial_tr()};

    for (int i = 0; i < SLEQP_NUM_ENUM_OPTIONS; ++i)
    {
      if (!option_enums[i])
      {
        continue;
      }

      generate_mex_enum(option_enums[i],
                        sleqp_options_enum_desc((SLEQP_OPTION_ENUM)i));
    }

    generate_mex_enum(sleqp_enum_active_state(),
                      "State of a variable or constrain in working set");

    generate_mex_enum(sleqp_enum_status(),
                      "Status of the solver after optimization");

    generate_mex_enum(sleqp_enum_steptype(),
                      "Type of steps taken during the solution process");
  }

  return EXIT_SUCCESS;
}
