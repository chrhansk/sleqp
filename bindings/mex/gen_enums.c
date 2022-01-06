#include <stdio.h>

#include "enum.h"
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
generate_mex_enum(const SleqpEnum* sleqp_enum)
{
  char buffer[BUF_SIZE];

  snprintf(buffer, BUF_SIZE, "%s.m", sleqp_enum->name);

  FILE* file = fopen(buffer, "w");

  fprintf(file, "%% Auto-generated file\n\n");

  fprintf(file, "classdef %s\n", sleqp_enum->name);

  fprintf(file, "  methods (Static = true)\n");

  generate_enum_entries(file, sleqp_enum->entries);

  fprintf(file, "  end\n");

  fprintf(file, "end\n");

  fclose(file);
}

int
main(int argc, char* argv[])
{
  {
    const SleqpEnum* sleqp_enums[] = {sleqp_enum_active_state(),
                                      sleqp_enum_deriv_check(),
                                      sleqp_enum_hess_eval(),
                                      sleqp_enum_bfgs_sizing(),
                                      sleqp_enum_steptype(),
                                      sleqp_enum_dual_estimation(),
                                      sleqp_enum_tr_solver(),
                                      sleqp_enum_polishing_type(),
                                      sleqp_enum_parametric_cauchy(),
                                      sleqp_enum_initial_tr(),
                                      sleqp_enum_linesearch(),
                                      sleqp_enum_status(),
                                      sleqp_enum_step_rule()};

    const int num_enums = sizeof(sleqp_enums) / sizeof(SleqpEnum*);

    for (int i = 0; i < num_enums; ++i)
    {
      generate_mex_enum(sleqp_enums[i]);
    }
  }

  return 0;
}
