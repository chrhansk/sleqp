#ifndef SLEQP_ENUM_H
#define SLEQP_ENUM_H

#include "pub_types.h"

typedef struct
{
  const char* name;
  int value;
} SleqpEnumEntry;

typedef struct
{
  const char* name;
  bool flags;

  SleqpEnumEntry entries[];

} SleqpEnum;

bool
sleqp_enum_member(const SleqpEnum* enum_struct, int value);

#endif /* SLEQP_ENUM_H */
