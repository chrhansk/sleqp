#include "enum.h"

#include "cmp.h"

#include <string.h>
#include <threads.h>

static bool
enum_member(const SleqpEnum* enum_struct, int value)
{
  const SleqpEnumEntry* entry = enum_struct->entries;

  for (; entry->name; ++entry)
  {
    if (entry->value == value)
    {
      return true;
    }
  }

  return false;
}

static bool
flags_member(const SleqpEnum* enum_struct, int value)
{
  const SleqpEnumEntry* entry = enum_struct->entries;

  int actual_value = 0;

  for (; entry->name; ++entry)
  {
    if (value & entry->value)
    {
      actual_value |= entry->value;
    }
  }

  return value == actual_value;
}

bool
sleqp_enum_member(const SleqpEnum* enum_struct, int value)
{
  if (enum_struct->flags)
  {
    return flags_member(enum_struct, value);
  }

  return enum_member(enum_struct, value);
}
