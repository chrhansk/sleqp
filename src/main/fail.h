#ifndef SLEQP_FAIL_H
#define SLEQP_FAIL_H

#include <assert.h>

#include "cmp.h"
#include "defs.h"
#include "error.h"
#include "log.h"
#include "types.h"

#define sleqp_log_assert_fail(assertion)                                       \
  sleqp_raise(SLEQP_FAILED_ASSERTION,                                          \
              "Assertion `%s' failed at %s:%d: %s",                            \
              #assertion,                                                      \
              __FILE__,                                                        \
              __LINE__,                                                        \
              __PRETTY_FUNCTION__);

#define sleqp_log_assert_fail_msg(format, ...)                                 \
  do                                                                           \
  {                                                                            \
    char message_buf[1024];                                                    \
    snprintf(message_buf, 1024, format, ##__VA_ARGS__);                        \
                                                                               \
    sleqp_raise(SLEQP_FAILED_ASSERTION,                                        \
                "Assertion `%s' failed at %s:%d: %s",                          \
                message_buf,                                                   \
                __FILE__,                                                      \
                __LINE__,                                                      \
                __PRETTY_FUNCTION__);                                          \
  } while (false)

#if defined(NDEBUG)

#define SLEQP_DEBUG 0

#define sleqp_assert(...)                                                      \
  do                                                                           \
  {                                                                            \
  } while (false)
#define sleqp_assert_msg(...)                                                  \
  do                                                                           \
  {                                                                            \
  } while (false)

#else

#define SLEQP_DEBUG 1

#define sleqp_assert(expr)                                                     \
  do                                                                           \
  {                                                                            \
    if (expr)                                                                  \
      ;                                                                        \
    else                                                                       \
    {                                                                          \
      sleqp_log_assert_fail(#expr);                                            \
    }                                                                          \
  } while (false)

#define sleqp_assert_msg(expr, format, ...)                                    \
  do                                                                           \
  {                                                                            \
    if (expr)                                                                  \
      ;                                                                        \
    else                                                                       \
    {                                                                          \
      sleqp_log_assert_fail_msg(format, ##__VA_ARGS__);                        \
    }                                                                          \
  } while (false)

#ifdef SLEQP_ENABLE_NUM_ASSERTS

#define SLEQP_HAS_NUM_ASSERT_DEFS

#define sleqp_num_assert(expr) sleqp_assert(expr)

#define sleqp_num_assert_msg(expr, format, ...)                                \
  sleqp_assert_msg(expr, format, ##__VAR_ARGS__)

#define sleqp_assert_is_eq(x, y, eps)                                          \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_eq(x, y, eps))                                               \
    {                                                                          \
      sleqp_log_assert_fail_msg(                                               \
        "%s (= %.14f) == %s (= %.14f) [tolerance: %g]",                        \
        #x,                                                                    \
        x,                                                                     \
        #y,                                                                    \
        y,                                                                     \
        eps);                                                                  \
    }                                                                          \
  } while (false)

#define sleqp_assert_is_lt(x, y, eps)                                          \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_lt(x, y, eps))                                               \
    {                                                                          \
      sleqp_log_assert_fail_msg("%s (= %.14f) < %s (= %.14f) [tolerance: %g]", \
                                #x,                                            \
                                x,                                             \
                                #y,                                            \
                                y,                                             \
                                eps);                                          \
    }                                                                          \
  } while (false)

#define sleqp_assert_is_gt(x, y, eps)                                          \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_gt(x, y, eps))                                               \
    {                                                                          \
      sleqp_log_assert_fail_msg("%s (= %.14f) > %s (= %.14f) [tolerance: %g]", \
                                #x,                                            \
                                x,                                             \
                                #y,                                            \
                                y,                                             \
                                eps);                                          \
    }                                                                          \
  } while (false)

#define sleqp_assert_is_leq(x, y, eps)                                         \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_leq(x, y, eps))                                              \
    {                                                                          \
      sleqp_log_assert_fail_msg(                                               \
        "%s (= %.14f) <= %s (= %.14f) [tolerance: %g]",                        \
        #x,                                                                    \
        x,                                                                     \
        #y,                                                                    \
        y,                                                                     \
        eps);                                                                  \
    }                                                                          \
  } while (false)

#define sleqp_assert_is_geq(x, y, eps)                                         \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_geq(x, y, eps))                                              \
    {                                                                          \
      sleqp_log_assert_fail_msg(                                               \
        "%s (= %.14f) >= %s (= %.14f) [tolerance: %g]",                        \
        #x,                                                                    \
        x,                                                                     \
        #y,                                                                    \
        y,                                                                     \
        eps);                                                                  \
    }                                                                          \
  } while (false)

#define sleqp_assert_is_neg(x, eps)                                            \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_neg(x, eps))                                                 \
    {                                                                          \
      sleqp_log_assert_fail_msg("%s (= %.14f) < 0 [tolerance: %g]",            \
                                #x,                                            \
                                x,                                             \
                                eps);                                          \
    }                                                                          \
  } while (false)

#define sleqp_assert_is_pos(x, eps)                                            \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_pos(x, eps))                                                 \
    {                                                                          \
      sleqp_log_assert_fail_msg("%s (= %.14f) > 0 [tolerance: %g]",            \
                                #x,                                            \
                                x,                                             \
                                eps);                                          \
    }                                                                          \
  } while (false)

#define sleqp_assert_is_zero(x, eps)                                           \
  do                                                                           \
  {                                                                            \
    if (!sleqp_is_zero(x, eps))                                                \
    {                                                                          \
      sleqp_log_assert_fail_msg("%s (= %.14f) = 0 [tolerance: %g]",            \
                                #x,                                            \
                                x,                                             \
                                eps);                                          \
    }                                                                          \
  } while (false)

#endif // SLEQP_ENABLE_NUM_ASSERTS

#endif // NDEBUG

#if defined(SLEQP_HAS_NUM_ASSERT_DEFS)

#define SLEQP_NUM_ASSERT_PARAM(x)

#else

#define sleqp_num_assert(expr)

#define sleqp_num_assert_msg(expr, format, args...)

#define sleqp_assert_is_eq(x, y, eps)

#define sleqp_assert_is_lt(x, y, eps)
#define sleqp_assert_is_gt(x, y, eps)

#define sleqp_assert_is_leq(x, y, eps)
#define sleqp_assert_is_geq(x, y, eps)

#define sleqp_assert_is_neg(x, eps)
#define sleqp_assert_is_pos(x, eps)

#define sleqp_assert_is_zero(x, eps)

#define SLEQP_NUM_ASSERT_PARAM(x) (void)(x)

#endif // SLEQP_HAS_NUM_ASSERT_DEFS

#endif /* SLEQP_FAIL_H */
