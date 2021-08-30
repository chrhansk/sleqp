#ifndef SLEQP_FAIL_H
#define SLEQP_FAIL_H

#include "cmp.h"
#include "defs.h"
#include "types.h"
#include "log.h"

#ifdef __cplusplus
extern "C" {
#endif

  void sleqp_log_assert_fail(const char *assertion,
                             const char *file,
                             unsigned int line,
                             const char *function);

  void sleqp_log_assert_fail_msg(const char *file,
                                 unsigned int line,
                                 const char *function,
                                 const char* format,
                                 ...);

#if defined(NDEBUG)

#define sleqp_assert(...) do{} while(false)
#define sleqp_assert_msg(...) do{} while(false)

#else
#define sleqp_assert(expr)                        \
  do                                              \
  {                                               \
    if(expr);                                     \
    else                                          \
    {                                             \
      sleqp_log_assert_fail(#expr,                \
                            __FILE__,             \
                            __LINE__,             \
                            __PRETTY_FUNCTION__); \
      return SLEQP_FAILED_ASSERTION;              \
    }                                             \
  }                                               \
  while(false)

#define sleqp_assert_msg(expr, format, ...)           \
  do                                                  \
  {                                                   \
    if(expr);                                         \
    else                                              \
    {                                                 \
      sleqp_log_assert_fail_msg(__FILE__,             \
                                __LINE__,             \
                                __PRETTY_FUNCTION__,  \
                                format,               \
                                ##__VA_ARGS__);       \
      return SLEQP_FAILED_ASSERTION;                  \
    }                                                 \
  }                                                   \
  while(false)

#ifdef SLEQP_ENABLE_NUM_ASSERTS

#define SLEQP_HAS_NUM_ASSERT_DEFS

#define sleqp_num_assert(expr)                  \
  sleqp_assert(expr)

#define sleqp_num_assert_msg(expr, format, ...) \
  sleqp_assert_msg(expr, format, ##__VAR_ARGS__)

#define sleqp_assert_is_eq(x, y, eps)                                              \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_eq(x, y, eps))                                                    \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) == %s (= %.14f) [tolerance: %.14f]", \
                                #x,                                                \
                                x,                                                 \
                                #y,                                                \
                                y,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)

#define sleqp_assert_is_lt(x, y, eps)                                              \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_lt(x, y, eps))                                                    \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) < %s (= %.14f) [tolerance: %.14f]",  \
                                #x,                                                \
                                x,                                                 \
                                #y,                                                \
                                y,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)

#define sleqp_assert_is_gt(x, y, eps)                                              \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_gt(x, y, eps))                                                    \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) > %s (= %.14f) [tolerance: %.14f]",  \
                                #x,                                                \
                                x,                                                 \
                                #y,                                                \
                                y,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)


#define sleqp_assert_is_leq(x, y, eps)                                             \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_leq(x, y, eps))                                                   \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) <= %s (= %.14f) [tolerance: %.14f]", \
                                #x,                                                \
                                x,                                                 \
                                #y,                                                \
                                y,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)

#define sleqp_assert_is_geq(x, y, eps)                                             \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_geq(x, y, eps))                                                   \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) >= %s (= %.14f) [tolerance: %.14f]", \
                                #x,                                                \
                                x,                                                 \
                                #y,                                                \
                                y,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)


#define sleqp_assert_is_neg(x, eps)                                                \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_neg(x, eps))                                                      \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) < 0 [tolerance: %.14f]",             \
                                #x,                                                \
                                x,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)


#define sleqp_assert_is_pos(x, eps)                                                \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_pos(x, eps))                                                      \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) > 0 [tolerance: %.14f]",             \
                                #x,                                                \
                                x,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)

#define sleqp_assert_is_zero(x, eps)                                               \
  do                                                                               \
  {                                                                                \
    if(!sleqp_is_zero(x, eps))                                                     \
    {                                                                              \
      sleqp_log_assert_fail_msg(__FILE__,                                          \
                                __LINE__,                                          \
                                __PRETTY_FUNCTION__,                               \
                                "%s (= %.14f) = 0 [tolerance: %.14f]",             \
                                #x,                                                \
                                x,                                                 \
                                eps);                                              \
    }                                                                              \
  } while(false)

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

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FAIL_H */
