#ifndef SLEQP_ASSERT_H
#define SLEQP_ASSERT_H

#include "sleqp_defs.h"
#include "sleqp_types.h"
#include "sleqp_log.h"

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

#define sleqp_assert()
#define sleqp_assert_msg()

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

#define sleqp_assert_msg(expr, format, args...)       \
  do                                                  \
  {                                                   \
    if(expr);                                         \
    else                                              \
    {                                                 \
      sleqp_log_assert_fail_msg(__FILE__,             \
                                __LINE__,             \
                                __PRETTY_FUNCTION__,  \
                                format,               \
                                args);                \
      return SLEQP_FAILED_ASSERTION;                  \
    }                                                 \
  }                                                   \
  while(false)

#ifdef SLEQP_ENABLE_NUM_ASSERTS

#define SLEQP_HAS_NUM_ASSERT_DEFS

#define sleqp_num_assert(expr)                  \
  sleqp_assert(expr)

#define sleqp_num_assert_msg(expr, format, args...) \
  sleqp_assert_msg(expr, format, args)

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

#if !defined(SLEQP_HAS_NUM_ASSERT_DEFS)

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


#endif // SLEQP_HAS_NUM_ASSERT_DEFS

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_ASSERT_H */
