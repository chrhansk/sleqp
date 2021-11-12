#ifndef SLEQP_MATH_ERROR_H
#define SLEQP_MATH_ERROR_H

#include <math.h>
#include <fenv.h>
#include "types.h"

#define SLEQP_INIT_MATH_CHECK                           \
        fenv_t fenv_current;                            \
        do                                              \
        {                                               \
                if(math_errhandling & MATH_ERREXCEPT)   \
                {                                       \
                        fegetenv(&fenv_current);        \
                        fesetenv(FE_DFL_ENV);           \
                }                                       \
        }                                               \
        while(false)


#define SLEQP_MATH_CHECK_ERRORS(error_flags, has_errors)                \
        do                                                              \
        {                                                               \
                if(math_errhandling & MATH_ERREXCEPT)                   \
                {                                                       \
                        *(has_errors) = fetestexcept(error_flags);      \
                                                                        \
                        if(*(has_errors))                               \
                        {                                               \
                                sleqp_log_error("Encountered floating point errors (%s, %s, %s, %s, %s)", \
                                                fetestexcept(FE_DIVBYZERO) ? "FE_DIVBYZERO" : "", \
                                                fetestexcept(FE_INEXACT) ? "FE_INEXACT" : "", \
                                                fetestexcept(FE_INVALID) ? "FE_INVALID" : "", \
                                                fetestexcept(FE_OVERFLOW) ? "FE_OVERFLOW" : "", \
                                                fetestexcept(FE_UNDERFLOW) ? "FE_UNDERFLOW" : ""); \
                        }                                               \
                }                                                       \
        }                                                               \
        while(false)

#define SLEQP_MATH_CHECK_WARNINGS(warn_flags)                           \
        do                                                              \
        {                                                               \
                if(math_errhandling & MATH_ERREXCEPT)                   \
                {                                                       \
                        const bool has_warning = fetestexcept(warn_flags); \
                                                                        \
                        if(has_warning)                                 \
                        {                                               \
                                sleqp_log_warn("Encountered floating point errors (%s, %s, %s, %s, %s)", \
                                               fetestexcept(FE_DIVBYZERO) ? "FE_DIVBYZERO" : "", \
                                               fetestexcept(FE_INEXACT) ? "FE_INEXACT" : "", \
                                               fetestexcept(FE_INVALID) ? "FE_INVALID" : "", \
                                               fetestexcept(FE_OVERFLOW) ? "FE_OVERFLOW" : "", \
                                               fetestexcept(FE_UNDERFLOW) ? "FE_UNDERFLOW" : ""); \
                        }                                               \
                }                                                       \
        }                                                               \
        while(false)

#define SLEQP_MATH_CHECK(error_flags, warn_flags)                       \
        do                                                              \
        {                                                               \
                if(math_errhandling & MATH_ERREXCEPT)                   \
                {                                                       \
                        bool has_errors = false;                        \
                        SLEQP_MATH_CHECK_ERRORS(error_flags, &has_errors); \
                        if(!has_errors)                                 \
                        {                                               \
                                SLEQP_MATH_CHECK_WARNINGS(warn_flags);  \
                        }                                               \
                        fesetenv(&fenv_current);                        \
                        if(has_errors)                                  \
                        {                                               \
                                return SLEQP_MATH_ERROR;                \
                        }                                               \
                }                                                       \
        }                                                               \
        while(false)


#endif /* SLEQP_MATH_ERROR_H */
