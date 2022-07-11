#ifndef SLEQP_LOG_H
#define SLEQP_LOG_H

#include "pub_log.h"

/**
 * @file log.h
 * @brief Definition of logging functions.
 **/

#if defined SLEQP_FORMAT_CODES

#define SLEQP_FORMAT_RESET "\x1B[0m"
#define SLEQP_FORMAT_RED "\x1B[31m"
#define SLEQP_FORMAT_GREEN "\x1B[32m"
#define SLEQP_FORMAT_YELLOW "\x1B[33m"
#define SLEQP_FORMAT_BLUE "\x1B[34m"
#define SLEQP_FORMAT_DARK "\x1b[90m"

#define SLEQP_FORMAT_BOLD "\x1B[1m"
#define SLEQP_FORMAT_NO_BOLD "\x1B[22m"

#else

#define SLEQP_FORMAT_RESET ""
#define SLEQP_FORMAT_RED ""
#define SLEQP_FORMAT_GREEN ""
#define SLEQP_FORMAT_YELLOW ""
#define SLEQP_FORMAT_BLUE ""
#define SLEQP_FORMAT_DARK ""

#define SLEQP_FORMAT_BOLD ""
#define SLEQP_FORMAT_NO_BOLD ""

#endif

// mu
#define SLEQP_SYMBOL_MU "\xc2\xb5"
// plus/minus
#define SLEQP_SYMBOL_PM "\xc2\xb1"

#endif /* SLEQP_LOG_H */
