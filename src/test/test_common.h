#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <check.h>

#include "log.h"
#include "types.h"

#define ASSERT_CALL(x)                                          \
        do                                                      \
        {                                                       \
                SLEQP_RETCODE _retcode_ = (x);                  \
                ck_assert_int_eq(_retcode_, SLEQP_OKAY);        \
        }                                                       \
        while(0)

int run_suite(int argc, char *argv[], Suite* suite);

#define TEST_MAIN(suite_factory)                        \
        int main(int argc, char *argv[])                \
        {                                               \
                Suite*suite = suite_factory();          \
                return run_suite(argc, argv, suite);    \
        }

#endif /* TEST_COMMON_H */
