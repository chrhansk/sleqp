#include <check.h>

#include <pthread.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"

#include "test_common.h"

#include "rosenbrock_fixture.h"

#define NUM_THREADS 8

bool start_solving          = false;
pthread_mutex_t start_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t start_cond   = PTHREAD_COND_INITIALIZER;

SLEQP_RETCODE
solve_problem()
{
  SleqpFunc* func;
  SleqpVec* var_lb;
  SleqpVec* var_ub;
  SleqpVec* cons_lb;
  SleqpVec* cons_ub;
  SleqpVec* initial;
  SleqpVec* opt;

  rosenbrock_create(&func,
                    &var_lb,
                    &var_ub,
                    &cons_lb,
                    &cons_ub,
                    &initial,
                    &opt);

  SleqpParams* params;
  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          func,
                                          params,
                                          var_lb,
                                          var_ub,
                                          cons_lb,
                                          cons_ub));

  ASSERT_CALL(
    sleqp_solver_create(&solver, problem, params, options, initial, NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_vec_free(&opt));
  ASSERT_CALL(sleqp_vec_free(&initial));

  ASSERT_CALL(sleqp_vec_free(&cons_ub));
  ASSERT_CALL(sleqp_vec_free(&cons_lb));

  ASSERT_CALL(sleqp_vec_free(&var_ub));
  ASSERT_CALL(sleqp_vec_free(&var_lb));

  ASSERT_CALL(sleqp_func_release(&func));

  return SLEQP_OKAY;
}

void*
solve_thread(void* _unused)
{
  pthread_mutex_lock(&start_mutex);
  while (!start_solving)
  {
    pthread_cond_wait(&start_cond, &start_mutex);
  }
  pthread_mutex_unlock(&start_mutex);

  ASSERT_CALL(solve_problem());

  return NULL;
}

START_TEST(test_thread)
{
  pthread_t threads[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; ++i)
  {
    pthread_create(threads + i, NULL, solve_thread, NULL);
  }

  pthread_mutex_lock(&start_mutex);
  start_solving = true;
  pthread_cond_broadcast(&start_cond);
  pthread_mutex_unlock(&start_mutex);

  for (int i = 0; i < NUM_THREADS; ++i)
  {
    pthread_join(threads[i], NULL);
  }
}
END_TEST

Suite*
thread_test_suite()
{
  Suite* suite;
  TCase* tc_thread;

  suite = suite_create("Thread tests");

  tc_thread = tcase_create("Thread test");

  tcase_add_test(tc_thread, test_thread);
  suite_add_tcase(suite, tc_thread);

  return suite;
}

TEST_MAIN(thread_test_suite)
