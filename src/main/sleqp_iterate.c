#include "sleqp_iterate.h"

#include "sleqp_mem.h"

SLEQP_RETCODE sleqp_iterate_create(SleqpIterate** star,
                                   SleqpProblem* problem,
                                   SleqpSparseVec* x)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpIterate* iterate = *star;

  size_t num_variables = problem->num_variables;
  size_t num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->func_grad,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->cons_val,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_matrix_create(&iterate->cons_jac,
                                        num_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->cons_dual,
                                        num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->vars_dual,
                                        num_variables,
                                        0));

  iterate->x = x;

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_iterate_free(SleqpIterate** star)
{
  SleqpIterate* iterate = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->vars_dual));
  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->cons_dual));

  SLEQP_CALL(sleqp_sparse_matrix_free(&iterate->cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->func_grad));

  sleqp_free(star);

  return SLEQP_OKAY;
}
