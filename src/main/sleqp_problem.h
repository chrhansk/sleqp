#ifndef SLEQP_PROBLEM_H
#define SLEQP_PROBLEM_H

#include "sleqp_func.h"
#include "sleqp_types.h"
#include "sleqp_params.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * @defgroup problem Problem definition
   * @{
   **/

  /**
   * A NLP is given as
   *
   * \f[
   * \begin{aligned}
   * \min \: & f(x)                       \\
   * \text{s.t. } \: & l \leq c(x) \leq u \\
   * & x_l \leq x \leq x_u
   * \end{aligned}
   * \f]
   *
   * where \f$ f : \mathbb{R}^{n} \to \mathbb{R} \f$, \f$ c : \mathbb{R}^{n} \to \mathbb{R}^{m} \f$
   * are functions, \f$ l, u \in \mathbb{R}^{m}, l \leq u \f$ are the constraint bounds, and
   * \f$ x_l, x_u \in \mathbb{R}^{n}, x_l \leq x_u \f$ are the variable bounds.
   *
   * @see Functions
   *
   **/
  typedef struct SleqpProblem
  {
    /**
     * The functions \f$ f\f$, and \f$ c \f$ combined.
     **/
    SleqpFunc* func;

    /**
     * lower variable bounds \f$ x_l \f$.
     **/
    SleqpSparseVec* var_lb;

    /**
     * upper variable bounds \f$ x_u \f$.
     **/
    SleqpSparseVec* var_ub;

    /**
     * lower constraint bounds \f$ l \f$.
     **/
    SleqpSparseVec* cons_lb;

    /**
     * upper constraint bounds \f$ u \f$.
     **/
    SleqpSparseVec* cons_ub;

    /**
     * number of variables \f$ n \f$.
     **/
    int num_variables;

    /**
     * number of constraints \f$ m \f$.
     **/
    int num_constraints;

  } SleqpProblem;


  /**
   * Creates a new problem, which will be saved
   * in the newly allocated pointer.
   *
   * The function pointer is kept as a reference,
   * The upper / lower bound vectors are copied.
   *
   **/
  SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                     SleqpFunc* func,
                                     SleqpParams* params,
                                     SleqpSparseVec* var_lb,
                                     SleqpSparseVec* var_ub,
                                     SleqpSparseVec* cons_lb,
                                     SleqpSparseVec* cons_ub);

  /**
   * Frees a previously created problem.
   **/
  SLEQP_RETCODE sleqp_problem_free(SleqpProblem** star);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PROBLEM_H */
