#ifndef SLEQP_H
#define SLEQP_H

/**
 * @file sleqp.h
 * @brief SLEQP callable libary.
 **/

#include "sleqp_aug_jacobian.h"
#include "sleqp_assert.h"
#include "sleqp_bfgs.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_defs.h"
#include "sleqp_deriv_check.h"
#include "sleqp_dual_estimation.h"
#include "sleqp_feas.h"
#include "sleqp_func.h"
#include "sleqp_hess_struct.h"
#include "sleqp_iterate.h"
#include "sleqp_linesearch.h"
#include "sleqp_log.h"
#include "sleqp_lsq.h"
#include "sleqp_options.h"
#include "sleqp_mem.h"
#include "sleqp_merit.h"
#include "sleqp_newton.h"
#include "sleqp_params.h"
#include "sleqp_problem.h"
#include "sleqp_problem_scaling.h"
#include "sleqp_scale.h"
#include "sleqp_soc.h"
#include "sleqp_solver.h"
#include "sleqp_sr1.h"
#include "sleqp_types.h"
#include "sleqp_working_set.h"
#include "sleqp_util.h"

#include "tr/sleqp_tr_solver.h"

#include "lp/sleqp_lpi.h"
#include "lp/sleqp_lpi_types.h"

#include "sparse/sleqp_sparse_factorization.h"
#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#endif /* SLEQP_H */
