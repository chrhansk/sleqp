#cython: language_level=3

cpdef enum Status:
  Optimal    = csleqp.SLEQP_OPTIMAL
  Feasible   = csleqp.SLEQP_FEASIBLE
  Infeasible = csleqp.SLEQP_INFEASIBLE
  Invalid    = csleqp.SLEQP_INVALID

cpdef enum DerivCheck:
  Skip  = csleqp.SLEQP_DERIV_CHECK_SKIP
  First = csleqp.SLEQP_DERIV_CHECK_FIRST
  Sec   = csleqp.SLEQP_DERIV_CHECK_SEC
  Both  = csleqp.SLEQP_DERIV_CHECK_BOTH

cpdef enum HessianEval:
  Exact      = csleqp.SLEQP_HESSIAN_EVAL_EXACT
  SR1        = csleqp.SLEQP_HESSIAN_EVAL_SR1
  SimpleBFGS = csleqp.SLEQP_HESSIAN_EVAL_SIMPLE_BFGS
  DampedBFGS = csleqp.SLEQP_HESSIAN_EVAL_DAMPED_BFGS
