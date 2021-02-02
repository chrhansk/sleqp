#cython: language_level=3

from enum import Enum, Flag

cpdef enum Status:
  Optimal    = csleqp.SLEQP_OPTIMAL
  Feasible   = csleqp.SLEQP_FEASIBLE
  Infeasible = csleqp.SLEQP_INFEASIBLE
  Invalid    = csleqp.SLEQP_INVALID

class DerivCheck(Flag):
  Skip             = csleqp.SLEQP_DERIV_CHECK_SKIP
  First            = csleqp.SLEQP_DERIV_CHECK_FIRST
  SecondExhaustive = csleqp.SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE
  SecondSimple     = csleqp.SLEQP_DERIV_CHECK_SECOND_SIMPLE

cpdef enum HessianEval:
  Exact      = csleqp.SLEQP_HESSIAN_EVAL_EXACT
  SR1        = csleqp.SLEQP_HESSIAN_EVAL_SR1
  SimpleBFGS = csleqp.SLEQP_HESSIAN_EVAL_SIMPLE_BFGS
  DampedBFGS = csleqp.SLEQP_HESSIAN_EVAL_DAMPED_BFGS

cpdef enum ValueReason:
  NoReason         = csleqp.SLEQP_VALUE_REASON_NONE
  Init             = csleqp.SLEQP_VALUE_REASON_INIT
  CheckingDeriv    = csleqp.SLEQP_VALUE_REASON_CHECKING_DERIV
  AcceptedIterate  = csleqp.SLEQP_VALUE_REASON_ACCEPTED_ITERATE
  TryingIterate    = csleqp.SLEQP_VALUE_REASON_TRYING_ITERATE
  TryingSOCIterate = csleqp.SLEQP_VALUE_REASON_TRYING_SOC_ITERATE
  RejectedIterate  = csleqp.SLEQP_VALUE_REASON_REJECTED_ITERATE

cpdef enum ActiveState:
  Inactive    = csleqp.SLEQP_INACTIVE
  ActiveLower = csleqp.SLEQP_ACTIVE_LOWER
  ActiveUpper = csleqp.SLEQP_ACTIVE_UPPER
  ActiveBoth  = csleqp.SLEQP_ACTIVE_BOTH

cpdef enum DualEstimationType:
  LP      = csleqp.SLEQP_DUAL_ESTIMATION_TYPE_LP
  LSQ     = csleqp.SLEQP_DUAL_ESTIMATION_TYPE_LSQ

class SolverEvent(Enum):
  AcceptedIterate = csleqp.SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE
