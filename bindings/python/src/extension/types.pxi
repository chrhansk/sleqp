#cython: language_level=3

from enum import Enum, Flag, auto


class _DocEnum(Enum):
    def __new__(cls, value, doc=None):
        self = object.__new__(cls)  # calling super().__new__(value) here would fail
        self._value_ = value
        if doc is not None:
            self.__doc__ = doc
            self.desc    = doc
        return self


class Status(_DocEnum):
  """
  The status resulting from the solution process
  """
  Unkown         = csleqp.SLEQP_STATUS_UNKNOWN, "Unknown state"
  Running        = csleqp.SLEQP_STATUS_RUNNING, "Solver is running"
  Optimal        = csleqp.SLEQP_STATUS_OPTIMAL, "An optimal solution was found"
  Infeasible     = csleqp.SLEQP_STATUS_INFEASIBLE, "Problem was detected to be infeasible"
  Unbounded      = csleqp.SLEQP_STATUS_UNBOUNDED, "Problem appears unbounded"
  AbortDeadpoint = csleqp.SLEQP_STATUS_ABORT_DEADPOINT, "Reached a dead point"
  AbortIter      = csleqp.SLEQP_STATUS_ABORT_ITER, "Aborted after reaching iteration limit"
  AbortManual    = csleqp.SLEQP_STATUS_ABORT_MANUAL, "Manual abort requested"
  AbortTime      = csleqp.SLEQP_STATUS_ABORT_TIME, "Aborted after reaching time limit"


class DerivCheck(Flag):
  Skip             = csleqp.SLEQP_DERIV_CHECK_SKIP
  FirstFunc        = csleqp.SLEQP_DERIV_CHECK_FIRST_OBJ
  FirstCons        = csleqp.SLEQP_DERIV_CHECK_FIRST_CONS
  First            = csleqp.SLEQP_DERIV_CHECK_FIRST
  SecondFunc       = csleqp.SLEQP_DERIV_CHECK_SECOND_OBJ
  SecondCons       = csleqp.SLEQP_DERIV_CHECK_SECOND_CONS
  SecondExhaustive = csleqp.SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE
  SecondSimple     = csleqp.SLEQP_DERIV_CHECK_SECOND_SIMPLE


class HessianEval(_DocEnum):
  """
  The evaluation method used for Hessian products
  """
  Exact      = csleqp.SLEQP_HESS_EVAL_EXACT, "Exact evaluation"
  SR1        = csleqp.SLEQP_HESS_EVAL_SR1, "The SR1 quasi-Newton method"
  SimpleBFGS = csleqp.SLEQP_HESS_EVAL_SIMPLE_BFGS, "The BFGS method for convex functions"
  DampedBFGS = csleqp.SLEQP_HESS_EVAL_DAMPED_BFGS, "A damped BFGS method "


class Sizing(Enum):
  NoSizing   = csleqp.SLEQP_BFGS_SIZING_NONE
  CenteredOL = csleqp.SLEQP_BFGS_SIZING_CENTERED_OL


class TRSolver(_DocEnum):
  """
  The trust-region solver used
  """
  TRlib = csleqp.SLEQP_TR_SOLVER_TRLIB, "The trlib implementation of the Generalized Lanczos method"
  CG    = csleqp.SLEQP_TR_SOLVER_CG, "Steihaug's conjugate gradient method"
  LSQR  = csleqp.SLEQP_TR_SOLVER_LSQR, "LSQR solver for LSQ functions"
  Auto  = csleqp.SLEQP_TR_SOLVER_AUTO, "Automatically chosen"


class PolishingType(_DocEnum):
  """
  The polishing methods use
  """
  NoPolishing = csleqp.SLEQP_POLISHING_NONE, "No polishing"
  ZeroDual    = csleqp.SLEQP_POLISHING_ZERO_DUAL, "Remove bounds and constraints with zero duals"
  Inactive    = csleqp.SLEQP_POLISHING_INACTIVE, "Remove inactive bounds and constraints"


class StepRule(_DocEnum):
  """
  The step rule used
  """
  Direct = csleqp.SLEQP_STEP_RULE_DIRECT, "Direct"
  Window = csleqp.SLEQP_STEP_RULE_WINDOW, "Sliding window"
  MinStep = csleqp.SLEQP_STEP_RULE_MINSTEP, "Min-step"


class LineSearch(_DocEnum):
  """
  The linesearch used
  """
  Exact  = csleqp.SLEQP_LINESEARCH_EXACT, "Exact line search"
  Approx = csleqp.SLEQP_LINESEARCH_APPROX, "Approximate line search"


class ParametricCauchy(_DocEnum):
  """
  The parametric Cauchy type to be used
  """
  Disabled = csleqp.SLEQP_PARAMETRIC_CAUCHY_DISABLED, "Disable parametric Cauchy"
  Coarse   = csleqp.SLEQP_PARAMETRIC_CAUCHY_COARSE, "Coarse parametric Cauchy"
  Fine     = csleqp.SLEQP_PARAMETRIC_CAUCHY_FINE, "Fine parametric Cauchy"


class InitialTRChoice(_DocEnum):
  """
  Choice for initial trust radius
  """
  Narrow = csleqp.SLEQP_INITIAL_TR_CHOICE_NARROW, "Narrow"
  Wide   = csleqp.SLEQP_INITIAL_TR_CHOICE_WIDE, "Wide"

class ValueReason(_DocEnum):
  """
  The reason for setting a new function value
  """
  NoReason         = csleqp.SLEQP_VALUE_REASON_NONE, "No reason given"
  Init             = csleqp.SLEQP_VALUE_REASON_INIT, "Initialization"
  CheckingDeriv    = csleqp.SLEQP_VALUE_REASON_CHECKING_DERIV, "Derivative check"
  AcceptedIterate  = csleqp.SLEQP_VALUE_REASON_ACCEPTED_ITERATE, "Accepted a trial iterate"
  TryingIterate    = csleqp.SLEQP_VALUE_REASON_TRYING_ITERATE, "Trying out a trial iterate"
  TryingSOCIterate = csleqp.SLEQP_VALUE_REASON_TRYING_SOC_ITERATE, "Trying a second-order correction"
  RejectedIterate  = csleqp.SLEQP_VALUE_REASON_REJECTED_ITERATE, "Rejected an iterate restoring the previous"


class ActiveState(_DocEnum):
  """
  The state of a variable or constraint in the working set
  """
  Inactive    = csleqp.SLEQP_INACTIVE, "Variable or constraint is inactive"
  ActiveLower = csleqp.SLEQP_ACTIVE_LOWER, "Variable or constraint is active at its lower bound"
  ActiveUpper = csleqp.SLEQP_ACTIVE_UPPER, "Variable or constraint is active at its upper bound"
  ActiveBoth  = csleqp.SLEQP_ACTIVE_BOTH, "Variable or constraint is active at both bounds (which must coincide)"


class DualEstimationType(_DocEnum):
  """
  The type of dual estimation used
  """
  LP      = csleqp.SLEQP_DUAL_ESTIMATION_TYPE_LP, "Dual variables of the LP"
  LSQ     = csleqp.SLEQP_DUAL_ESTIMATION_TYPE_LSQ, "Least-squares estimation"


class SolverEvent(Enum):
  AcceptedIterate    = csleqp.SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE
  PerformedIteration = csleqp.SLEQP_SOLVER_EVENT_PERFORMED_ITERATION
  Finished           = csleqp.SLEQP_SOLVER_EVENT_FINISHED


class SolverState(Enum):
  TrustRadius              = auto()
  LPTrustRadius            = auto()
  ScaledFuncVal            = auto()
  ScaledMeritVal           = auto()
  ScaledFeasRes            = auto()
  ScaledStatRes            = auto()
  ScaledSlackRes           = auto()
  PenaltyParameter         = auto()
  MinRayleigh              = auto()
  MaxRayleigh              = auto()
  LastStepOnBoundary       = auto()
  Iteration                = auto()
  LastStepType             = auto()
  ScaledStatResiduals      = auto()
  ScaledFeasResiduals      = auto()
  ScaledConsSlackResiduals = auto()
  ScaledVarSlackResiduals  = auto()


class StepType(Enum):
  NoStep       = csleqp.SLEQP_STEPTYPE_NONE
  Accepted     = csleqp.SLEQP_STEPTYPE_ACCEPTED
  AcceptedFull = csleqp.SLEQP_STEPTYPE_ACCEPTED_FULL
  AcceptedSOC  = csleqp.SLEQP_STEPTYPE_ACCEPTED_SOC
  Rejected     = csleqp.SLEQP_STEPTYPE_REJECTED
