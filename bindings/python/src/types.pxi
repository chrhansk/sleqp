class Status(enum.Enum):
  Optimal = csleqp.SLEQP_OPTIMAL
  Feasible = csleqp.SLEQP_FEASIBLE
  Infeasible = csleqp.SLEQP_INFEASIBLE
  Invalid = csleqp.SLEQP_INVALID
