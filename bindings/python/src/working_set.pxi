#cython: language_level=3

cdef class WorkingSet:
  cdef csleqp.SleqpWorkingSet* working_set

  cdef dict __dict__

  cdef _release(self):
    if self.working_set:
      csleqp_call(csleqp.sleqp_working_set_release(&self.working_set))

  cdef _set_working_set(self, csleqp.SleqpWorkingSet* working_set):
    self._release()

    csleqp_call(csleqp.sleqp_working_set_capture(working_set))

    self.working_set = working_set

  @property
  def num_variables(self) -> int:
    cdef csleqp.SleqpProblem* problem

    problem = csleqp.sleqp_working_set_get_problem(self.working_set)

    return problem.num_variables

  @property
  def num_constraints(self) -> int:
    cdef csleqp.SleqpProblem* problem

    problem = csleqp.sleqp_working_set_get_problem(self.working_set)

    return problem.num_constraints

  @property
  def num_active_variables(self) -> int:
    assert self.working_set

    return csleqp.sleqp_working_set_num_active_vars(self.working_set)

  @property
  def num_active_constraints(self) -> int:
    assert self.working_set

    return csleqp.sleqp_working_set_num_active_cons(self.working_set)

  def variable_state(self, int i) -> ActiveState:
    assert self.working_set

    state = csleqp.sleqp_working_set_get_variable_state(self.working_set,
                                                        i)

    return ActiveState(state)

  def constraint_state(self, int i) -> ActiveState:
    assert self.working_set

    state = csleqp.sleqp_working_set_get_constraint_state(self.working_set,
                                                          i)

    return ActiveState(state)

  def __str__(self):
    str_val = "Active set, variables: {0}, constraints: {1}\n".format(self.num_variables,
                                                                      self.num_constraints)

    str_val += "\n".join(("State of variable {0}: {1}".format(j, str(self.variable_state(j)))
                          for j in range(self.num_variables)))

    str_val += "\n"

    str_val += "\n".join(("State of constraint {0}: {1}".format(i, str(self.constraint_state(i)))
                          for i in range(self.num_constraints)))

    return str_val

  @property
  def size(self) -> int:
    assert self.working_set

    return csleqp.sleqp_working_set_size(self.working_set)

  def __dealloc__(self):
    self._release()
