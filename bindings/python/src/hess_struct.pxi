cdef class HessianStruct:
  cdef csleqp.SleqpHessianStruct* hess_struct
  cdef dict __dict__

  def __cinit__(self, Func func):
    self._func = func
    self.hess_struct = csleqp.sleqp_func_get_hess_struct(func.func)
    csleqp_call(csleqp.sleqp_hessian_struct_capture(self.hess_struct))

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_hessian_struct_release(&self.hess_struct))

  def clear(self) -> None:
    csleqp_call(csleqp.sleqp_hessian_struct_clear(self.hess_struct))

  def push(self, int end) -> None:
    csleqp_call(csleqp.sleqp_hessian_struct_push_block(self.hess_struct, end))

  @property
  def num_blocks(self) -> int:
    return csleqp.sleqp_hessian_struct_get_num_blocks(self.hess_struct)

  def block_range(self, int block) -> typing.Tuple[int, int]:
    cdef int begin = 0
    cdef int end = 0

    csleqp_call(csleqp.sleqp_hessian_struct_get_block_range(self.hess_struct,
                                                            block,
                                                            &begin,
                                                            &end))

    return (begin, end)

  @property
  def block_ranges(self) -> typing.Iterable[typing.Tuple[int, int]]:

    num_blocks = self.num_blocks

    for block in range(num_blocks):
        yield self.block_range(block)

  @property
  def linear_range(self) -> typing.Tuple[int, int]:
    cdef int begin = 0
    cdef int end = 0

    csleqp_call(csleqp.sleqp_hessian_struct_get_linear_range(self.hess_struct,
                                                             &begin,
                                                             &end))

    return (begin, end)

  def __str__(self):
    blocks = ', '.join(('[{0}, {1})'.format(*block_range)
                        for block_range in self.block_ranges))

    return 'HessianStruct: [{0}]'.format(blocks)
