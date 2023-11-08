cdef class HessianStruct:
    cdef csleqp.SleqpHessStruct * hess_struct
    cdef dict __dict__
    cdef _Func _func

    def __cinit__(self, _Func func):
        self._func = func
        self.hess_struct = csleqp.sleqp_func_hess_struct(func.cfunc)
        csleqp_call(csleqp.sleqp_hess_struct_capture(self.hess_struct))

    def __dealloc__(self):
        csleqp_call(csleqp.sleqp_hess_struct_release( & self.hess_struct))

    def clear(self) -> None:
        csleqp_call(csleqp.sleqp_hess_struct_clear(self.hess_struct))

    def push(self, int end) -> None:
        csleqp_call(csleqp.sleqp_hess_struct_push_block(self.hess_struct, end))

    @property
    def num_blocks(self) -> int:
        return csleqp.sleqp_hess_struct_num_blocks(self.hess_struct)

    def block_range(self, int block) -> collections.abc.Iterable[int]:
        cdef int begin = 0
        cdef int end = 0

        csleqp_call(csleqp.sleqp_hess_struct_block_range(self.hess_struct,
                                                         block,
                                                         & begin,
                                                         & end))

        return range(begin, end)

    @property
    def block_ranges(self) -> typing.Iterable[collections.abc.Iterable[int]]:

        num_blocks = self.num_blocks

        for block in range(num_blocks):
            yield self.block_range(block)

    @property
    def linear_range(self) -> collections.abc.Iterable[int]:
        cdef int begin = 0
        cdef int end = 0

        csleqp_call(csleqp.sleqp_hess_struct_lin_range(self.hess_struct,
                                                       & begin,
                                                       & end))

        return range(begin, end)

    def __str__(self):
        blocks = ', '.join(('[{0}, {1})'.format(block_range.start,
                                                block_range.stop)
                            for block_range in self.block_ranges))

        return 'HessianStruct: [{0}]'.format(blocks)
