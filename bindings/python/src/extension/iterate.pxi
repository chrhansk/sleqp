# cython: language_level=3

cdef class Iterate:
    cdef csleqp.SleqpIterate * iterate

    cdef dict __dict__

    def __cinit__(self, **properties):
        if not "_create" in properties:
            raise Exception("Manually called constructor")

    cdef _release(self):
        if self.iterate:
            csleqp_call(csleqp.sleqp_iterate_release( & self.iterate))

    cdef _set_iterate(self, csleqp.SleqpIterate * iterate):
        self._release()

        csleqp_call(csleqp.sleqp_iterate_capture(iterate))

        self.iterate = iterate

    def __dealloc__(self):
        self._release()

    @property
    def primal(self) -> np.array:
        assert self.iterate

        return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_primal(self.iterate))

    @property
    def vars_dual(self) -> np.array:
        assert self.iterate

        return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_vars_dual(self.iterate))

    @property
    def cons_dual(self) -> np.array:
        assert self.iterate

        return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_cons_dual(self.iterate))

    @property
    def working_set(self) -> WorkingSet:
        cdef csleqp.SleqpWorkingSet * _working_set
        cdef WorkingSet working_set = WorkingSet()

        assert self.iterate

        _working_set = csleqp.sleqp_iterate_working_set(self.iterate)

        working_set._set_working_set(_working_set)

        return working_set

    @property
    def cons_jac(self) -> scipy.sparse.csc_matrix:
        assert self.iterate

        return sleqp_sparse_matrix_to_scipy(csleqp.sleqp_iterate_cons_jac(self.iterate))

    @property
    def obj_val(self) -> float:
        assert self.iterate

        return csleqp.sleqp_iterate_obj_val(self.iterate)

    @property
    def obj_grad(self) -> np.array:
        assert self.iterate

        return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_obj_grad(self.iterate))

    @property
    def cons_val(self) -> np.array:
        assert self.iterate

        return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_cons_val(self.iterate))
