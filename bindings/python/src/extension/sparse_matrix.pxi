#cython: language_level=3

cimport numpy as np

from scipy.sparse import csc_matrix

cdef class MatrixRef:
  cdef csleqp.SleqpMat* matrix

  cdef _release(self):
    if self.matrix:
      csleqp_call(csleqp.sleqp_mat_release(&self.matrix))

  cdef _set_matrix(self, csleqp.SleqpMat* matrix):
    self._release()

    csleqp_call(csleqp.sleqp_mat_capture(matrix))

    self.matrix = matrix

  def __dealloc__(self):
    self._release()

cdef object sleqp_sparse_matrix_to_scipy(csleqp.SleqpMat* _matrix):
  cdef MatrixRef matrix_ref = MatrixRef()

  matrix_ref._set_matrix(_matrix)

  num_cols = csleqp.sleqp_mat_num_cols(_matrix)

  num_rows = csleqp.sleqp_mat_num_rows(_matrix)

  nnz = csleqp.sleqp_mat_nnz(_matrix)

  data = np.asarray(<double[:nnz]> csleqp.sleqp_mat_data(_matrix))

  indices = np.asarray(<int[:nnz]> csleqp.sleqp_mat_rows(_matrix))

  indptr = np.asarray(<int[:num_cols + 1]> csleqp.sleqp_mat_cols(_matrix))

  matrix = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))

  matrix._ref = matrix_ref

  return matrix
