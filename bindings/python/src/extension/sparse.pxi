#cython: language_level=3


cdef object sleqp_sparse_vec_to_array(const csleqp.SleqpVec* vec):
  assert vec
  values = np.zeros((vec.dim,), dtype=np.float64)

  for k in range(vec.nnz):
    values[vec.indices[k]] = vec.data[k]

  return values


@cython.boundscheck(False)
@cython.wraparound(False)
cdef csleqp.SLEQP_RETCODE array_set_values(np.ndarray array,
                                           csleqp.SleqpVec* vec) except csleqp.SLEQP_ERROR:
  cdef Py_ssize_t i
  cdef Py_ssize_t n = array.shape[0]
  cdef double v

  for i in range(n):
    v = array[i]
    if v != 0.:
      csleqp_call(csleqp.sleqp_vec_push(vec, i, v))


cdef csleqp.SLEQP_RETCODE array_to_sleqp_sparse_vec(np.ndarray array,
                                                    csleqp.SleqpVec* vec) except csleqp.SLEQP_ERROR:
  assert vec != NULL

  csleqp_call(csleqp.sleqp_vec_clear(vec))

  if array is None:
    return csleqp.SLEQP_OKAY

  assert array.ndim == 1

  cdef int dim = array.shape[0]

  csleqp_call(csleqp.sleqp_vec_reserve(vec, dim))

  csleqp_call(array_set_values(array, vec))

  return csleqp.SLEQP_OKAY

class DenseIterator:
  def __init__(self, matrix):
    self.matrix = matrix
    self.iter = np.nditer(self.matrix,
                          flags=['multi_index'],
                          order='F')

  def __iter__(self):
    while not self.iter.finished:
      data = self.iter[0]
      if not data:
        self.iter.iternext()
        continue

      yield self.iter.multi_index, data

      self.iter.iternext()

  def length_bound(self):
    return np.prod(self.matrix.shape)


class CSCIterator:
  def __init__(self, matrix):
    assert scipy.sparse.isspmatrix_csc(matrix)
    self.matrix = matrix
    self._length_bound = -1

  def __iter__(self):
    for col in range(self.matrix.shape[1]):
      for ind in range(self.matrix.indptr[col], self.matrix.indptr[col+1]):
        data = self.matrix.data[ind]
        if data:
          yield (self.matrix.indices[ind], col), self.matrix.data[ind]

  def length_bound(self):
    if self._length_bound != -1:
      return self._length_bound

    self._length_bound = self.matrix.count_nonzero()

    return self._length_bound

def matrix_iterator(matrix):
  assert not (matrix is None)
  assert matrix.ndim == 2

  if not scipy.sparse.issparse(matrix):
    return DenseIterator(matrix)

  assert scipy.sparse.issparse(matrix)

  if not scipy.sparse.isspmatrix_csc(matrix):
    matrix = scipy.sparse.csc_matrix(matrix)

  assert scipy.sparse.isspmatrix_csc(matrix)

  return CSCIterator(matrix)

cdef csleqp.SLEQP_RETCODE matrix_to_sleqp_sparse_matrix(object mat,
                                                        csleqp.SleqpMat* matrix) \
                                                        except csleqp.SLEQP_ERROR:
  assert matrix

  if mat is None:
    return csleqp.SLEQP_OKAY

  assert mat.ndim == 2

  num_rows = csleqp.sleqp_mat_num_rows(matrix)
  num_cols = csleqp.sleqp_mat_num_cols(matrix)

  assert mat.shape == (num_rows, num_cols)

  if num_rows == 0 or num_cols == 0:
    return csleqp.SLEQP_OKAY

  cdef int last_col = -1

  matrix_iter = matrix_iterator(mat)

  for (row, col), data in matrix_iter:
    while last_col != col:
      last_col += 1
      csleqp_call(csleqp.sleqp_mat_push_col(matrix, last_col))

    nnz = csleqp.sleqp_mat_nnz(matrix)
    nnz_max = csleqp.sleqp_mat_nnz_max(matrix)

    if nnz == nnz_max:
      csleqp_call(csleqp.sleqp_mat_reserve(matrix,
                                           matrix_iter.length_bound()))

    csleqp_call(csleqp.sleqp_mat_push(matrix,
                                      row,
                                      col,
                                      data))

  while last_col < csleqp.sleqp_mat_num_cols(matrix) - 1:
    last_col += 1
    csleqp_call(csleqp.sleqp_mat_push_col(matrix, last_col))

  return csleqp.SLEQP_OKAY
