#cython: language_level=3

cdef object sleqp_sparse_vec_to_array(csleqp.SleqpSparseVec* vec):
  assert vec
  values = np.zeros([vec.dim], dtype=np.float64)

  for k in range(vec.nnz):
    data = vec.data[k]
    if data:
      values[vec.indices[k]] = vec.data[k]

  return values

cdef csleqp.SLEQP_RETCODE array_to_sleqp_sparse_vec(np.ndarray array,
                                                    csleqp.SleqpSparseVec* vec) except csleqp.SLEQP_INTERNAL_ERROR:
  assert vec != NULL

  if array is None:
    return csleqp.SLEQP_OKAY

  cdef int dim = array.shape[0]

  csleqp_call(csleqp.sleqp_sparse_vector_reserve(vec, dim))
  csleqp_call(csleqp.sleqp_sparse_vector_clear(vec))

  for i in range(dim):
    csleqp_call(csleqp.sleqp_sparse_vector_push(vec, i, array[i]))

  return csleqp.SLEQP_OKAY

def iter_matrix_entries(matrix):

  if matrix is None:
    return

  m = matrix

  assert matrix.ndim == 2

  if not scipy.sparse.issparse(matrix):
    it = np.nditer(matrix, flags=['multi_index'], order='F')
    while not it.finished:
      if it[0] != 0:
        yield *it.multi_index, it[0]
      it.iternext()

    return

  if not scipy.sparse.isspmatrix_csc(m):
    m = scipy.sparse.csc_matrix(matrix)

  assert scipy.sparse.isspmatrix_csc(m)

  for col in range(m.shape[1]):
    for ind in range(m.indptr[col], m.indptr[col+1]):
      data = m.data[ind]
      if data:
        yield m.indices[ind], col, m.data[ind]

cdef csleqp.SLEQP_RETCODE matrix_to_sleqp_sparse_matrix(object mat,
                                                        csleqp.SleqpSparseMatrix* matrix) except csleqp.SLEQP_INTERNAL_ERROR:

  assert matrix

  cdef int last_col = -1

  for (row, col, data) in iter_matrix_entries(mat):
    while last_col != col:
      last_col += 1
      csleqp_call(csleqp.sleqp_sparse_matrix_push_column(matrix, last_col))

    csleqp_call(csleqp.sleqp_sparse_matrix_push(matrix,
                                                row,
                                                col,
                                                data))

  while last_col < matrix.num_cols - 1:
    last_col += 1
    csleqp_call(csleqp.sleqp_sparse_matrix_push_column(matrix, last_col))

  return csleqp.SLEQP_OKAY
