#!/usr/bin/python
#cython: language_level=3

cdef object sleqp_sparse_vec_to_array(csleqp.SleqpSparseVec* vec):
  assert vec
  values = np.zeros([vec.dim], dtype=np.float64)

  for k in range(vec.nnz):
    data = vec.data[k]
    if data:
      values[vec.indices[k]] = vec.data[k]

  return values

cdef void array_to_sleqp_sparse_vec(np.ndarray array, csleqp.SleqpSparseVec* vec):
  assert vec

  if array is None:
    return

  cdef int dim = array.shape[0]

  csleqp_call(csleqp.sleqp_sparse_vector_reserve(vec, dim))
  csleqp_call(csleqp.sleqp_sparse_vector_clear(vec))

  for i in range(dim):
    csleqp_call(csleqp.sleqp_sparse_vector_push(vec, i, array[i]))


def iter_matrix_entries(matrix):

  m = matrix

  if not scipy.sparse.issparse(matrix):
    m = scipy.sparse.csc_matrix(matrix)

  assert scipy.sparse.isspmatrix_csc(m)

  #if not scipy.sparse.isspmatrix_csc(m):
  #  m = m.tocss()

  for col in range(m.shape[1]):
    for ind in range(m.indptr[col], m.indptr[col+1]):
      data = m.data[ind]
      if data:
        yield m.indices[ind], col, m.data[ind]

cdef void matrix_to_sleqp_sparse_matrix(object mat, csleqp.SleqpSparseMatrix* matrix):

  assert matrix

  cold = -1

  for (row, col, data) in iter_matrix_entries(mat):
    if cold != col:
      cold = col
      csleqp_call(csleqp.sleqp_sparse_matrix_add_column(matrix, col))

    csleqp_call(csleqp.sleqp_sparse_matrix_push(matrix,
                                                row,
                                                col,
                                                data))
