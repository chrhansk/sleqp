#!/usr/bin/python
#cython: language_level=3

cdef object sleqp_sparse_vec_to_array(csleqp.SleqpSparseVec* vec):
  assert vec
  values = np.zeros([vec.dim], dtype=np.float64)

  for k in range(vec.nnz):
    values[vec.indices[k]] = vec.data[k]

  return values

cdef void array_to_sleqp_sparse_vec(np.ndarray array, csleqp.SleqpSparseVec* vec):
  assert vec

  if array is None:
    return

  cdef int dim = array.shape[0]

  csleqp_call(csleqp.sleqp_sparse_vector_reserve(vec, dim))

  for i in range(dim):
    csleqp_call(csleqp.sleqp_sparse_vector_push(vec, i, array[i]))

# cdef void sparse_to_sleqp_sparse_matrix(np.ndarray mat, csleqp.SleqpSparseMatrix* matrix):
#   assert matrix
#   mcsc = mat.tocsc()

#   if not mcsc.has_sorted_indices():
#     mcsc = mcsc.sort_indices()

#   assert mcsc.shape == (matrix.num_cols, matrix.num_rows)

#   csleqp.sleqp_sparse_matrix_reserve(matrix, mcsc.nnz)

#   for c in range(matrix.num_cols):
#     matrix.cols[c] = mcsc.indptr[c]

#   for i in range(mcsc.nnz):
#     matrix.rows[i] = mcsc.indices[i]
#     matrix.data[i] = mcsc.data[i]
