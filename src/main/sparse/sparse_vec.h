#ifndef SLEQP_SPARSE_VEC_H
#define SLEQP_SPARSE_VEC_H

#include "pub_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Dumps this vector to the given file. The vector
   * is printed out as a set of lines each consisting
   * one entry
   *
   * @param[in]  vec     A pointer to the vector
   * @param[in]  output  A pointer to an output `FILE*`
   **/
  SLEQP_RETCODE sleqp_sparse_vector_dump(const SleqpSparseVec* vec,
                                         FILE* output);

  SLEQP_RETCODE sleqp_sparse_vector_dump_to_file(const SleqpSparseVec* vec,
                                                 const char* name);

  SLEQP_RETCODE sleqp_sparse_vector_remove_entries(const SleqpSparseVec* source,
                                                   SleqpSparseVec* target,
                                                   const int* entry_indices,
                                                   int num_entries);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_VEC_H */
