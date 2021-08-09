#ifndef SLEQP_SPARSE_VEC_H
#define SLEQP_SPARSE_VEC_H

#include "pub_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_sparse_vector_remove_entries(const SleqpSparseVec* source,
                                                   SleqpSparseVec* target,
                                                   const int* entry_indices,
                                                   int num_entries);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_VEC_H */
