#ifndef SLEQP_VEC_H
#define SLEQP_VEC_H

#include "pub_vec.h"

/**
 * Dumps this vector to the given file. The vector
 * is printed out as a set of lines each consisting
 * one entry
 *
 * @param[in]  vec     A pointer to the vector
 * @param[in]  output  A pointer to an output `FILE*`
 **/
SLEQP_RETCODE
sleqp_vec_dump(const SleqpVec* vec, FILE* output);

SLEQP_RETCODE
sleqp_vec_dump_to_file(const SleqpVec* vec, const char* name);

SLEQP_RETCODE
sleqp_vec_remove_entries(const SleqpVec* source,
                         SleqpVec* target,
                         const int* entry_indices,
                         int num_entries);

#endif /* SLEQP_VEC_H */
