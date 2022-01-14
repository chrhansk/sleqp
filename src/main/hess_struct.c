#include "hess_struct.h"

#include "error.h"
#include "log.h"
#include "mem.h"

struct SleqpHessStruct
{
  int refcount;
  int* block_ends;
  int num_blocks;
  int dimension;
};

SLEQP_RETCODE
sleqp_hess_struct_create(SleqpHessStruct** star, int dimension, bool empty)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpHessStruct* hessian_struct = *star;

  *hessian_struct = (SleqpHessStruct){0};

  hessian_struct->refcount = 1;

  SLEQP_CALL(sleqp_alloc_array(&hessian_struct->block_ends, dimension));

  hessian_struct->dimension = dimension;

  if (!empty && dimension > 0)
  {
    SLEQP_CALL(sleqp_hess_struct_push_block(hessian_struct, dimension));
  }

  return SLEQP_OKAY;
}

int
sleqp_hess_struct_num_blocks(const SleqpHessStruct* hessian_struct)
{
  return hessian_struct->num_blocks;
}

SLEQP_RETCODE
sleqp_hess_struct_push_block(SleqpHessStruct* hessian_struct, int end)
{
  if (end > hessian_struct->dimension)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid block index (%d)", end);
  }

  const int num_blocks = sleqp_hess_struct_num_blocks(hessian_struct);

  if (num_blocks >= hessian_struct->dimension)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid block index (%d)", end);
  }

  if (num_blocks > 0 && end <= hessian_struct->block_ends[num_blocks - 1])
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid block index (%d)", end);
  }

  hessian_struct->block_ends[num_blocks] = end;
  ++(hessian_struct->num_blocks);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_hess_struct_clear(SleqpHessStruct* hessian_struct)
{
  hessian_struct->num_blocks = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_hess_struct_block_range(const SleqpHessStruct* hessian_struct,
                              int block,
                              int* begin,
                              int* end)
{
  if (block > hessian_struct->num_blocks)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid block index (%d)", block);
  }

  *begin = (block > 0) ? hessian_struct->block_ends[block - 1] : 0;

  *end = hessian_struct->block_ends[block];

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_hess_struct_lin_range(const SleqpHessStruct* hessian_struct,
                            int* begin,
                            int* end)
{
  const int num_blocks = sleqp_hess_struct_num_blocks(hessian_struct);

  *begin = (num_blocks == 0) ? 0 : hessian_struct->block_ends[num_blocks - 1];

  *end = hessian_struct->dimension;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_hess_struct_fprintf(SleqpHessStruct* hessian_struct, FILE* output)
{
  fprintf(output,
          "Hessian structure, dimension: %d, nonlinear blocks: %d\n",
          hessian_struct->dimension,
          hessian_struct->num_blocks);

  int begin, end;

  const int num_blocks = sleqp_hess_struct_num_blocks(hessian_struct);

  for (int block = 0; block < num_blocks; ++block)
  {

    SLEQP_CALL(
      sleqp_hess_struct_block_range(hessian_struct, block, &begin, &end));

    fprintf(output, "Block %d: [%d, %d)\n", block, begin, end);
  }

  SLEQP_CALL(sleqp_hess_struct_lin_range(hessian_struct, &begin, &end));

  if (begin < end)
  {
    fprintf(output, "Linear range: [%d, %d)\n", begin, end);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
hessian_struct_free(SleqpHessStruct** star)
{
  SleqpHessStruct* hessian_struct = *star;

  if (!hessian_struct)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&hessian_struct->block_ends);

  sleqp_free(&hessian_struct);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_hess_struct_copy(const SleqpHessStruct* source, SleqpHessStruct* target)
{
  target->dimension = source->dimension;

  SLEQP_CALL(sleqp_hess_struct_clear(target));

  const int num_blocks = sleqp_hess_struct_num_blocks(source);

  for (int block = 0; block < num_blocks; ++block)
  {
    int begin;
    int end;

    SLEQP_CALL(sleqp_hess_struct_block_range(source, block, &begin, &end));

    SLEQP_CALL(sleqp_hess_struct_push_block(target, end));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_hess_struct_capture(SleqpHessStruct* hessian_struct)
{
  ++hessian_struct->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_hess_struct_release(SleqpHessStruct** star)
{
  SleqpHessStruct* hessian_struct = *star;

  if (!hessian_struct)
  {
    return SLEQP_OKAY;
  }

  if (--hessian_struct->refcount == 0)
  {
    SLEQP_CALL(hessian_struct_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
