#include "hess_struct.h"

#include "log.h"
#include "mem.h"

struct SleqpHessianStruct
{
  int refcount;
  int* block_ends;
  int num_blocks;
  int dimension;
};

SLEQP_RETCODE sleqp_hessian_struct_create(SleqpHessianStruct** star,
                                          int dimension,
                                          bool empty)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpHessianStruct* hessian_struct = *star;

  *hessian_struct = (SleqpHessianStruct) {0};

  hessian_struct->refcount = 1;

  SLEQP_CALL(sleqp_alloc_array(&hessian_struct->block_ends, dimension));

  hessian_struct->dimension = dimension;

  if(!empty && dimension > 0)
  {
    SLEQP_CALL(sleqp_hessian_struct_push_block(hessian_struct,
                                               dimension));
  }

  return SLEQP_OKAY;
}

int sleqp_hessian_struct_get_num_blocks(const SleqpHessianStruct* hessian_struct)
{
  return hessian_struct->num_blocks;
}

SLEQP_RETCODE sleqp_hessian_struct_push_block(SleqpHessianStruct* hessian_struct,
                                              int end)
{
  if(end > hessian_struct->dimension)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  const int num_blocks = sleqp_hessian_struct_get_num_blocks(hessian_struct);

  if(num_blocks >= hessian_struct->dimension)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  if(num_blocks > 0 && end <= hessian_struct->block_ends[num_blocks - 1])
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  hessian_struct->block_ends[num_blocks] = end;
  ++(hessian_struct->num_blocks);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hessian_struct_clear(SleqpHessianStruct* hessian_struct)
{
  hessian_struct->num_blocks = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hessian_struct_get_block_range(const SleqpHessianStruct* hessian_struct,
                                                   int block,
                                                   int* begin,
                                                   int* end)
{
  if(block > hessian_struct->num_blocks)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  *begin = (block > 0) ? hessian_struct->block_ends[block - 1] : 0;

  *end = hessian_struct->block_ends[block];

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hessian_struct_get_linear_range(const SleqpHessianStruct* hessian_struct,
                                                    int* begin,
                                                    int* end)
{
  const int num_blocks = sleqp_hessian_struct_get_num_blocks(hessian_struct);

  *begin = (num_blocks == 0) ? 0 : hessian_struct->block_ends[num_blocks - 1];

  *end = hessian_struct->dimension;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hessian_struct_fprintf(SleqpHessianStruct* hessian_struct,
                                           FILE* output)
{
  fprintf(output,
          "Hessian structure, dimension: %d, nonlinear blocks: %d\n",
          hessian_struct->dimension,
          hessian_struct->num_blocks);

  int begin, end;

  const int num_blocks = sleqp_hessian_struct_get_num_blocks(hessian_struct);

  for(int block = 0; block < num_blocks; ++block)
  {

    SLEQP_CALL(sleqp_hessian_struct_get_block_range(hessian_struct,
                                                    block,
                                                    &begin,
                                                    &end));

    fprintf(output,
            "Block %d: [%d, %d)\n",
            block,
            begin,
            end);
  }

  SLEQP_CALL(sleqp_hessian_struct_get_linear_range(hessian_struct,
                                                   &begin,
                                                   &end));

  if(begin < end)
  {
    fprintf(output,
            "Linear range: [%d, %d)\n",
            begin,
            end);
  }


  return SLEQP_OKAY;
}

static SLEQP_RETCODE hessian_struct_free(SleqpHessianStruct** star)
{
  SleqpHessianStruct* hessian_struct = *star;

  if(!hessian_struct)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&hessian_struct->block_ends);

  sleqp_free(&hessian_struct);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hessian_struct_copy(const SleqpHessianStruct* source,
                                        SleqpHessianStruct* target)
{
  target->dimension = source->dimension;

  SLEQP_CALL(sleqp_hessian_struct_clear(target));

  const int num_blocks = sleqp_hessian_struct_get_num_blocks(source);

  for(int block = 0; block < num_blocks; ++block)
  {
    int begin;
    int end;

    SLEQP_CALL(sleqp_hessian_struct_get_block_range(source,
                                                    block,
                                                    &begin,
                                                    &end));

    SLEQP_CALL(sleqp_hessian_struct_push_block(target, end));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hessian_struct_capture(SleqpHessianStruct* hessian_struct)
{
  ++hessian_struct->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hessian_struct_release(SleqpHessianStruct** star)
{
  SleqpHessianStruct* hessian_struct = *star;

  if(!hessian_struct)
  {
    return SLEQP_OKAY;
  }

  if(--hessian_struct->refcount == 0)
  {
    SLEQP_CALL(hessian_struct_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
