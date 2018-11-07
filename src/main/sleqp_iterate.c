#include "sleqp_iterate.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

SLEQP_RETCODE sleqp_iterate_create(SleqpIterate** star,
                                   SleqpProblem* problem,
                                   SleqpSparseVec* x)
{
  SLEQP_CALL(sleqp_malloc(star));

  assert(sleqp_sparse_vector_valid(x));

  SleqpIterate* iterate = *star;

  int num_variables = problem->num_variables;
  int num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->x,
                                        num_variables,
                                        x->nnz));

  SLEQP_CALL(sleqp_sparse_vector_copy(x, iterate->x));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->func_grad,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->cons_val,
                                        num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_matrix_create(&iterate->cons_jac,
                                        num_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_active_set_create(&iterate->active_set,
                                     problem));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->cons_dual,
                                        num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->vars_dual,
                                        num_variables,
                                        0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_active_set_contains(SleqpIterate* iterate,
                                                SleqpProblem* problem,
                                                SleqpSparseVec* d,
                                                double eps,
                                                double* cache,
                                                bool* contained)
{
  *contained = true;

  SleqpActiveSet* active_set = iterate->active_set;

  const int dim = iterate->x->dim;

  {
    SleqpSparseVec* ub = problem->var_ub;
    SleqpSparseVec* x = iterate->x;

    int k_x = 0, k_d = 0, k_ub = 0;

    while(k_x < x->nnz || k_d < d->nnz || k_ub < ub->nnz)
    {
      bool valid_ub = (k_ub < ub->nnz);
      bool valid_d = (k_d < d->nnz);
      bool valid_x = (k_x < x->nnz);

      int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
      int i_d = valid_d ? d->indices[k_d] : dim + 1;
      int i_x = valid_x ? x->indices[k_x] : dim + 1;

      int i_combined = SLEQP_MIN(i_ub, i_d);
      i_combined = SLEQP_MIN(i_combined, i_x);

      valid_ub = valid_ub && i_combined == i_ub;
      valid_d = valid_d && i_combined == i_d;
      valid_x = valid_x && i_combined == i_x;

      double ubval = valid_ub ? ub->data[k_ub] : 0.;
      double dval = valid_d ? d->data[k_d] : 0.;
      double xval = valid_x ? x->data[k_x] : 0.;

      SLEQP_ACTIVE_STATE var_state = sleqp_active_set_get_variable_state(active_set,
                                                                         i_combined);

      if(var_state == SLEQP_ACTIVE_UPPER)
      {
        if(!sleqp_eq(xval + dval, ubval, eps))
        {
          *contained = false;
          return SLEQP_OKAY;
        }
      }

      if(valid_ub)
      {
        ++k_ub;
      }

      if(valid_d)
      {
        ++k_d;
      }

      if(valid_x)
      {
        ++k_x;
      }

    }
  }

  {
    SleqpSparseVec* lb = problem->var_lb;
    SleqpSparseVec* x = iterate->x;

    int k_x = 0, k_d = 0, k_lb = 0;

    while(k_x < x->nnz || k_d < d->nnz || k_lb < lb->nnz)
    {
      bool valid_lb = (k_lb < lb->nnz);
      bool valid_d = (k_d < d->nnz);
      bool valid_x = (k_x < x->nnz);

      int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
      int i_d = valid_d ? d->indices[k_d] : dim + 1;
      int i_x = valid_x ? x->indices[k_x] : dim + 1;

      int i_combined = SLEQP_MIN(i_lb, i_d);
      i_combined = SLEQP_MIN(i_combined, i_x);

      valid_lb = valid_lb && i_combined == i_lb;
      valid_d = valid_d && i_combined == i_d;
      valid_x = valid_x && i_combined == i_x;

      double lbval = valid_lb ? lb->data[k_lb] : 0.;
      double dval = valid_d ? d->data[k_d] : 0.;
      double xval = valid_x ? x->data[k_x] : 0.;

      SLEQP_ACTIVE_STATE var_state = sleqp_active_set_get_variable_state(active_set,
                                                                         i_combined);

      if(var_state == SLEQP_ACTIVE_LOWER)
      {
        if(!sleqp_eq(xval + dval, lbval, eps))
        {
          *contained = false;
          return SLEQP_OKAY;
        }
      }

      if(valid_lb)
      {
        ++k_lb;
      }

      if(valid_d)
      {
        ++k_d;
      }

      if(valid_x)
      {
        ++k_x;
      }

    }
  }

  {
    SLEQP_CALL(sleqp_sparse_matrix_vector_product(iterate->cons_jac,
                                                  d,
                                                  cache));

    {
      int k_lb = 0, k_ub = 0, k_v = 0;

      SleqpSparseVec* lb = problem->cons_lb;
      SleqpSparseVec* ub = problem->cons_ub;
      SleqpSparseVec* v = iterate->cons_val;

      while(k_lb < lb->nnz ||
            k_ub < ub->nnz ||
            k_v < v->nnz)
      {
        bool valid_lb = (k_lb < lb->nnz);
        bool valid_ub = (k_ub < ub->nnz);
        bool valid_v = (k_v < v->nnz);

        int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        int i_v = valid_v ? v->indices[k_v] : dim + 1;

        int i_combined = SLEQP_MIN(i_lb, i_ub);
        i_combined = SLEQP_MIN(i_combined, i_v);

        valid_lb = valid_lb && i_combined == i_lb;
        valid_ub = valid_ub && i_combined == i_ub;
        valid_v = valid_v && i_combined == i_v;

        double lbval = valid_lb ? lb->data[k_lb] : 0.;
        double ubval = valid_ub ? ub->data[k_ub] : 0.;
        double vval = valid_v ? v->data[k_v] : 0.;

        SLEQP_ACTIVE_STATE cons_state = sleqp_active_set_get_constraint_state(active_set,
                                                                              i_combined);

        if(cons_state == SLEQP_ACTIVE_UPPER)
        {
          if(!sleqp_eq(vval + cache[i_combined], ubval, eps))
          {
            *contained = false;
            return SLEQP_OKAY;
          }
        }

        if(cons_state == SLEQP_ACTIVE_LOWER)
        {
          if(!sleqp_eq(vval + cache[i_combined], lbval, eps))
          {
            *contained = false;
            return SLEQP_OKAY;
          }
        }

        if(valid_lb)
        {
          ++k_lb;
        }

        if(valid_ub)
        {
          ++k_ub;
        }

        if(valid_v)
        {
          ++k_v;
        }
      }

    }


  }

  return SLEQP_OKAY;
}

double sleqp_iterate_slackness_residuum(SleqpIterate* iterate,
                                        SleqpProblem* problem)
{
  double residuum = 0.;

  {
    SleqpSparseVec* x = iterate->x;

    SleqpSparseVec* lb = problem->var_lb;
    SleqpSparseVec* ub = problem->var_ub;

    SleqpSparseVec* d = iterate->vars_dual;

    {
      int k_x = 0, k_ub = 0, k_d = 0;

      int dim = x->dim;

      while((k_x < x->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_x = k_x < x->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_d = k_d < d->nnz;

        int i_x = valid_x ? x->indices[k_x] : dim + 1;
        int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_x, i_ub);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_x = valid_x && (i_x == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double x_val = valid_x ? x->data[k_x] : 0.;
        double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MAX(ub_val - x_val, 0.) * SLEQP_MAX(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_x)
        {
          ++k_x;
        }

        if(valid_ub)
        {
          ++k_ub;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }

    {
      int k_x = 0, k_lb = 0, k_d = 0;

      int dim = x->dim;

      while((k_x < x->nnz || k_lb < lb->nnz) && k_d < d->nnz)
      {
        bool valid_x = k_x < x->nnz;
        bool valid_lb = k_lb < lb->nnz;
        bool valid_d = k_d < d->nnz;

        int i_x = valid_x ? x->indices[k_x] : dim + 1;
        int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_x, i_lb);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_x = valid_x && (i_x == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double x_val = valid_x ? x->data[k_x] : 0.;
        double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MIN(lb_val - x_val, 0.) * SLEQP_MIN(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_x)
        {
          ++k_x;
        }

        if(valid_lb)
        {
          ++k_lb;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  {
    SleqpSparseVec* v = iterate->cons_val;

    SleqpSparseVec* lb = problem->cons_lb;
    SleqpSparseVec* ub = problem->cons_ub;

    SleqpSparseVec* d = iterate->cons_dual;

    {
      int k_v = 0, k_ub = 0, k_d = 0;

      int dim = v->dim;

      while((k_v < v->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_v = k_v < v->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_d = k_d < d->nnz;

        int i_v = valid_v ? v->indices[k_v] : dim + 1;
        int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_v, i_ub);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_v = valid_v && (i_v == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double v_val = valid_v ? v->data[k_v] : 0.;
        double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MAX(ub_val - v_val, 0.) * SLEQP_MAX(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_v)
        {
          ++k_v;
        }

        if(valid_ub)
        {
          ++k_ub;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }


    {
      int k_v = 0, k_lb = 0, k_d = 0;

      int dim = v->dim;

      while((k_v < v->nnz || k_lb < lb->nnz) && k_d < d->nnz)
      {
        bool valid_v = k_v < v->nnz;
        bool valid_lb = k_lb < lb->nnz;
        bool valid_d = k_d < d->nnz;

        int i_v = valid_v ? v->indices[k_v] : dim + 1;
        int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_v, i_lb);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_v = valid_v && (i_v == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double v_val = valid_v ? v->data[k_v] : 0.;
        double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MIN(lb_val - v_val, 0.) * SLEQP_MIN(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_v)
        {
          ++k_v;
        }

        if(valid_lb)
        {
          ++k_lb;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  return residuum;
}


SLEQP_RETCODE sleqp_iterate_free(SleqpIterate** star)
{
  SleqpIterate* iterate = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->vars_dual));
  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->cons_dual));

  SLEQP_CALL(sleqp_active_set_free(&iterate->active_set));

  SLEQP_CALL(sleqp_sparse_matrix_free(&iterate->cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->func_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->x));

  sleqp_free(star);

  return SLEQP_OKAY;
}
