/**
 * @file bqpd.h
 * @brief Glue for interfacing the Fortran BQPD solver
 * @author Christian Kirches
 * @date 2017-02-05
 */

#pragma once
#ifndef _BQPD_H_INCLUDED_
#define _BQPD_H_INCLUDED_

#ifdef __cplusplus
extern "C" {
#endif

#include "util_types.h"     // fint_t, freal_t

extern void
bqpd_ (
    fint_t *n,
    fint_t *m,
    fint_t *k,
    fint_t *kmax,
    freal_t *a,
    fint_t *la,
    freal_t *x,
    freal_t *bl,
    freal_t *bu,
    freal_t *f,
    freal_t *fmin,
    freal_t *g,
    freal_t *r,
    freal_t *w,
    freal_t *e,
    fint_t *ls,
    fint_t *alp,
    fint_t *lp,
    fint_t *mlp,
    fint_t *peq,
    freal_t *ws,
    fint_t *lws,
    fint_t *mode,
    fint_t *ifail,
    fint_t *info,
    fint_t *iprint,
    fint_t *nout
);

#ifdef __cplusplus
}
#endif

#endif // _BQPD_H_INCLUDED_
