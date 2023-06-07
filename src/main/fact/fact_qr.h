#ifndef SLEQP_FACT_QR_H
#define SLEQP_FACT_QR_H

#include "fact_qr_types.h"
#include "pub_settings.h"
#include "settings.h"

/**
 * Create data structure modeling a (full) QR
 * factorization of a matrix \f$ A \in \mathbb{R}^{m \times n} \f$
 * with \f$ m \geq n\f$,
 * i.e., \f$ A = QR \f$,
 * where \f$ Q \in \mathbb{R}^{m \times m} \f$ and
 * \f$ R \in \mathbb{R}^{m \times n} \f$,
 * where \f$ R = \begin{pmatrix} \overline{R} \\ 0 \end{pmatrix}\f$
 * with \f$ \overline{R} \in \mathbb{R}^{n \times n}\f$
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_qr_create(SleqpFactQR** star,
                const char* name,
                const char* version,
                SleqpSettings* settings,
                SleqpQRCallbacks* callbacks,
                void* fact_data);

const char*
sleqp_qr_name(SleqpFactQR* fact);

const char*
sleqp_qr_version(SleqpFactQR* fact);

/**
 * Sets the matrix \f$ A \f$
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_qr_set_matrix(SleqpFactQR* fact, SleqpMat* matrix);

/**
 * Solve the system \f$ \overline{R} x = y \f$
 *
 * @param[in]  rhs The right-hand side \f$ y \in \mathbb{R}^{m} \f$
 * @param[out] sol The solution \f$ x \in \mathbb{R}^{m} \f$
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_qr_solve_tri(SleqpFactQR* fact, const SleqpVec* rhs, SleqpVec* sol);

/**
 * Solve the system \f$ \overline{R} x = y \f$
 *
 * @param[in]  rhs The right-hand side \f$ y \in \mathbb{R}^{m} \f$
 * @param[out] sol The solution \f$ x \in \mathbb{R}^{m} \f$
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_qr_solve_tri_trans(SleqpFactQR* fact, const SleqpVec* rhs, SleqpVec* sol);

/**
 * Compute \f$y = Q x \f$
 *
 * @param[in]  direction The direction \f$ x \in \mathbb{R}^{n} \f$
 * @param[out] product   The product \f$ y \in \mathbb{R}^{n} \f$
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_qr_mult_orth(SleqpFactQR* fact,
                   const SleqpVec* direction,
                   SleqpVec* product);

/**
 * Compute \f$y = Q^{T} x \f$
 *
 * @param[in]  direction The direction \f$ x \in \mathbb{R}^{n} \f$
 * @param[out] product   The product \f$ y \in \mathbb{R}^{n} \f$
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_qr_mult_orth_trans(SleqpFactQR* fact,
                         const SleqpVec* direction,
                         SleqpVec* product);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_qr_capture(SleqpFactQR* fact);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_qr_release(SleqpFactQR** star);

#ifdef SLEQP_HAVE_QR_FACT

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_qr_create_default(SleqpFactQR** star, SleqpSettings* settings);

#endif

#endif /* SLEQP_FACT_QR_H */
