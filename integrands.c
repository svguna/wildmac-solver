/*
 * wildmac-solver - returns the proper configuration of the wildmac protocol,
 * given a desired detection latency and probability.
 * Copyright (C) 2010  Stefan Guna, http://disi.unitn.it/~guna
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see 
 * http://www.gnu.org/licenses/gpl-3.0-standalone.html.
 */

#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "integrands.h"

static double pdf_uniform(double x, double begin, double end)
{
    if (x < begin || x > end)
        return 0;
    return 1. / (end - begin);
}


double integrand_n_n(double *x, size_t dim, void *params)
{
    protocol_params_t *p = (protocol_params_t *) params;

    return pdf_uniform(x[0] - x[1] + x[2], 0, 2 * M_PI) *
        pdf_uniform(x[1], 0, 2 * M_PI - p->on) *
        pdf_uniform(x[2], 0, 2 * M_PI - p->on);
}


double integrand_n_n1(double *x, size_t dim, void *params)
{
    protocol_params_t *p = (protocol_params_t *) params;

    return pdf_uniform(x[0] - x[1] + x[2], 0, 2 * M_PI) *
        pdf_uniform(x[1], 0, 2 * M_PI - p->on) *
        pdf_uniform(x[2], 2 * M_PI, 4 * M_PI - p->on);
}


inline double pdfx_n_n(double *x, int n, void *params)
{
    protocol_params_t *p = (protocol_params_t *) params;

    return pdf_uniform(x[0] - x[1] + x[2], 0, 2 * M_PI) *
        pdf_uniform(x[1], 2 * n * M_PI, 2 * (n + 1) * M_PI - p->on) *
        pdf_uniform(x[2], 2 * n * M_PI, 2 * (n + 1) * M_PI - p->on);
}


inline double pdfx_n_n1(double *x, int n, void *params)
{
    protocol_params_t *p = (protocol_params_t *) params;

    return pdf_uniform(x[0] - x[1] + x[2], 0, 2 * M_PI) *
        pdf_uniform(x[1], 2 * (n - 1) * M_PI, 2 * n * M_PI - p->on) *
        pdf_uniform(x[2], 2 * n * M_PI, 2 * (n + 1) * M_PI - p->on);
}


double integrand_chain_bn(double *x, size_t dim, void *params)
{
    chain_params_t *cp = (chain_params_t *) params;
    protocol_params_t *p = cp->protocol;
    double res, xt[3];

    res = pdfx_n_n(x, cp->n, p);
    if (cp->k == 1)
        return res;

    res *= pdfx_n_n1(x + 3, cp->n, p);
    xt[0] = x[0] - x[6] + x[3];
    xt[1] = x[7];
    xt[2] = x[8];
    res *= pdfx_n_n(xt, cp->n, p);
    if (cp->k == 2)
        return res;

    res *= pdfx_n_n(x + 9, cp->n - 1, p);
    res *= pdfx_n_n1(x + 12, cp->n, p);
    xt[0] = x[0] - x[6] + x[15] + x[12] - x[9];
    xt[1] = x[16];
    xt[2] = x[17];
    res *= pdfx_n_n(xt, cp->n, p);
    if (cp->k == 3)
        return res;

    res *= pdfx_n_n1(x + 18, cp->n - 1, p);
    res *= pdfx_n_n(x + 21, cp->n - 1, p);
    res *= pdfx_n_n1(x + 24, cp->n, p);
    xt[0] = x[0] - x[6] + x[15] - x[27] + x[24] - x[21] + x[18];
    xt[1] = x[28];
    xt[2] = x[29];
    res *= pdfx_n_n(xt, cp->n, p);
    if (cp->k == 4)
        return res;

    res *= pdfx_n_n(x + 30, cp->n - 2, p);
    res *= pdfx_n_n1(x + 33, cp->n - 1, p);
    res *= pdfx_n_n(x + 36, cp->n - 1, p);
    res *= pdfx_n_n1(x + 39, cp->n, p);
    xt[0] = x[0] - x[6] + x[15] - x[27] + x[42] + x[39] - x[36] + x[33] - x[30];
    xt[1] = x[43];
    xt[2] = x[44];
    res *= pdfx_n_n(xt, cp->n, p);
    return res;
}


double integrand_chain_an(double *x, size_t dim, void *params)
{
    chain_params_t *cp = (chain_params_t *) params;
    protocol_params_t *p = cp->protocol;
    double res, xt[3];

    res = pdfx_n_n1(x, cp->n, p);
    if (cp->k == 1) 
        return res;
    
    res *= pdfx_n_n(x + 3, cp->n - 1, p);
    xt[0] = x[0] - x[6] + x[3];
    xt[1] = x[7];
    xt[2] = x[8];
    res *= pdfx_n_n1(xt, cp->n, p);
    if (cp->k == 2)
        return res;

    res *= pdfx_n_n1(x + 9, cp->n - 1, p);
    res *= pdfx_n_n(x + 12, cp->n - 1, p);
    xt[0] = x[0] - x[6] + x[15] + x[12] - x[9];
    xt[1] = x[16];
    xt[2] = x[17];
    res *= pdfx_n_n1(xt, cp->n, p);
    if (cp->k == 3)
        return res;

    res *= pdfx_n_n(x + 18, cp->n - 2, p);
    res *= pdfx_n_n1(x + 21, cp->n - 1, p);
    res *= pdfx_n_n(x + 24, cp->n - 1, p);
    xt[0] = x[0] - x[6] + x[15] - x[27] + x[24] - x[21] + x[18];
    xt[1] = x[28];
    xt[2] = x[29];
    res *= pdfx_n_n1(xt, cp->n, p);
    if (cp->k == 4)
        return res;

    res *= pdfx_n_n1(x + 30, cp->n - 2, p);
    res *= pdfx_n_n(x + 33, cp->n - 2, p);
    res *= pdfx_n_n1(x + 36, cp->n - 1, p);
    res *= pdfx_n_n(x + 39, cp->n - 1, p);
    xt[0] = x[0] - x[6] + x[15] - x[27] + x[42] + x[39] - x[36] + x[33] - x[30];
    xt[1] = x[43];
    xt[2] = x[44];
    res *= pdfx_n_n1(xt, cp->n, p);
    return res;
}

