/*
 * wildmac-solver - returns the proper configuration of the wildmac protocol,
 * given a desired detection latency and probability.
 * Copyright (C) 2010  Stefan Guna
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
#ifndef __INTEGRANDS_H
#define __INTEGRANDS_H
#include "wildmac.h"

struct chain_params {
    int n;
    int k;
    protocol_params_t *protocol;
};
typedef struct chain_params chain_params_t;

double integrand_n_n(double *x, size_t dim, void *params);
double integrand_n_n1(double *x, size_t dim, void *params);
double integrand_chain_bn(double *x, size_t dim, void *params);
double integrand_chain_an(double *x, size_t dim, void *params);

#endif
