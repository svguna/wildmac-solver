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
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <pthread.h>

#include "integrands.h"
#include "wildmac.h"
#include "probability.h"
#include "hashtable.h"
#include "hashkeys.h"

#define CALLS 500000

double probability_an_bn(protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
   
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, 1, 0); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    double xl[] = {
        p->tau,
        0,
        0
    };
    double xu[] = {
        p->on - p->lambda,
        2 * M_PI - p->on,
        2 * M_PI - p->on
    };
    double res, err;
    gsl_monte_function F ={
        .f = &integrand_n_n,
        .dim = 3,
        .params = p
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_a0_b0(protocol_params_t *p)
{
#ifdef CONTACT_VARIABLE
    double result = (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;
#else
    double result = 1.;
#endif
    return result * probability_an_bn(p);
}


double probability_an_bn1(protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
   
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, 1, 0); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    double xl[] = {
        p->tau,
        0,
        2 * M_PI
    };
    double xu[] = {
        p->on - p->lambda,
        2 * M_PI - p->on,
        4 * M_PI - p->on
    };
    double res, err;
    gsl_monte_function F ={
        .f = &integrand_n_n1,
        .dim = 3,
        .params = p
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_a0_bm1(protocol_params_t *p)
{
#ifdef CONTACT_VARIABLE
    double result = (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;
#else
    double result = 1.;
#endif
    return result * probability_an_bn1(p);
}


double probability_bn_an(protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
   
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, 1, 0); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    double xl[] = {
        p->lambda - p->on,
        0,
        0
    };
    double xu[] = {
        -p->tau,
        2 * M_PI - p->on,
        2 * M_PI - p->on
    };
    double res, err;
    gsl_monte_function F ={
        .f = &integrand_n_n,
        .dim = 3,
        .params = p
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_b0_a0(protocol_params_t *p)
{
#ifdef CONTACT_VARIABLE
    double result = (2 * M_PI + p->on - 2 * p->lambda - 
            (p->lambda + p->on) / (2 * M_PI - p->on)) / 4 / M_PI;
#else
    double result = 1.;
#endif
    return result * probability_bn_an(p);
}


double probability_bn1_an(protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
   
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, 1, 0); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    double xl[] = {
        p->lambda - p->on,
        0,
        2 * M_PI
    };
    double xu[] = {
        -p->tau,
        2 * M_PI - p->on,
        4 * M_PI - p->on
    };
    double res, err;
    gsl_monte_function F ={
        .f = &integrand_n_n1,
        .dim = 3,
        .params = p
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_bm1_a0(protocol_params_t *p)
{
#ifdef CONTACT_VARIABLE
    double result = (2 * M_PI - p->lambda) / (2 * M_PI - p->on) / 2 / M_PI;
#else
    double result = 1.;
#endif
    return result * probability_bn1_an(p);
}


double probability_slotn(protocol_params_t *p)
{
    return probability_an_bn(p) + probability_bn_an(p);
}


double probability_slot0(protocol_params_t *p)
{
    return probability_b0_a0(p) + probability_a0_b0(p);
}


double probability_slotn1(protocol_params_t *p)
{
    return probability_an_bn1(p) + probability_bn1_an(p);
}


double probability_slotm1(protocol_params_t *p)
{
    return probability_bm1_a0(p) + probability_a0_bm1(p);
}

