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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <pthread.h>

#include "wildmac.h"
#include "probability.h"
#include "hashtable.h"
#include "hashkeys.h"
#include "integrands.h"

#define CALLS 500000
#define CONSEC5(p) (3 * p->tau * (p->samples + 1) - p->lambda)


static double probability_chain_an(int n, int k, protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;

    int i, j, diff;
    double xl[45], xu[45];
    double res, err;
    chain_params_t chain_params = {
        .n = n,
        .k = k,
        .protocol = p
    };
    gsl_monte_function F = {
        .f = &integrand_chain_an,
        .dim = 2 * k + 1,
        .params = &chain_params
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    assert(k > 0);
    assert(k < 6);

    if (k > n * 2 + 1)
        return 0;

    if (k > 3 && CONSEC5(p) < 2 * M_PI)
        return 0; 
    
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, n, k); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    F.dim = 0;

    for (i = 0; i < k; i++) {
        diff = i / 2 + 1;
        
        for (j = 0; j < i; j++) {
            if (j % 2 == 0) { 
                xl[F.dim] = p->on - 4 * M_PI;
                xu[F.dim++] = 2 * M_PI - p->on;
            } else {
                xl[F.dim] = p->on - 2 * M_PI;
                xu[F.dim++] = 4 * M_PI - p->on;
            }
            
            xl[F.dim] = 2 * (n - diff) * M_PI;
            xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
            
            if ((i + j) % 2 == 0)
                diff--;

            xl[F.dim] = 2 * (n - diff) * M_PI;
            xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
        }

        if (i % 2 == 1) {
            xl[F.dim] = p->tau;
            xu[F.dim++] = p->on - p->lambda;
        } else {
            xl[F.dim] = p->lambda - p->on;
            xu[F.dim++] = -p->tau;
        }
        
        xl[F.dim] = 2 * (n - diff) * M_PI;
        xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
        
        diff--;
        
        xl[F.dim] = 2 * (n - diff) * M_PI;
        xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
    }
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);

    gsl_rng_free(r);
#ifdef CONTACT_VARIABLE
    if (n * 2 + 1 - k == 1)
        res *= (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;
    if (n * 2 + 1 - k == 0)
        res *= (2 * M_PI - p->lambda) / (2 * M_PI - p->on) / 2 / M_PI;
#endif

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


static double probability_chain_bn(int n, int k, protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;

    int i, j, diff;
    double xl[45], xu[45];
    double res, err;
    chain_params_t chain_params = {
        .n = n,
        .k = k,
        .protocol = p
    };
    gsl_monte_function F = {
        .f = &integrand_chain_bn,
        .dim = 2 * k + 1,
        .params = &chain_params
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    assert(k > 0);
    assert(k < 6);

    if (k > 2 * (n + 1))
        return 0;

    if (k > 3 && CONSEC5(p) < 2 * M_PI)
        return 0; 
    
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, n, k); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    F.dim = 0;

    for (i = 0; i < k; i++) {
        diff = (i + 1) / 2;

        for (j = 0; j < i; j++) {
            if (j % 2 == 1) { 
                xl[F.dim] = p->on - 4 * M_PI;
                xu[F.dim++] = 2 * M_PI - p->on;
            } else {
                xl[F.dim] = p->on - 2 * M_PI;
                xu[F.dim++] = 4 * M_PI - p->on;
            }
            
            xl[F.dim] = 2 * (n - diff) * M_PI;
            xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
            
            if ((i + j) % 2 == 1)
                diff--;

            xl[F.dim] = 2 * (n - diff) * M_PI;
            xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
        }

        if (i % 2 == 0) {
            xl[F.dim] = p->tau;
            xu[F.dim++] = p->on - p->lambda;
        } else {
            xl[F.dim] = p->lambda - p->on;
            xu[F.dim++] = -p->tau;
        }

        xl[F.dim] = 2 * (n - diff) * M_PI;
        xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
        
        xl[F.dim] = 2 * (n - diff) * M_PI;
        xu[F.dim++] = 2 * (n + 1 - diff) * M_PI - p->on;
    }

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);

    gsl_rng_free(r);
#ifdef CONTACT_VARIABLE
    if (2 * (n + 1) - k == 1)
        res *= (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;
    if (2 * (n + 1) - k == 0)
        res *= (2 * M_PI - p->lambda) / (2 * M_PI - p->on) / 2 / M_PI;
#endif

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_ank_bn(int n, int k, protocol_params_t *p)
{
    return probability_chain_bn(n, 2 * k + 1, p);
}


double probability_bnk_bn(int n, int k, protocol_params_t *p)
{
    return probability_chain_bn(n, 2 * k, p);
}


double probability_ank_an(int n, int k, protocol_params_t *p)
{
    return probability_chain_an(n, 2 * k, p);
}


double probability_bnk_an(int n, int k, protocol_params_t *p)
{
    return probability_chain_an(n, 2 * k - 1, p);
}

