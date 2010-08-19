#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>

#include "wildmac.h"
#include "probability.h"
#include "pdf_chain.h"
#include "hashtable.h"
#include "hashkeys.h"

#define CALLS 500000


double probability_ank_bn(int n, int k, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;

    int i;
    double xl[7], xu[7];
    double res, err;
    chain_params_t chain_params = {
        .n = n,
        .k = k,
        .protocol = p
    };
    gsl_monte_function F = {
        .f = &pdf_chain_ank_bn,
        .dim = 2 * k + 1,
        .params = &chain_params
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    assert(k >= 0);
    assert(k < 3);

    if (n - k < 0)
        return 0;
    
    if (k == 0) {
        if (n == 0)
            return probability_a0_b0(p);
        return probability_an_bn(p);
    }

    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, k); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    for (i = 0; i < k; i++) {
        xl[2 * i] = p->tau;
        xu[2 * i] = p->on - p->lambda;
        xl[2 * i + 1] = p->lambda - p->on;
        xu[2 * i + 1] = -p->tau;
    }
    xl[2 * k] = p->tau;
    xu[2 * k] = p->on - p->lambda;
   
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);
    gsl_rng_free(r);
    
    if (n - k == 0)
        res *= (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_bnk_bn(int n, int k, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;

    int i;
    double xl[6], xu[6];
    double res, err;
    chain_params_t chain_params = {
        .n = n,
        .k = k,
        .protocol = p
    };
    gsl_monte_function F = {
        .f = &pdf_chain_bnk_bn,
        .dim = 2 * k,
        .params = &chain_params
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    assert(k > 0);
    assert(k < 3);

    if (n - k < -1) 
        return 0;

    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, k); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }
    
    for (i = 0; i < k; i++) {
        xl[2 * i] = p->tau;
        xu[2 * i] = p->on - p->lambda;
        xl[2 * i + 1] = p->lambda - p->on;
        xu[2 * i + 1] = -p->tau;
    }
   
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);
    gsl_rng_free(r);
  
    if (n - k == -1) 
        res *= (2 * M_PI + p->on - 2 * p->lambda -
                (p->lambda + p->on) / (2 * M_PI - p->on)) / 4 / M_PI;

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_ank_an(int n, int k, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;

    int i;
    double xl[6], xu[6];
    double res, err;
    chain_params_t chain_params = {
        .n = n,
        .k = k,
        .protocol = p
    };
    gsl_monte_function F = {
        .f = &pdf_chain_ank_an,
        .dim = 2 * k,
        .params = &chain_params
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    assert(k > 0);
    assert(k < 3);

    if (n - k < 0) 
        return 0;

    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, k); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    for (i = 0; i < k; i++) {
        xl[2 * i] = p->lambda - p->on;
        xu[2 * i] = -p->tau;
        xl[2 * i + 1] = p->tau;
        xu[2 * i + 1] = p->on - p->lambda;
    }
   
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);
    gsl_rng_free(r);

    if (n - k == 0)
        res *= (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;

    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}


double probability_bnk_an(int n, int k, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;

    int i;
    double xl[7], xu[7];
    double res, err;
    chain_params_t chain_params = {
        .n = n,
        .k = k,
        .protocol = p
    };
    gsl_monte_function F = {
        .f = &pdf_chain_bnk_an,
        .dim = 2 * k - 1,
        .params = &chain_params
    };
    gsl_monte_plain_state *s;
    const gsl_rng_type *T;
    gsl_rng *r;

    assert(k > 0);
    assert(k < 4);

    if (n - k < -1) 
        return 0; 

    if (k == 1) {
        if (n == 0)
            return probability_bm1_a0(p);
        return probability_bn1_an(p);
    }

    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, k); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    for (i = 0; i < k - 1; i++) {
        xl[2 * i] = p->lambda - p->on;
        xu[2 * i] = -p->tau;
        xl[2 * i + 1] = p->tau;
        xu[2 * i + 1] = p->on - p->lambda;
    }
    xl[2 * k - 2] = p->lambda - p->on;
    xu[2 * k - 2] = -p->tau;
   
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    s = gsl_monte_plain_alloc(F.dim);
    gsl_monte_plain_integrate(&F, xl, xu, F.dim, CALLS, r, s, &res, &err);
    gsl_monte_plain_free(s);
    gsl_rng_free(r);
    
    if (n - k == -1) 
        res *= (2 * M_PI + p->on - 2 * p->lambda -
                (p->lambda + p->on) / (2 * M_PI - p->on)) / 4 / M_PI;
    
    hash_res = malloc(sizeof(double));
    *hash_res = res;
    hashtable_insert(hash_table, hash_key, hash_res);

    return res;
}



