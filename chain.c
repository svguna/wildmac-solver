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
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <pthread.h>

#include "wildmac.h"
#include "probability.h"
#include "probability_chain.h"
#include "hashtable.h"
#include "hashkeys.h"

static double union_funcg(int n, protocol_params_t *p);
static double intersect_funcg(int n, int s, protocol_params_t *p);


double probability_contact(int n, protocol_params_t *p)
{
    double res = 0;
    if (n == 0) {
        res += probability_slot0(p);
        res += probability_slotm1(p);
        res -= probability_bnk_bn(0, 1, p);
        return res;
    }
    res += probability_slotn(p);
    res += probability_slotn1(p);
    res -= probability_bnk_bn(n, 1, p);
    return res;
}


double contact_union(int n, protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i;

    if (n < 0) 
        return 0;
    
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, n, n); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    r += union_funcg(n, p);

    if (n == 0)
        r += probability_b0_a0(p);
    else
        r += probability_bn_an(p) * (1 - contact_union(n - 1, p)); 
    
    for (i = 0; i <= 2; i++) 
        r += (1 - contact_union(n - i - 1, p)) * probability_ank_bn(n, i, p);
    
    for (i = 1; i <= 2; i++) 
        r += (union_funcg(n - i, p) - 1) * probability_bnk_bn(n, i, p);
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}


static double union_funcg(int n, protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i;

    if (n < 0) 
        return 0;
    
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, n, n); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    r += contact_union(n - 1, p);

    if (n == 0)
        r += probability_a0_bm1(p);
    else
        r += probability_an_bn1(p) * (1 - union_funcg(n - 1, p));

    for (i = 1; i <= 2; i++)
        r += probability_ank_an(n, i, p) * (contact_union(n - i - 1, p) - 1);
    
    for (i = 1; i <= 3; i++) 
        r += probability_bnk_an(n, i, p) * (1 - union_funcg(n - i, p));
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}


double contact_intersect(int n, int s, protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i, sign;

    if (n < s) 
        return 0;

    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, n, n); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    if (n == 0)
        r += probability_slot0(p);
    else
        if (n == s)
            r += probability_slotn(p);
        else
            r += contact_intersect(n - 1, s, p) * probability_slotn(p);
  
    r += intersect_funcg(n, s, p);
    
    for (i = 1, sign = -1; n - i >= s && i <= 2; i++, sign *= -1) 
        r += sign * probability_ank_bn(n, i, p) * 
            contact_intersect(n - i - 1, s, p);
    
    for (i = 1, sign = -1; n - i >= s - 1 && i <= 2; i++, sign *= -1) 
        r += sign * probability_bnk_bn(n, i, p) * 
            intersect_funcg(n - i, s, p);
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}


static double intersect_funcg(int n, int s, protocol_params_t *p)
{
    static pthread_mutex_t hash_mutex = PTHREAD_MUTEX_INITIALIZER;
    static struct hashtable *hash_table = NULL;
    hashkey_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i, sign;

    if (n < s - 1) 
        return 0;
    if (n == s - 1)
        return 1;
    
    pthread_mutex_lock(&hash_mutex);
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal, &hash_mutex);
    pthread_mutex_unlock(&hash_mutex);
    
    hash_key = create_key_protocol_nk(p, n, n); 
    hash_res = hashtable_search(hash_table, hash_key);
    if (hash_res != NULL) {
        free(hash_key);
        return *hash_res;
    }

    if (n == 0)
        r += probability_slotm1(p); 
    else
        r += intersect_funcg(n - 1, s, p) * probability_slotn1(p);

    for (i = 1, sign = 1; n - i >= s && i <= 2; i++, sign *= -1)
        r += probability_ank_an(n, i, p) * contact_intersect(n - i - 1, s, p);

    for (i = 2, sign = -1; n - i >= s - 1&& i <= 3; i++, sign *= -1)
        r += probability_bnk_an(n, i, p) * intersect_funcg(n - i, s, p);
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}

