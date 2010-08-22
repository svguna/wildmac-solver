#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <assert.h>

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
        res -= probability_bnk_bn(0, 1, 0, p);
        return res;
    }
    res += probability_slotn(p);
    res += probability_slotn1(p);
    res -= probability_bnk_bn(n, 1, 0, p);
    return res;
}


double contact_union(int n, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i;

//    printf("slot %d, protocol params: tau = %f, lambda = %f, samples = %d\n",
//            n, p->tau, p->lambda, p->samples);
    if (n < 0) 
        return 0;
    
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, n, 0); 
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
        r += (1 - contact_union(n - i - 1, p)) * probability_ank_bn(n, i, 0, p);
    
    for (i = 1; i <= 2; i++) 
        r += (union_funcg(n - i, p) - 1) * probability_bnk_bn(n, i, 0, p);
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}


static double union_funcg(int n, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i;

    if (n < 0) 
        return 0;
    
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, n, 0); 
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
        r += probability_ank_an(n, i, 0, p) * (contact_union(n - i - 1, p) - 1);
    
    for (i = 1; i <= 3; i++) 
        r += probability_bnk_an(n, i, 0, p) * (1 - union_funcg(n - i, p));
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}


double contact_intersect(int n, int s, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i, sign;

    if (n < s) 
        return 0;

    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, n, 0); 
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
        r += sign * probability_ank_bn(n, i, 0, p) * 
            contact_intersect(n - i - 1, s, p);
    
    for (i = 1, sign = -1; n - i >= s - 1 && i <= 2; i++, sign *= -1) 
        r += sign * probability_bnk_bn(n, i, 0, p) * 
            intersect_funcg(n - i, s, p);
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}


static double intersect_funcg(int n, int s, protocol_params_t *p)
{
    static struct hashtable *hash_table = NULL;
    key_t *hash_key;
    double *hash_res;
    
    double r = 0;
    int i, sign;

    if (n < s - 1) 
        return 0;
    if (n == s - 1)
        return 1;
    
    if (hash_table == NULL)
        hash_table = create_hashtable(16, key_hash, key_equal);
    
    hash_key = create_key_protocol_nk(p, n, 0, n); 
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
        r += probability_ank_an(n, i, 0, p) * 
            contact_intersect(n - i - 1, s, p);

    for (i = 2, sign = -1; n - i >= s - 1&& i <= 3; i++, sign *= -1)
        r += probability_bnk_an(n, i, 0, p) * intersect_funcg(n - i, s, p);
    
    hash_res = malloc(sizeof(double));
    *hash_res = r;
    hashtable_insert(hash_table, hash_key, hash_res);

    return r;
}

