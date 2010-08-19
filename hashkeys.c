#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>

#include "wildmac.h"
#include "hashkeys.h"

key_t *create_key_protocol_nk(protocol_params_t *p, int n, int k, int negate)
{
    key_t *res = malloc(sizeof(key_t));
    switch (n - k) {
        case -1:
            n = k - 1;
            k = -1;
            break;
        case 0:
            n = k;
            k = 0;
            break;
        default:
            n = 1 + k;
            k = 1;
    }

    res->n = n;
    res->k = k;
    res->negate = negate;
    memcpy(&res->p, p, sizeof(protocol_params_t));
    return res;
}


unsigned int key_hash(void *k)
{
    unsigned int result = 0;
    key_t *key = (key_t *) k;
    result = ((unsigned int) (key->p.tau / M_PI * 100) << 24);
    result |= (key->p.samples << 16);
    result |= (key->n - key->k) << 8;
    result |= key->n;
    if (key->negate)
        result &= 0x7fffffff;
    else
        result |= 0x80000000; 
    return result;
}

int key_equal(void *k1, void *k2)
{
    return memcmp(k1, k2, sizeof(key_t)) == 0;
}


