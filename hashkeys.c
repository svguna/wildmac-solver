#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>

#include "wildmac.h"
#include "hashkeys.h"

key_t *create_key_protocol_nk(protocol_params_t *p, int n, int k)
{
    key_t *res = malloc(sizeof(key_t));
    if (n - k == 0) {
        n = k;
        k = 0;
    } else {
        n = 1 + k;
        k = 1;
    }
    res->n = n;
    res->k = k;
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
    return result;
}

int key_equal(void *k1, void *k2)
{
    return memcmp(k1, k2, sizeof(key_t)) == 0;
}


