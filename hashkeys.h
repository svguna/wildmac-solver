#ifndef __HASHKEYS_H
#define __HASHKEYS_H

#include "wildmac.h"

struct key {
    int n;
    int k;
    protocol_params_t p;
};
typedef struct key key_t;


key_t *create_key_protocol_nk(protocol_params_t *p, int n, int k);
unsigned int key_hash(void *k);
int key_equal(void *k1, void *k2);

#endif
