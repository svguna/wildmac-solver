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
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>

#include "wildmac.h"
#include "hashkeys.h"

hashkey_t *create_key_protocol_nk(protocol_params_t *p, int n, int k)
{
    hashkey_t *res = malloc(sizeof(hashkey_t));
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
    memcpy(&res->p, p, sizeof(protocol_params_t));
    return res;
}


unsigned int key_hash(void *k)
{
    unsigned int result = 0;
    hashkey_t *key = (hashkey_t *) k;
    result = ((unsigned int) (key->p.tau / M_PI * 100) << 24);
    result |= (key->p.samples << 16);
    result |= (key->n - key->k) << 8;
    result |= key->n;
    return result;
}

int key_equal(void *k1, void *k2)
{
    return memcmp(k1, k2, sizeof(hashkey_t)) == 0;
}


