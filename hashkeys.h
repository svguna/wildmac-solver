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
#ifndef __HASHKEYS_H
#define __HASHKEYS_H

#include "wildmac.h"

struct key {
    int n;
    int k;
    int negate;
    protocol_params_t p;
};
typedef struct key key_t;


key_t *create_key_protocol_nk(protocol_params_t *p, int n, int k, int negate);
unsigned int key_hash(void *k);
int key_equal(void *k1, void *k2);

#endif
