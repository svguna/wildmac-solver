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
#ifndef __WILDMAC_H
#define __WILDMAC_H

struct protocol_params {
    double tau;
    double lambda;
    int samples;
    double on; // solely to improve performance
    double active; // solely to improve performance
    double active2; // solely to improve performance

};
typedef struct protocol_params protocol_params_t;

// time values expressed in tens of us
#define Itx 1652

#define Iup 908
#define Tup 445

#define Idown 950
#define Tdown 110

#define Isample 977
#define Tsample 406

#define Ioff 9
#define trx 400

#define MAXLATENCY 360000000

extern int battery;
extern int min_ttx;

#define get_lambda(T) (trx * M_PI * 2 / T)
#define SET_ACTIVE(p) (p)->active = 2 * M_PI - (p)->on;\
                      (p)->active2 = (p)->active * (p)->active
#define SET_ON(p) (p)->on = ((p)->samples + 1) * (p)->tau + (p)->lambda

#define HASH_SIZE 16

#define MAX_EVAL 100
#define MAX_CALLS 10
#define TOL_REL 1e-2

#endif

