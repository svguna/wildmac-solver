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

// curent values expressed in one tenth of a mA, i.e., 1569 = 15.69 mA
// time values expressed in one tenth of a ms, i.e., 1478 = 14.78 mS

// current consumed while transmitting
#define Itx 1652

// average current consumed while ramping up
#define Iup 908
// duration of ramp up
#define Tup 445

// average current consumed while ramping down
#define Idown 950
// duration of ramp down
#define Tdown 110

// average current consumed while sampling the channel
#define Isample 1569
// duration of sample
#define Tsample 1478

// current consumed when the radio if off
#define Ioff 9

// minimum overlap between beacon and reception to register a contact 
// required for 4 11B packets and 3 .48ms gaps
#define trx 285

// maximum latency for a lifetime problem
#define MAXLATENCY 360000000

extern int battery;
extern int min_ttx;

#define get_lambda(T) (trx * M_PI * 2 / T)
#define SET_ACTIVE(p) (p)->active = 2 * M_PI - (p)->on;\
                      (p)->active2 = (p)->active * (p)->active
#define SET_ON(p) (p)->on = ((p)->samples + 1) * (p)->tau + (p)->lambda

// maximum number of calls to energy function
#define MAX_CALLS 10
// accuracy of the computation
#define TOL_REL 1e-2

#endif

