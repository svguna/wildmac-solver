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
#ifndef __PROBABILITY_H
#define __PROBABILITY_H

#include "wildmac.h"

double probability_slotn(protocol_params_t *p);
double probability_an_bn(protocol_params_t *p);
double probability_bn_an(protocol_params_t *p);

double probability_slotn1(protocol_params_t *p);
double probability_an_bn1(protocol_params_t *p);
double probability_bn1_an(protocol_params_t *p);

double probability_slot0(protocol_params_t *p);
double probability_a0_b0(protocol_params_t *p);
double probability_b0_a0(protocol_params_t *p);

double probability_slotm1(protocol_params_t *p);
double probability_a0_bm1(protocol_params_t *p);
double probability_bm1_a0(protocol_params_t *p);

#endif
