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
#ifndef __SOLVER_H
#define __SOLVER_H

#include "wildmac.h"

// Currents are given in mA.
// Time is given in tens of us.

double energy(double period, double tau, int samples);
double get_latency_params(double latency, double probability, double *period, 
        protocol_params_t *params);
double get_lifetime_params(double lifetime, double probability, double *period,
        protocol_params_t *params);

#endif

