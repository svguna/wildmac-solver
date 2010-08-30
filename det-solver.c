/*
 * wildmac-prob-solver - returns the proper configuration of the wildmac 
 * protocol, for a deterministic latency.
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

#include <stdio.h>
#include <math.h>
#include "wildmac.h"

static double tau_min, lambda, lambda1, lambda2;

static double slot_energy(double tau, int s)
{
    double w;
    w = (s + 1) * Irx * lambda;
    w += Irx * (tau + lambda);
    w += Ioff * (2 * M_PI - tau - 2 * lambda - s * lambda);
    return w;
}


static double get_lambda1()
{
    double sq = lambda * (Irx - Ioff) * (Itx - Ioff) /  M_PI;
    return sqrt(sq);
}


static double get_lambda2()
{
    double sq = lambda * (Irx - Ioff) * (Itx - Ioff) / (2 * M_PI - lambda);
    return -sqrt(sq);
}


int assert_constraints(double tau, int s)
{
    printf("%f\n", (s + 1) * tau + lambda);
    if ((s + 1) * tau > M_PI) 
        return 0;
    if ((s + 1) * tau + lambda > 2 * M_PI) 
        return 0;
    if (tau < tau_min)
        return 0;
    return 1;
}


int main()
{
    double tau, st;

    lambda = get_lambda(100000);
    tau_min = 2 * M_PI * MINttx / 100000; 
    lambda1 = get_lambda1();
    lambda2 = get_lambda2();

    printf("%f\n", tau_min);

    tau = lambda2 * (2 * M_PI - lambda) / (Ioff - Itx);
    st = (2 * M_PI - lambda) / tau - 1;
    printf("tau = %f, st = %f\n", tau, st);
    if (assert_constraints(tau, floor(st)))
            printf("\tW(floor) = %f\n", slot_energy(tau, floor(st)));
    if (assert_constraints(tau, ceil(st)))
        printf("\tW(ceil) = %f\n", slot_energy(tau, ceil(st)));

    tau = lambda1 * M_PI / (Itx - Ioff);
    st = M_PI / tau - 1;
    printf("tau = %f, st = %f\n", tau, st);
    if (assert_constraints(tau, floor(st)))
        printf("\tW(floor) = %f\n", slot_energy(tau, floor(st)));
    if (assert_constraints(tau, ceil(st)))
        printf("\tW(ceil) = %f\n", slot_energy(tau, ceil(st)));

    tau = tau_min;
    st = (2 * M_PI - lambda) / tau - 1;
    printf("tau = %f, st = %f\n", tau, st);
    if (assert_constraints(tau, floor(st)))
        printf("\tW(floor) = %f\n", slot_energy(tau, floor(st)));
    if (assert_constraints(tau, ceil(st)))
        printf("\tW(ceil) = %f\n", slot_energy(tau, ceil(st)));

    st = M_PI / tau - 1;
    printf("tau = %f, st = %f\n", tau, st);
    if (assert_constraints(tau, floor(st)))
        printf("\tW(floor) = %f\n", slot_energy(tau, floor(st)));
    if (assert_constraints(tau, ceil(st)))
        printf("\tW(ceil) = %f\n", slot_energy(tau, ceil(st)));

    return 0;
}

