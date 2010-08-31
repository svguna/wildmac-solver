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

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <assert.h>

#include "common-prints.h"
#include "wildmac.h"


static double tau_min, lambda, period;


static double energy_per_time(double tau, int s)
{
    double w = 0;
    
    w += Itx * (tau + lambda);
    w += (s + 1) * Irx * lambda;
    w += Ioff * (2 * M_PI - tau - 2 * lambda - s * lambda);
    return w / 2 / M_PI;
}


static double get_protocol_parameters(double *tau, int *s)
{
    int smax, i;
    double w_min = DBL_MAX;

    assert(tau != NULL);
    assert(s != NULL);

    smax = floor((2 * M_PI - 2 * lambda) / tau_min - 1);
   
    for (i = 1; i <= smax; i++) {
        double w, t;
        t = M_PI / (i + 1);
        if (t < tau_min)
            t = tau_min;

        if ((i + 1) * t + lambda > 2 * M_PI - lambda)
            continue;
        w = energy_per_time(t, i);

        if (w > w_min)
            continue;
        w_min = w;
        *tau = t;
        *s = i;
    }
    return w_min;
}


int main(int narg, char *varg[])
{
    double tau = 0, w;
    int s = 0;

    if (narg != 2) {
        print_boilerplate();
        printf("Invalid arguments. Please run the solver as follows:\n\n"
                "\t%s PERIOD\n\n"
                "where:\n"
                "\tPERIOD must be provided in ms.\n\n",
                varg[0]);
        return -1;
    }
    sscanf(varg[1], "%lf", &period);
    assert(period * 100 > (MINttx * 2 + trx) * 2);
    period *= 100; 

    lambda = get_lambda(period);
    tau_min = 2 * M_PI * MINttx / period; 

    print_boilerplate();
    
    w = get_protocol_parameters(&tau, &s);
    if (w == DBL_MAX) {
        printf("No suitable configuration found.\n");
        return -1;
    }
    period /= 100;

    printf("energy per sec: %f\n", w);
    printf("        period: %.2fms\n", period);
    printf("        beacon: %.2fms\n", period * tau / 2 / M_PI + trx / 100.);
    printf("    CCA period: %.2fms\n", period * tau / 2 / M_PI);
    printf("       samples: %d\n\n", s);
    
    return 0;
}

