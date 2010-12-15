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

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "common-prints.h"
#include "wildmac.h"


static double tau_min, lambda, period;
static double w_max, lifetime;
int battery;
int min_ttx;

static double energy_per_time(double tau, int s)
{
    double w = 0;
    double beacon = tau * period / 2 / M_PI;
    
    w += Isample * Tsample * s;
    w += Iup * Tup + Idown * Tdown;
    w += Itx * beacon;
    w += Ioff * (period - s * Tsample - Tdown - Tup - beacon);
    if (w < 0)
        return DBL_MAX;
    return w / period;
}


static double latency_params(double *tau, int *s)
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

        w = energy_per_time(t, i);

        if (w == DBL_MAX || w > w_min)
            continue;
        w_min = w;
        *tau = t;
        *s = i;
    }
    return w_min;
}


static int check_args(int narg, char *varg[])
{
    if (narg == 5 && strlen(varg[1]) == 1 && (varg[1][0] == 'l' || 
            varg[1][0] == 'e'))
        return 0;

    print_boilerplate();
    printf("Invalid arguments. Please run the solver as follows:\n\n"
            "\t%s (l LATENCY) | (e LIFETIME) BATTERY MIN_TTX\n\n"
            "where:\n"
            "\t `l' gives the best configuration to meet the latency "
            "requirements\n"
            "\t     (LATENCY must be provided in ms)\n"
            "\t `e' gives the best configuration to meet the lifetime "
            "requirements\n"
            "\t     (LIFETIME must be provided in hours)\n"
            "\t BATTERY capacity is expressed in mAh\n"
            "\t MIN_TTX minimum beacon duration (in us).\n\n",
            varg[0]);
    return 1;
}


static void solve_latency()
{
    double tau = 0, w;
    int s = 0;

    period *= 100; 

    lambda = get_lambda(period);
    tau_min = 2 * M_PI * min_ttx / period; 

    print_boilerplate();
    
    w = latency_params(&tau, &s);
    if (w == DBL_MAX) {
        printf("No suitable configuration found.\n");
        return;
    }
    period /= 100;

    printf("For the desired latency of %.2f ms, battery capacity of %d mAh, "
            "and minimum beacon duration of %d us, use the following "
            "configuration:\n", period, battery / 100, min_ttx / 100);
    printf("      lifetime: %.2f h\n", battery / w);
    printf("   avg current: %f mA\n", w / 100);
    printf("        period: %.2f ms\n", period);
    printf("        beacon: %.2f ms\n", period * tau / 2 / M_PI);
    printf("       samples: %d\n\n", s);
} 


static void solve_lifetime()
{
    double lb, ub, middle;
    double last_latency;
    unsigned long calls;
    double tau;
    int s;
    double w;

    w_max = battery / lifetime;

    lb = 4 * min_ttx;
    ub = MAXLATENCY;
    middle = (ub - lb) / 2 + lb;
 
    print_boilerplate();
    
    period = last_latency = ub;
    lambda = get_lambda(ub);
    tau_min = 2 * M_PI * min_ttx / ub;
    w = latency_params(&tau, &s);

    if (w > w_max) {
        printf("No suitable configuration found.\n");
        return;
    }

    for (calls = 0; calls < MAX_CALLS * 100; calls++) {
        lambda = get_lambda(middle);
        period = middle;
        tau_min = 2 * M_PI * min_ttx / middle;
        w = latency_params(&tau, &s);

        if (w <= w_max) {
            double delta;

            delta = fabs(middle - last_latency);
            if (delta / last_latency < TOL_REL) {
                break;
            }
            last_latency = ub = middle;
        } else {
            lb = middle;
        }

        middle = (ub - lb) / 2 + lb;
    }
    middle /= 100;

    printf("For the desired lifetime of %.2f h, battery capacity of %d mAh, "
            "and minimum beacon duration of %d us, use the following "
            "configuration:\n", lifetime, battery / 100, min_ttx / 100);
    printf("        lifetime: %.2f h\n", battery / w);
    printf("latency & period: %.2f ms\n", middle);
    printf("     avg current: %f mA\n", w / 100);
    printf("          beacon: %.2f ms\n", middle * tau / 2 / M_PI);
    printf("         samples: %d\n\n", s);
}


int main(int narg, char *varg[])
{
    if (check_args(narg, varg))
        return -1;

    sscanf(varg[3], "%d", &battery);
    assert(battery > 10);
    assert(battery < 10000);
    battery *= 100;

    sscanf(varg[4], "%d", &min_ttx);
    assert(min_ttx >= 20);
    min_ttx *= 100;

    switch(varg[1][0]) {
        case 'l':
            sscanf(varg[2], "%lf", &period);
            assert(period * 100 > (min_ttx * 2 + trx) * 2);
            solve_latency();
            break;
        case 'e':
            sscanf(varg[2], "%lf", &lifetime);
            assert(lifetime >= 1.);
            solve_lifetime();
            break;
    }
   
    return 0;
}

