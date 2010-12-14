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
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "common-prints.h"
#include "probability.h"
#include "probability_chain.h"
#include "chain.h"
#include "solver.h"

int battery;
int min_ttx;

static void solve_latency(double latency, double probability)
{
    protocol_params_t params;
    double energy, period;
    
    print_boilerplate();

    energy = get_latency_params(latency, probability, &period, &params);

    if (energy == DBL_MAX) {
        printf("No suitable configuration found.\n");
        return;
    }

    period /= 100;

    printf("For the desired latency of %.2f ms, battery capacity of %d mAh, "
            "minimum beacon duration of %d us, and minimum probability %.2f "
            "use the following configuration:\n", period, battery / 100, 
            min_ttx / 100, probability);
    printf("      lifetime: %.2f h\n", battery / energy);
    printf("   avg current: %.2f mA\n", energy / 100);
    printf("        period: %.2f ms\n", period);
    printf("        beacon: %.2f ms\n", period * params.tau / 2 / M_PI);
    printf("       samples: %d\n\n", params.samples);
}


static void solve_lifetime(double lifetime, double probability)
{
    protocol_params_t params;
    double latency, period;
    
    print_boilerplate();

    latency = get_lifetime_params(lifetime, probability, &period, &params);

    if (latency == DBL_MAX) {
        printf("No suitable configuration found.\n");
        return;
    }

    latency /= 100;

    printf("For the desired lifetime of %.2f h, battery capacity of %d mAh, "
            "minimum beacon duration of %d us, and minimum probability %.2f, "
            "use the following configuration:\n", lifetime, battery / 100,
            min_ttx / 100, probability);
    printf("      lifetime: %.2f h\n", 
            battery / energy(period, params.tau, params.samples));
    printf("   avg current: %.2f mA\n", 
            energy(period, params.tau, params.samples) / 100);
    
    period /= 100;
    
    printf("       latency: %.2f ms\n", latency);
    printf("        period: %.2f ms\n", period);
    printf("        beacon: %.2f ms\n", period * params.tau / 2 / M_PI);
    printf("       samples: %d\n\n", params.samples);

}

static int check_args(int narg, char *varg[])
{
    if (narg == 6 && strlen(varg[1]) == 1 && (varg[1][0] == 'l' || 
            varg[1][0] == 'e'))
        return 0;

    print_boilerplate();
    printf("Invalid arguments. Please run the solver as follows:\n\n"
            "\t%s (l LATENCY) | (e LIFETIME) PROBABILITY BATTERY MIN_TTX\n\n"
            "where:\n"
            "\t `l' gives the best configuration to meet the latency "
            "requirements\n"
            "\t     (LATENCY must be provided in ms)\n"
            "\t `e' gives the best configuration to meet the lifetime "
            "requirements\n"
            "\t     (LIFETIME must be provided in hours)\n"
            "\t BATTERY capacity is expressed in mAh.\n"
            "\t MIN_TTX minimum beacon duration (in us).\n\n",
            varg[0]);
    return 1;
}



int main(int narg, char *varg[])
{
    double latency, probability, lifetime;

    if (check_args(narg, varg))
        return -1;

    sscanf(varg[3], "%lf", &probability);
    assert(probability < 1);
    assert(probability > 0);

    sscanf(varg[4], "%d", &battery);
    assert(battery > 10);
    assert(battery < 10000);
    battery *= 100;
    
    sscanf(varg[5], "%d", &min_ttx);
    assert(min_ttx >= 20);
    min_ttx *= 100;

    switch(varg[1][0]) {
        case 'l':
            sscanf(varg[2], "%lf", &latency);
            assert(latency * 100 > (min_ttx * 2 + trx) * 2);

            solve_latency(latency, probability);
            break;
        case 'e':
            sscanf(varg[2], "%lf", &lifetime);
            assert(lifetime >= 1.);
            
            solve_lifetime(lifetime, probability);
            break;
    }
    
    return 0;
}

