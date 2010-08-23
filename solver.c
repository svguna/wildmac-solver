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
#include <assert.h>

#include "probability.h"
#include "probability_chain.h"
#include "chain.h"
#include "energy.h"


static void print_boilerplate()
{
    printf("wildmac-solver  Copyright (C) 2010  Stefan Guna\n"
           "This program comes with ABSOLUTELY NO WARRANTY.\n"
           "This is free software, and you are welcome to redistribute it\n"
           "under certain conditions; see the accompanying license file\n"
           "for details.\n\n");
}


int main(int narg, char *varg[])
{
    protocol_params_t params;
    double energy, period;
    double latency, probability;

    if (narg != 3) {
        print_boilerplate();
        printf("Invalid arguments. Please run the solver as follows:\n\n"
                "\t%s LATENCY PROBABILITY\n\n"
                "where:\n"
                "\tLATENCY must be provided in ms, and\n"
                "\tPROBABILITY is a float between 0 and 1.\n\n",
                varg[0]);
        return -1;
    }
    sscanf(varg[1], "%lf", &latency);
    assert(latency * 100 > (MINttx * 2 + Irx) * 2);
    
    sscanf(varg[2], "%lf", &probability);
    assert(probability < 1);
    assert(probability > 0);
    
    energy = get_protocol_parameters(latency, probability, &period, &params);

    print_boilerplate();
    if (energy == DBL_MAX) {
        printf("No suitable configuration found.\n");
        return -1;
    }

    printf("energy per time: %f\n", energy);
    printf("         period: %.2fms\n", period);
    printf("         beacon: %.2fms\n", period * params.tau / 2 / M_PI);
    printf("        samples: %d\n\n", params.samples);
    
    return 0;
}

