#include <stdio.h>
#include <gsl/gsl_math.h>
#include <assert.h>

#include "probability.h"
#include "probability_chain.h"
#include "chain.h"
#include "energy.h"

int main(int narg, char *varg[])
{
    protocol_params_t params;
    double energy, period;
    double latency, probability;

    if (narg != 3) {
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

