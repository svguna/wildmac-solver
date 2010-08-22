#include <gsl/gsl_math.h>
#include <stdio.h>

#include "probability.h"
#include "probability_chain.h"
#include "chain.h"
#include "energy.h"

int main()
{
    int n;
    protocol_params_t params;
    double energy, period;
    
    
    energy = get_protocol_parameters(2000., .90, &period, &params);

    printf("energy per time: %f\n", energy);
    printf("         period: %.2fms\n", period);
    printf("         beacon: %.2fms\n", period * params.tau / 2 / M_PI);
    printf("        samples: %d\n", params.samples);
    
    return 0;
}

