#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>

#include "wildmac.h"
#include "probability.h"


double probability_an_bn(protocol_params_t *p)
{
    double active = ACTIVE(p);
    return (p->on - p->tau - p->lambda)  *
        (4 * M_PI - p->on + p->tau - p->lambda) / 
        4 / M_PI / active;
}


double probability_a0_b0(protocol_params_t *p)
{
    double result = (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;
    return result * probability_an_bn(p);
}


double probability_bn_an(protocol_params_t *p)
{
    double active = ACTIVE(p);
    return (p->on - p->tau - p->lambda) / 4 / M_PI / active *
        (4 * M_PI - 3 * p->on - p->tau + p->lambda);
}


double probability_b0_a0(protocol_params_t *p)
{
    double result = (2 * M_PI + p->on - 2 * p->lambda - 
            (p->lambda + p->on) / (2 * M_PI - p->on)) / 4 / M_PI;
    return result * probability_bn_an(p);
}


double probability_slotn(protocol_params_t *p)
{
    return (p->samples * p->tau) / M_PI;
}


double probability_slot0(protocol_params_t *p)
{
    return probability_b0_a0(p) + probability_a0_b0(p);
}

