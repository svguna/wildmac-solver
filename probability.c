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
#include <assert.h>
#include <gsl/gsl_math.h>

#include "wildmac.h"
#include "probability.h"


double probability_an_bn(protocol_params_t *p)
{
    return (p->on - p->tau - p->lambda)  *
        (4 * M_PI - p->on + p->tau - p->lambda) / 
        4 / M_PI / p->active;
}


double probability_a0_b0(protocol_params_t *p)
{
    double result = (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;
    return result * probability_an_bn(p);
}


double probability_an_bn1(protocol_params_t *p)
{
    return probability_bn_an(p);
}


double probability_a0_bm1(protocol_params_t *p)
{
    double result = (2 * M_PI + p->on - 2 * p->lambda) / 4 / M_PI;
    return result * probability_bn_an(p);
}


double probability_bn_an(protocol_params_t *p)
{
    return (p->on - p->tau - p->lambda) / 4 / M_PI / p->active *
        (4 * M_PI - 3 * p->on - p->tau + p->lambda);
}


double probability_b0_a0(protocol_params_t *p)
{
    double result = (2 * M_PI + p->on - 2 * p->lambda - 
            (p->lambda + p->on) / (2 * M_PI - p->on)) / 4 / M_PI;
    return result * probability_bn_an(p);
}


double probability_bn1_an(protocol_params_t *p)
{
    return probability_an_bn(p);
}


double probability_bm1_a0(protocol_params_t *p)
{
    double result = (2 * M_PI + p->on - 2 * p->lambda - 
            (p->lambda + p->on) / (2 * M_PI - p->on)) / 4 / M_PI;
    return result * probability_an_bn(p);
}


double probability_slotn(protocol_params_t *p)
{
    return (p->samples * p->tau) / M_PI;
}


double probability_slot0(protocol_params_t *p)
{
    return probability_b0_a0(p) + probability_a0_b0(p);
}


double probability_slotn1(protocol_params_t *p)
{
    return probability_slotn(p);
}


double probability_slotm1(protocol_params_t *p)
{
    return probability_bm1_a0(p) + probability_a0_bm1(p);
}

