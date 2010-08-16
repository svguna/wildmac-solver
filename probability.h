#ifndef __PROBABILITY_H
#define __PROBABILITY_H

#include "wildmac.h"

double probability_slotn(protocol_params_t *p);
double probability_an_bn(protocol_params_t *p);
double probability_bn_an(protocol_params_t *p);

double probability_slot0(protocol_params_t *p);
double probability_a0_b0(protocol_params_t *p);
double probability_b0_a0(protocol_params_t *p);

#endif
