#ifndef __CHAIN_H
#define __CHAIN_H

#include "wildmac.h"

double probability_contact(int n, protocol_params_t *p);
double union_funcf(int n, protocol_params_t *p);
double intersect_funcf(int n, int s, protocol_params_t *p);

#endif
