#ifndef __CHAIN_H
#define __CHAIN_H

#include "wildmac.h"

double probability_contact(int n, protocol_params_t *p);
double contact_union(int n, protocol_params_t *p);
double contact_intersect(int n, int s, protocol_params_t *p);

#endif
