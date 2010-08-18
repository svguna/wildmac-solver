#ifndef __PROBABILITY_CHAIN
#define __PROBABILITY_CHAIN

double probability_ank_bn(int n, int k, protocol_params_t *p);
double probability_bnk_bn(int n, int k, protocol_params_t *p);
double probability_ank_an(int n, int k, protocol_params_t *p);
double probability_bnk_an(int n, int k, protocol_params_t *p);

#endif
