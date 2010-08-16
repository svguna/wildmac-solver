#ifndef __PDF_CHAIN_H
#define __PDF_CHAIN_H

#include <stdlib.h>

#include "wildmac.h"

struct chain_params {
    int n;
    int k;
    protocol_params_t *protocol;
};
typedef struct chain_params chain_params_t;


double pdf_an_bn1(double x, void *params);
double pdf_an_bn(double x, void *params);
double pdf_chain_ank_an(double *x, size_t dim, void *params);
double pdf_chain_bnk_an(double *x, size_t dim, void *params);
double pdf_chain_bnk_bn(double *x, size_t dim, void *params);
double pdf_chain_ank_bn(double *x, size_t dim, void *params);



#endif
