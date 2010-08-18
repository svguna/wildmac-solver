#include <stdio.h>
#include <gsl/gsl_math.h>

#include "wildmac.h"
#include "probability.h"
#include "probability_chain.h"

static double union_funcg(int n, protocol_params_t *p);

double union_funcf(int n, protocol_params_t *p)
{
    double r = 0;
    int i;

    if (n < 0)
        return 0;
    if (n == 0)
        return probability_slot0(p);

    r += union_funcg(n, p);

    if (n == 0)
        r += probability_b0_a0(p);
    else
        r += probability_bn_an(p) * (1 - union_funcf(n - 1, p)); 
    
    for (i = 0; i <= 2; i++) 
        r += (1 - union_funcf(n - i - 1, p)) * probability_ank_bn(n, i, p);
    
    for (i = 1; i <= 2; i++) 
        r += (union_funcg(n - i, p) - 1) * probability_bnk_bn(n, i, p);

    return r;
}


static double union_funcg(int n, protocol_params_t *p)
{
    double r = 0;
    int i;

    if (n <= 0)
        return 0;

    r += union_funcf(n - 1, p);

    r += probability_bn_an(p) * (1 - union_funcg(n - 1, p));

    for (i = 1; i <= 2; i++)
        r += probability_ank_an(n, i, p)  * (union_funcf(n - i - 1, p) - 1);
    
    for (i = 1; i <= 3; i++) 
        r += probability_bnk_an(n, i, p) * (1 - union_funcg(n - i, p));
    return r;
}


