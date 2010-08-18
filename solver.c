#include <gsl/gsl_math.h>
#include <stdio.h>

#include "probability.h"
#include "probability_chain.h"
#include "chain.h"

int main()
{
    protocol_params_t params = {
        .tau = .1 * M_PI,
        .lambda = .05 * M_PI,
        .samples = 6
    };
    SET_ON(&params);

    printf("slot n: %f\n", probability_slotn(&params));
    printf("An <- Bn: %f\n", probability_an_bn(&params));
    printf("Bn <- An: %f\n", probability_bn_an(&params));
   
    printf("\n");
    printf("slot 0: %f\n", probability_slot0(&params));
    printf("A0 <- B0: %f\n", probability_a0_b0(&params));
    printf("B0 <- A0: %f\n", probability_b0_a0(&params));

    int n = 1;
    printf("\n");
    printf("An-0 <- Bn: %f\n", probability_ank_bn(n, 0, &params));
    printf("An-1 <- Bn: %f\n", probability_ank_bn(n, 1, &params));
    printf("An-2 <- Bn: %f\n", probability_ank_bn(n, 2, &params));
    
    printf("\n");
    printf("Bn-1 <- Bn: %f\n", probability_bnk_bn(n, 1, &params));
    printf("Bn-2 <- Bn: %f\n", probability_bnk_bn(n, 2, &params));
    
    printf("\n");
    printf("An-1 <- An: %f\n", probability_ank_an(n, 1, &params));
    printf("An-2 <- An: %f\n", probability_ank_an(n, 2, &params));

    
    printf("\n");
    printf("Bn-1 <- An: %f\n", probability_bnk_an(n, 1, &params));
    printf("Bn-2 <- An: %f\n", probability_bnk_an(n, 2, &params));
    printf("Bn-3 <- An: %f\n", probability_bnk_an(n, 3, &params));

    printf("\n");
    
    printf("union_f(0): %f\n", union_funcf(0, &params));
    printf("union_f(1): %f\n", union_funcf(1, &params));
    printf("union_f(2): %f\n", union_funcf(2, &params));
    printf("union_f(3): %f\n", union_funcf(3, &params));
    printf("union_f(4): %f\n", union_funcf(4, &params));
    printf("union_f(5): %f\n", union_funcf(5, &params));
    printf("union_f(6): %f\n", union_funcf(6, &params));
    printf("union_f(7): %f\n", union_funcf(7, &params));
    return 0;
}
