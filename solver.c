#include <gsl/gsl_math.h>
#include <stdio.h>

#include "probability.h"

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
    return 0;
}
