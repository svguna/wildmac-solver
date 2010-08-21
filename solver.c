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

    printf("\n");
    printf("slot n n1: %f\n", probability_slotn1(&params));
    printf("An <- Bn1: %f\n", probability_an_bn1(&params));
    printf("Bn1 <- An: %f\n", probability_bn1_an(&params));

    printf("\n");
    printf("slot 0 -1: %f\n", probability_slotm1(&params));
    printf("A0 <- B-1: %f\n", probability_a0_bm1(&params));
    printf("B-1 <- A0: %f\n", probability_bm1_a0(&params));

    int n = 5;
    printf("\n");
    printf("An-0 <- Bn: %f\n", probability_ank_bn(n, 0, 0, &params));
    printf("An-1 <- Bn: %f\n", probability_ank_bn(n, 1, 0, &params));
    printf("An-2 <- Bn: %f\n", probability_ank_bn(n, 2, 0, &params));
    
    printf("\n");
    printf("An-1 <- !Bn: %f\n", probability_ank_bn(n, 1, 1, &params));
    printf("An-2 <- !Bn: %f\n", probability_ank_bn(n, 2, 1, &params));
    printf("An-3 <- !Bn: %f\n", probability_ank_bn(n, 3, 1, &params));
    
    printf("\n");
    printf("Bn-2 <- !Bn: %f\n", probability_bnk_bn(n, 2, 1, &params));
    printf("Bn-3 <- !Bn: %f\n", probability_bnk_bn(n, 3, 1, &params));
    printf("Bn-4 <- !Bn: %f\n", probability_bnk_bn(n, 3, 1, &params));
    
    printf("\n");
    printf("Bn-1 <- Bn: %f\n", probability_bnk_bn(n, 1, 0, &params));
    printf("Bn-2 <- Bn: %f\n", probability_bnk_bn(n, 2, 0, &params));
    
    printf("\n");
    printf("An-1 <- An: %f\n", probability_ank_an(n, 1, 0, &params));
    printf("An-2 <- An: %f\n", probability_ank_an(n, 2, 0, &params));
    
    printf("\n");
    printf("An-2 <- !An: %f\n", probability_ank_an(n, 2, 1, &params));
    printf("An-3 <- !An: %f\n", probability_ank_an(n, 3, 1, &params));

    printf("\n");
    printf("Bn-1 <- An: %f\n", probability_bnk_an(n, 1, 0, &params));
    printf("Bn-2 <- An: %f\n", probability_bnk_an(n, 2, 0, &params));
    printf("Bn-3 <- An: %f\n", probability_bnk_an(n, 3, 0, &params));

    printf("\n");
    printf("Bn-2 <- !An: %f\n", probability_bnk_an(n, 2, 1, &params));
    printf("Bn-3 <- !An: %f\n", probability_bnk_an(n, 3, 1, &params));
    printf("Bn-4 <- !An: %f\n", probability_bnk_an(n, 4, 1, &params));

    printf("\n");

    printf("contact(0): %f\n", probability_contact(0, &params));
    printf("contact(1): %f\n", probability_contact(1, &params));
    printf("contact(2): %f\n", probability_contact(2, &params));
    
    printf("\n");
    
    for (n = 0; n < 0; n++) {
        printf("    union_f(%d): %f\n", n, union_funcf(n, &params));
        printf("intersect_f(%d): %f\n", n, intersect_funcf(n, 0, &params));
    }
    return 0;
}
