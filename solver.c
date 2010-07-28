#include <gsl/gsl_math.h>
#include <stdio.h>

#include "pdf.h"
#include "probability.h"

int main()
{
    printf("slot n: %f\n", probability_slotn());
    printf("slot 0: %f\n", probability_slot0());
    return 0;
}
