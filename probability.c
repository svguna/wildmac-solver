#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "pdf.h"

double probability_slotn()
{
    double i;
    double lambda = .01 * M_PI;
    int s = 5;
    double tau = .07 * M_PI;
    double on = tau * (s + 1) + lambda;
    double resultAB, resultBA, error;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(WSIZE);
    gsl_function F;

    F.function = &pdf_slotn;
    F.params = NULL;
    
    gsl_integration_qags(&F, tau - lambda, on - lambda, EPSABS, EPSREL, WSIZE, 
            w, &resultAB, &error);
    gsl_integration_workspace_free(w);
    w = gsl_integration_workspace_alloc(WSIZE);
    gsl_integration_qags(&F, lambda - on, lambda - tau, EPSABS, EPSREL, WSIZE, 
            w, &resultBA, &error);
    
    gsl_integration_workspace_free(w);
    return resultAB + resultBA;
}

