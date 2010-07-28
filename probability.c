#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "pdf.h"

double probability_slotn()
{
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


static double slot0_integrand1(double x, void *p)
{
    func_params_t *params = (func_params_t *) p;
    assert(params->type == OFFSET);

    return pdf_after_contact(params->data.offset - x, NULL) * 
        pdf_slotn(x, NULL);
}


static double slot0_integrand2(double x, void *p)
{
    func_params_t integrand_params;
    func_params_t *params = (func_params_t *) p;
    gsl_function F;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(WSIZE);
    assert(params->type == INTERVAL);
    double result, error;

    fill_offset(&integrand_params, x);

    F.function = &slot0_integrand1;
    F.params = &integrand_params;

    gsl_integration_qags(&F, params->data.domain.begin, params->data.domain.end,
            EPSABS, EPSREL, WSIZE, w, &result, &error);
    
    gsl_integration_workspace_free(w);
    return 0;
}


double probability_slot0()
{
    double lambda = .01 * M_PI;
    int s = 5;
    double tau = .07 * M_PI;
    double on = tau * (s + 1) + lambda;
    double resultAB, resultBA, error;
    func_params_t limits;

    gsl_integration_workspace *w;
    gsl_function F;

    F.function = &slot0_integrand2;
    F.params = &limits;

    w = gsl_integration_workspace_alloc(WSIZE);
    fill_interval(&limits, tau - lambda, on - lambda);
    gsl_integration_qags(&F, 0, 2 * M_PI, EPSABS, EPSREL, WSIZE, w, &resultAB, 
            &error);
    gsl_integration_workspace_free(w);
    
    w = gsl_integration_workspace_alloc(WSIZE);
    fill_interval(&limits, lambda - on, lambda - tau);
    gsl_integration_qags(&F, 0, 2 * M_PI, EPSABS, EPSREL, WSIZE, w, &resultBA, 
            &error);
    gsl_integration_workspace_free(w);
    return resultAB + resultBA;
}

