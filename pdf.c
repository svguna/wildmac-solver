#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "pdf.h"


struct conv_int_params {
    double x;
    struct conv_params *params;
};


void fill_domain_pdf_beacon(func_params_t *params, double beacon, double cca, 
        int cca_samples, int slot)
{
    params->type = SIMPLE;
    params->data.domain.begin = 2 * M_PI * slot;
    params->data.domain.end  = 2 * M_PI * (slot + 1) - 
        beacon * (cca_samples + 1) - cca;
}


void fill_domain_pdf_phase(func_params_t *params, int slot)
{
    params->type = SIMPLE;
    params->data.domain.begin = 2 * M_PI * slot;
    params->data.domain.end = 2 * M_PI * (slot + 1);
}


static void fill_composed_params(int type, func_params_t *params,
        basic_function pdf1, void *params1, basic_function pdf2, void *params2)
{
    params->type = type;
    params->data.conv.pdf1.function = pdf1;
    params->data.conv.pdf2.function = pdf2;
    params->data.conv.pdf1.params = params1;
    params->data.conv.pdf2.params = params2;
}


void fill_composed_sum(func_params_t *params, basic_function pdf1,
        void *params1, basic_function pdf2, void *params2)
{
    fill_composed_params(SUM, params, pdf1, params1, pdf2, params2);
}


void fill_composed_difference(func_params_t *params, basic_function pdf1,
        void *params1, basic_function pdf2, void *params2)
{
    fill_composed_params(DIFFERENCE, params, pdf1, params1, pdf2, params2);
}


double pdf_uniform(double x, void *p)
{
    func_params_t *params = (func_params_t *) p;
    assert(params->type == SIMPLE);
    
    if (params->data.domain.begin <= x && x <= params->data.domain.end) 
        return 1. / (params->data.domain.end - params->data.domain.begin);
    return 0; 
}


static double convolution_integrand(double x, void *p)
{
    struct conv_int_params *int_params = (struct conv_int_params *) p;
    double x_out = int_params->x;
    struct conv_params *conv_params = int_params->params;
    gsl_function *pdf1 = &conv_params->pdf1;
    gsl_function *pdf2 = &conv_params->pdf2;
    return pdf1->function(x_out - x, pdf1->params) * 
        pdf2->function(x, pdf2->params);     
}


static double deconvolution_integrand(double x, void *p)
{
    struct conv_int_params *int_params = (struct conv_int_params *) p;
    double x_out = int_params->x;
    struct conv_params *conv_params = int_params->params;
    gsl_function *pdf1 = &conv_params->pdf1;
    gsl_function *pdf2 = &conv_params->pdf2;

    return pdf1->function(x_out + x, pdf1->params) * 
        pdf2->function(x, pdf2->params);     
}


// TODO recursive formula, switch to iterative
static void fill_domain(struct domain *domain, func_params_t *params)
{
    struct domain d1, d2;
    switch (params->type) {
        case SIMPLE:
            domain->begin = params->data.domain.begin;
            domain->end= params->data.domain.end;
            return;

        case SUM:
            fill_domain(&d1, params->data.conv.pdf1.params);
            fill_domain(&d2, params->data.conv.pdf2.params);
            domain->begin = d1.begin + d2.begin;
            domain->end = d1.end + d2.end;
            return;

        case DIFFERENCE:
            fill_domain(&d1, params->data.conv.pdf1.params);
            fill_domain(&d2, params->data.conv.pdf2.params);
            domain->begin = d1.begin - d2.end;
            domain->end = d1.end - d2.begin;
            return;
    }
    assert(0);
}


double pdf_composed(double x, void *p)
{
    func_params_t *params = (func_params_t *) p;
    assert(params->type == SUM || params->type == DIFFERENCE);
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(WSIZE);
    struct domain limits;
    
    gsl_function F;
    struct conv_int_params int_params;
    double result, error;
    size_t neval;
    static size_t repeat = 0;

    fill_domain(&limits, params->data.conv.pdf2.params);

    int_params.x = x;
    int_params.params = &params->data.conv;

    if (params->type == SUM) 
        F.function = &convolution_integrand;
    else 
        F.function = &deconvolution_integrand;
    F.params = &int_params;

    
    gsl_integration_qags(&F, limits.begin, limits.end, EPSABS, EPSREL, WSIZE, w,
            &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}


double pdf_slotn(double x, void *params)
{
    func_params_t dom_beacon, dom_phase;
    func_params_t sum_params, diff_params;

    fill_domain_pdf_beacon(&dom_beacon, .07 * M_PI, .01 * M_PI, 5, 0);
    fill_domain_pdf_phase(&dom_phase, 0);

    fill_composed_sum(&sum_params, pdf_uniform, &dom_phase, pdf_uniform,
            &dom_beacon);
    
    fill_composed_difference(&diff_params, pdf_composed, &sum_params, 
            pdf_uniform, &dom_beacon);
    return pdf_composed(x, &diff_params);
}

