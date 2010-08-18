#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "pdf_chain.h"

#define CONSEC5(p) (3 * p->tau * (p->samples + 1) - p->lambda)


static double pdf_an_bn1(double x, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    double tmp = ACTIVE(p->protocol);

    if (p->protocol->on - 4 * M_PI < x && x < p->protocol->on - 2 * M_PI)
        return (x + 4 * M_PI - p->protocol->on) / 2 / M_PI / tmp;

    if (p->protocol->on - 2 * M_PI < x && x < -p->protocol->on)
        return 1 / tmp;
    if (-p->protocol->on < x && x < 2 * M_PI - p->protocol->on)
        return (2 * M_PI - x - p->protocol->on) / 2 / M_PI / tmp;
        
    return 0;
}


static double pdf_an_bn(double x, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    double tmp = ACTIVE(p->protocol);

    if (p->protocol->on - 2 * M_PI < x && x < p->protocol->on)
        return (x + 2 * M_PI - p->protocol->on) / 2 / M_PI / tmp;
    if (p->protocol->on < x && x < 2 * M_PI - p->protocol->on)
        return 1 / tmp;
    if (2 * M_PI - p->protocol->on < x && x < 4 * M_PI - p->protocol->on)
        return (4 * M_PI - x - p->protocol->on) / 2 / M_PI / tmp;

    return 0;
}


static double pdf_an_ank(double x, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    double tmp = ACTIVE(p->protocol);

    if (p->k == 0)
        return 1;

    if (2 * (p->k - 1) * M_PI + p->protocol->on < x &&
            x < 2 * p->k * M_PI)
        return (x + 2 * (1 - p->k) * M_PI - p->protocol->on) / tmp / tmp;
    
    if (2 * p->k * M_PI < x && x < 2 * (p->k + 1) * M_PI - p->protocol->on)
        return (2 * (1 + p->k) * M_PI - x - p->protocol->on) / tmp / tmp;
    
    if (x == 2 * p->k * M_PI)
        return 1 / tmp;

    return 0;
}


static double pdf_an_bnk(double x, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    double tmp = ACTIVE(p->protocol);

    if (p->k == 1)
        return pdf_an_bn1(x, params);

    if (2 * (p->k - 2) * M_PI + p->protocol->on < x &&
            x < 2 * (p->k - 1) * M_PI)
        return (x + 2 * (2 - p->k) * M_PI - p->protocol->on) / 2 / M_PI / tmp;
    
    if (2 * (p->k - 1) * M_PI < x && x < 2 * p->k * M_PI)
        return 1 / 2 / M_PI;

    if (2 * p->k * M_PI < x && x < 2 * (p->k + 1) * M_PI - p->protocol->on)
        return (2 * (p->k + 1) * M_PI - x - p->protocol->on) / 2 / M_PI / tmp;

    return 0;
}



static double pdf_bnk_bn(double x, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    double tmp = ACTIVE(p->protocol);

    if (p->k == 0)
        return 1;

    if (p->protocol->on - 2 * (1 + p->k) * M_PI < x && x < -2 * p->k * M_PI)
        return (x + 2 * (1 + p->k) * M_PI - p->protocol->on) / tmp / tmp;
    
    if (-2 * p->k * M_PI < x && x < 2 * (1 - p->k) * M_PI - p->protocol->on)
        return (2 * (1 - p->k) * M_PI - x - p->protocol->on) / tmp / tmp;

    if (x == 2 * p->k * M_PI)
        return 1 / tmp;
    
    return 0;
}


static double pdf_ank_bn(double x, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    double tmp = ACTIVE(p->protocol);

    if (p->k == 0)
        return pdf_an_bn(x, params);

    if (p->protocol->on - 2 * (2 + p->k) * M_PI < x && 
            x < -2 * (1 + p->k) - p->protocol->on)
        return (x + 4 * M_PI - p->protocol->on) / 2 / M_PI / tmp;

    if (-2 * (1 + p->k) - p->protocol->on < x &&
            x < -2 * p->k * M_PI)
        return 1 / 2 / M_PI;

    if (-2 * p->k * M_PI < x && 2 * (1 - p->k) * M_PI - p->protocol->on)
        return (2 * M_PI - p->protocol->on - x) / 2 / M_PI / tmp;

    return 0;
}


double pdf_chain_ank_an(double *x, size_t dim, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    int i, j;
    double result = 1;
    assert(dim == p->k * 2);
    assert(p->k > 0);

    assert(p->k < 3);
    if (p->k == 2 && CONSEC5(p->protocol) < 2 * M_PI)
        return 0;

    for (i = p->k; i > 0; i--) {
        double y;
        chain_params_t p2 = {
            .n = p->n,
            .k = i,
            .protocol = p->protocol
        };
        
        y = 0;
        for (j = 0; j < i; j++)
            y += x[2 * j + 1] - x[2 * j];
        result *= pdf_an_ank(y, &p2);

        y = 0;
        for (j = 0; j < i - 1; j++)
            y += x[2 * j + 1] - x[2 * j];
        y -= x[2 * i - 2];
        result *= pdf_an_bnk(y, &p2);
    }
    return result;
}


double pdf_chain_bnk_an(double *x, size_t dim, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    int i, j;
    double result = 1;
    assert(dim == p->k * 2 - 1);
    assert(p->k > 0);

    assert(p->k < 4);
    if (p->k == 3 && CONSEC5(p->protocol) < 2 * M_PI)
        return 0;

    for (i = p->k; i > 0; i--) {
        double y;
        chain_params_t p2 = {
            .n = p->n,
            .k = i,
            .protocol = p->protocol
        };
        
        y = 0;
        for (j = 0; j < i - 1; j++)
            y += x[2 * j + 1] - x[2 * j];
        y -= x[2 * i - 2];
        result *= pdf_an_bnk(y, &p2);
    
        p2.k = i - 1;
        y = 0;
        for (j = 0; j < i - 1; j++)
            y += x[2 * j + 1] - x[2 * j];
        result *= pdf_an_ank(y, &p2);
    }

    return result;
}


double pdf_chain_bnk_bn(double *x, size_t dim, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    int i, j;
    double result = 1;
    assert(dim == p->k * 2);
    assert(p->k > 0);

    assert(p->k < 3); // for k > 2 pdf is 0, but we don't use it
    if (p->k == 2 && CONSEC5(p->protocol) < 2 * M_PI)
        return 0;

    for (i = p->k; i > 0; i--) {
        double y;
        chain_params_t p2 = {
            .n = p->n,
            .k = i,
            .protocol = p->protocol
        };
        
        y = 0;
        for (j = 0; j < i; j++)
            y += x[2 * j + 1] - x[2 * j];
        result *= pdf_bnk_bn(y, &p2);

        p2.k = i - 1;
        y = 0;
        for (j = 0; j < i - 1; j++)
            y += x[2 * j + 1] - x[2 * j];
        y -= x[2 * i - 2];
        result *= pdf_ank_bn(y, &p2);
    }
    return result;
}


double pdf_chain_ank_bn(double *x, size_t dim, void *params)
{
    chain_params_t *p = (chain_params_t *) params;
    int i, j;
    double result = 1;
    assert(dim == p->k * 2 + 1);
    assert(p->k >= 0);

    assert(p->k < 3); // for k > 2 pdf is 0, but we don't use it
    
    if (p->k == 2 && CONSEC5(p->protocol) < 2 * M_PI)
        return 0;

    for (i = p->k; i >= 0; i--) {
        double y;
        chain_params_t p2 = {
            .n = p->n,
            .k = i,
            .protocol = p->protocol
        };
        
        y = 0;
        for (j = 0; j < i; j++)
            y += x[2 * j + 1] - x[2 * j];
        y -= x[2 * i];
        result *= pdf_ank_bn(y, &p2);
    
        y = 0;
        for (j = 0; j < i; j++)
            y += x[2 * j + 1] - x[2 * j];
        result *= pdf_bnk_bn(y, &p2);
    }
    return result;
}



