#ifndef __PDF_H
#define __PDF_H

#include <gsl/gsl_math.h>

#define EPSABS 0
#define EPSREL 1e-3
#define WSIZE 5000 

typedef double (* basic_function) (double x, void * params);

enum {
    SIMPLE,
    CONVOLUTION,
    DECONVOLUTION
};


struct domain {
    double begin;
    double end;
};


struct conv_params {
    gsl_function pdf1;
    gsl_function pdf2;
};


struct func_params {
    int type;
    union {
        struct domain domain;
        struct conv_params conv;
    } data;
};
typedef struct func_params func_params_t;


void fill_domain_pdf_beacon(func_params_t *params, double beacon, double cca, 
        int cca_samples, int slot);
void fill_domain_pdf_phase(func_params_t *params, int slot);

void fill_convolution_params(func_params_t *params, basic_function pdf1,
        void *params1, basic_function pdf2, void *params2);
void fill_deconvolution_params(func_params_t *params, basic_function pdf1,
        void *params1, basic_function pdf2, void *params2);

double pdf_uniform(double x, void *p);
double convolution(double x, void *p);
double pdf_slotn(double x, void *params);
#endif
