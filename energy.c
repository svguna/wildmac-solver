#include <assert.h>
#include <nlopt.h>
#include <gsl/gsl_math.h>
#include <stdio.h>
#include <string.h>

#include "energy.h"
#include "chain.h"

enum {
    PROBABILITY = 0,
    CONSTRAINTS_NO
};


struct opt_data {
    int slot_count;
    double lambda;
    double T;
    int samples;
};
typedef struct opt_data opt_data_t;


struct constraint_data {
    int type;
    double bound;
    opt_data_t *opt_data;
};
typedef struct constraint_data constraint_data_t;


static inline double slot_energy(double tau, double T, int samples)
{
    double res = 0;

    res += tau * T / 2 / M_PI * Itx;
    res += trx *(Itx + samples * Irx);
    res += (T - tau * T / 2 / M_PI - (samples + 1) * trx) * Ioff;
    printf(" %f\n", res); 
    return res;
}


static double obj_energy(unsigned n, const double *x, double *grad, void *data)
{
    opt_data_t *opt_data = (opt_data_t *) data;
    assert(n == 1);

    if (grad) {
        grad[0] = (Itx - Ioff) / 2 / M_PI;
    }

    printf("tau=%f T=%.2f n=%d samples=%d", x[0], opt_data->T / 100, 
            opt_data->slot_count, opt_data->samples);
    return slot_energy(x[0], opt_data->T, opt_data->samples) / opt_data->T;
}


static double constraint(unsigned n, const double *x, double *grad, void *data)
{
    constraint_data_t *constraint_data = (constraint_data_t *) data;
    protocol_params_t params = {
        .tau = x[0],
        .lambda = constraint_data->opt_data->lambda,
        .samples = constraint_data->opt_data->samples
    };
    SET_ON(&params);
    int slot_count = constraint_data->opt_data->slot_count;

    assert(constraint_data->type >= 0 || 
            constraint_data->type < CONSTRAINTS_NO);
    assert(grad == NULL);
    assert(n == 1);
 
    printf("slot %d, tau %f, lambda %f, samples %d, duration %f\n",
            slot_count, params.tau, params.lambda, params.samples, x[2]);

    switch (constraint_data->type) {
        case PROBABILITY:
            printf("\tprobability: %f\n", constraint_data->bound - contact_union(slot_count, &params));
            return constraint_data->bound - contact_union(slot_count, &params);
    }

    return 0;
}


double get_protocol_parameters(double latency, double probability,
        double *period, protocol_params_t *params)
{
    int i, max_slots; 
    
    latency *= 100;
    max_slots = latency / 6 / trx;

    printf("%d\n", max_slots);
    
    for (i = 1; i < max_slots; i++) {
        double T = latency / (i + 1);
        double lambda = get_lambda(T);
        int j, max_samples;
        double tx_lb = 2 * lambda;

        if (2 * M_PI * MINttx / T > tx_lb)
            tx_lb = 2 * M_PI * MINttx / T;

        max_samples = M_PI / lambda - 2;

        for (j = 1; j <= max_samples; j++) {
            double lb[1] = {tx_lb};
            double ub[1] = {(M_PI - lambda) / (j + 1)};
            double x[1] = {(ub[0] - lb[0]) / 2 + lb[0]};
            int opt_result;
            nlopt_opt opt;
            double minf;
            
            opt_data_t opt_data = {
                .slot_count = i,
                .lambda = lambda,
                .T = T,
                .samples = j
            };

            constraint_data_t constraint_data[1] = {
                {.type = PROBABILITY, 
                 .bound = probability, 
                 .opt_data = &opt_data}
            };

//            opt = nlopt_create(NLOPT_GN_ISRES, 1);
            opt = nlopt_create(NLOPT_LD_MMA, 1);
            nlopt_set_lower_bounds(opt, lb);
            nlopt_set_upper_bounds(opt, ub);
            nlopt_set_min_objective(opt, obj_energy, &opt_data);
            nlopt_add_inequality_constraint(opt, constraint, constraint_data, 
                    1e-8);
            nlopt_set_ftol_rel(opt, 1e-2);
            nlopt_set_maxeval(opt, MAX_EVAL);

            printf("%d %f %d %d %f-%f\n", i + 1, T, max_samples, j, lb[0], ub[0]);
            printf("%d\n", nlopt_optimize(opt, x, &minf));
            nlopt_destroy(opt);
            return -1;
        }
    }

    
    return -1;





    printf("done\n");
}

