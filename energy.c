#include <assert.h>
#include <gsl/gsl_math.h>
#include <stdio.h>
#include <string.h>

#include "energy.h"
#include "chain.h"


enum {
    NO_SOLUTION = -1,
    TRIVIAL,
    TOL_REACHED,
    MAXCALL_REACHED
};


static inline double slot_energy(double tau, double T, int samples)
{
    double res = 0;

    res += tau * T / 2 / M_PI * Itx;
    res += trx *(Itx + samples * Irx);
    res += (T - tau * T / 2 / M_PI - (samples + 1) * trx) * Ioff;
    res /= T;
    return res;
}


static int find_optimal(double prob_bound, double lb, double ub, double T, 
        int slot, protocol_params_t *params, double *energy)
{
    unsigned long calls;
    double middle = (ub - lb) / 2 + lb;
    double last_energy, new_energy;
    
    assert(energy != NULL);
    assert(params != NULL);
    
    *energy = DBL_MAX;

    params->tau = ub;
    SET_ON(params);

    if (contact_union(slot, params) < prob_bound)
        return NO_SOLUTION;
    last_energy = slot_energy(params->tau, T, params->samples);

    params->tau = lb; 
    SET_ON(params);
    
    if (contact_union(slot, params) > prob_bound) {
        *energy = last_energy;
        return TRIVIAL;
    }

    for (calls = 0; calls < MAX_CALLS; calls++) {
        double prob;
        
        params->tau = middle;
        SET_ON(params);
        
        prob = contact_union(slot, params);

        if (prob >= prob_bound) { 
            double delta;
            
            new_energy = slot_energy(params->tau, T, params->samples);
            delta = fabs(new_energy - last_energy);
            if (delta / last_energy < TOL_REL) {
                *energy = new_energy;
                return TOL_REACHED;
            }
            last_energy = new_energy;

            ub = middle;
        }
        else 
            lb = middle;
        
        middle = (ub - lb) / 2 + lb;
    }
    *energy = new_energy;
    return MAXCALL_REACHED;
}


double get_protocol_parameters(double latency, double probability,
        double *period, protocol_params_t *params)
{
    int i, max_slots; 
    double min_energy = DBL_MAX;

    assert(period != NULL);
    assert(params != NULL);

    latency *= 100;
    max_slots = latency / 2 / (2 * MINttx + trx);

    for (i = 0; i < max_slots; i++) {
        double T = latency / (i + 1);
        double lambda = get_lambda(T);
        int j, max_samples;
        double tx_lb = 2 * lambda;

        if (2 * M_PI * MINttx / T > tx_lb)
            tx_lb = 2 * M_PI * MINttx / T;

        max_samples = (M_PI - lambda) / tx_lb - 1;

        for (j = 1; j <= max_samples; j++) {
            double lb = tx_lb;
            double ub = (M_PI - lambda) / (j + 1);
            double energy;
            protocol_params_t pc = {
                .lambda = lambda,
                .samples = j
            };
            int res;
            
            res = find_optimal(probability, lb, ub, T, i, &pc, &energy);
            if (res == NO_SOLUTION || energy > min_energy)
                continue;

            min_energy = energy;
            memcpy(params, &pc, sizeof(protocol_params_t));
            *period = T;
        }
    }

    return min_energy;
}

