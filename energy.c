#include <assert.h>
#include <gsl/gsl_math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h>

#include "energy.h"
#include "chain.h"
#include "wildmac.h"
#include "pthread_sem.h"

enum {
    NO_SOLUTION = -1,
    TRIVIAL,
    TOL_REACHED,
    MAXCALL_REACHED
};


struct worker_task {
    double lb, ub;
    double T;
    int slot;
    protocol_params_t pc;
};


struct worker_data {
    double probability;
    struct worker_task *task;

    int *finish;

    /* result is outputed in the following */
    double *min_energy;
    protocol_params_t *params;
    double *period;

    pthread_sem_t *sem_new_task;
    pthread_sem_t *sem_task_buffered;
    pthread_sem_t *sem_cpu_available;
    pthread_mutex_t *task_mutex;
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


static void *worker_thread(void *data)
{
    int res;
    struct worker_data *wd = (struct worker_data *) data;
    struct worker_task task;
    double energy;

    pthread_sem_up(1, wd->sem_cpu_available);
    while (!*wd->finish) {
        pthread_mutex_lock(wd->task_mutex);
        pthread_sem_down(1, wd->sem_new_task, wd->task_mutex); 
        
        if (*wd->finish) {
            pthread_mutex_unlock(wd->task_mutex);
            break;
        }

        memcpy(&task, wd->task, sizeof(struct worker_task));

        pthread_sem_up(1, wd->sem_task_buffered);
        pthread_mutex_unlock(wd->task_mutex);

        res = find_optimal(wd->probability, task.lb, task.ub, task.T, task.slot,
                &task.pc, &energy);
        if (res == NO_SOLUTION) {
            pthread_sem_up(1, wd->sem_cpu_available);
            continue;
        }
        
        pthread_mutex_lock(wd->task_mutex);

        if (energy > *wd->min_energy) {
            pthread_sem_up(1, wd->sem_cpu_available);
            pthread_mutex_unlock(wd->task_mutex);
            continue;
        }

        *wd->min_energy = energy;
        memcpy(wd->params, &task.pc, sizeof(protocol_params_t));
        *wd->period = task.T;
        
        pthread_sem_up(1, wd->sem_cpu_available);
        pthread_mutex_unlock(wd->task_mutex);
    }
    return NULL;
}


double get_protocol_parameters(double latency, double probability,
        double *period, protocol_params_t *params)
{
    int i, j, max_slots, max_samples;
    double min_energy = DBL_MAX;
    double lambda;
    pthread_t *threads;
    int thread_num = sysconf(_SC_NPROCESSORS_ONLN);

    int finish = 0;
    
    pthread_sem_t sem_cpu_available;
    pthread_sem_t sem_new_task;
    pthread_sem_t sem_task_buffered;
    pthread_mutex_t task_mutex = PTHREAD_MUTEX_INITIALIZER;

    struct worker_task task;
    struct worker_data worker_data = {
        .probability = probability,
        .task = &task,
        
        .finish = &finish,
        
        .min_energy = &min_energy,
        .params = params,
        .period = period,

        .sem_new_task = &sem_new_task,
        .sem_task_buffered = &sem_task_buffered,
        .sem_cpu_available = &sem_cpu_available,
        .task_mutex = &task_mutex,
    };

    assert(period != NULL);
    assert(params != NULL);
    
    pthread_sem_init(0, &sem_cpu_available);
    pthread_sem_init(0, &sem_new_task);
    pthread_sem_init(0, &sem_task_buffered);

    latency *= 100;
    max_slots = latency / 2 / (2 * MINttx + trx);

    threads = malloc(thread_num * sizeof(pthread_t));
    for (i = 0; i < thread_num; i++)
        pthread_create(threads + i, NULL, worker_thread, &worker_data);

    pthread_mutex_lock(&task_mutex);

    pthread_sem_down(1, &sem_cpu_available, &task_mutex);
    
    for (i = 0; i < max_slots; i++) {
        task.slot = i;
        task.T = latency / (i + 1);
        lambda = get_lambda(task.T);
        task.lb = 2 * lambda;

        if (2 * M_PI * MINttx / task.T > task.lb)
            task.lb = 2 * M_PI * MINttx / task.T;

        max_samples = (M_PI - lambda) / task.lb - 1;

        for (j = 1; j <= max_samples; j++) {
            task.ub = (M_PI - lambda) / (j + 1);
            task.pc.lambda = lambda;
            task.pc.samples = j;

            pthread_sem_up(1, &sem_new_task);
            pthread_sem_down(1, &sem_task_buffered, &task_mutex);
            pthread_sem_down(1, &sem_cpu_available, &task_mutex);
        }
    }
    pthread_mutex_unlock(&task_mutex);
    
    finish = 1;
    pthread_mutex_lock(&task_mutex);
    pthread_sem_up(thread_num, &sem_new_task);
    pthread_mutex_unlock(&task_mutex);

    for (i = 0; i < thread_num; i++)
        pthread_join(threads[i], NULL);
    
    free(threads);
    return min_energy;
}

