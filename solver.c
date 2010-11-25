/*
 * wildmac-solver - returns the proper configuration of the wildmac protocol,
 * given a desired detection latency and probability.
 * Copyright (C) 2010  Stefan Guna
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see 
 * http://www.gnu.org/licenses/gpl-3.0-standalone.html.
 */
#include <sys/time.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h>

#include "solver.h"
#include "chain.h"
#include "wildmac.h"
#include "pthread_sem.h"


static int thread_cnt = 0;


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
    double *energy;
    int *slots;
    protocol_params_t *params;
    double *period;

    pthread_sem_t *sem_new_task;
    pthread_sem_t *sem_task_buffered;
    pthread_sem_t *sem_worker_available;
    pthread_mutex_t *task_mutex;
};


static inline double energy_per_time(double period, double tau, int samples)
{
    double w = 0;
    double beacon = tau * period / 2 / M_PI;
    
    w += Isample * Tsample * samples;
    w += Iup * Tup + Idown * Tdown;
    w += Itx * beacon;
    w += Ioff * (period - samples * Tsample - Tdown - Tup - beacon);
    return w / period;
}


double energy(double period, double tau, int samples)
{
    return energy_per_time(period, tau, samples);
}


static int find_optimal(double prob_bound, double lb, double ub, double T, 
        int slot, protocol_params_t *params, double *energy)
{
    unsigned long calls;
    double middle = (ub - lb) / 2 + lb;
    double last_energy, new_energy = DBL_MAX;
    
    assert(energy != NULL);
    assert(params != NULL);
    
    *energy = DBL_MAX;

    params->tau = ub;
    SET_ON(params);
    SET_ACTIVE(params);

    if (contact_union(slot, params) < prob_bound)
        return NO_SOLUTION;
    last_energy = energy_per_time(T, params->tau, params->samples);

    params->tau = lb; 
    SET_ON(params);
    SET_ACTIVE(params);
    
    if (contact_union(slot, params) > prob_bound) {
        *energy = last_energy;
        return TRIVIAL;
    }

    for (calls = 0; calls < MAX_CALLS; calls++) {
        double prob;
        
        params->tau = middle;
        SET_ON(params);
        SET_ACTIVE(params);
        
        prob = contact_union(slot, params);

        if (prob >= prob_bound) { 
            double delta;
            
            new_energy = energy_per_time(T, params->tau, params->samples);
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
    int thread_id;
    
    pthread_mutex_lock(wd->task_mutex);
    thread_cnt++;
    thread_id = thread_cnt;
    pthread_mutex_unlock(wd->task_mutex);

    pthread_sem_up(1, wd->sem_worker_available);
    printf("[%d] online\n", thread_id);
    
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
            printf("[%d] finished %dx%.2fms samples=%d no solution\n", 
                    thread_id, task.slot + 1, task.T / 100, task.pc.samples); 
            pthread_sem_up(1, wd->sem_worker_available);
            continue;
        }
        
        printf("[%d] finished %dx%.2fms samples=%d tau=%.2fms I=%.2f mA\n", 
                thread_id, task.slot + 1, task.T / 100, task.pc.samples,
                task.pc.tau * task.T / 100 / 2 / M_PI, energy / 100); 
        pthread_mutex_lock(wd->task_mutex);

        if (energy < *wd->energy) {
            if (wd->slots != NULL) {
                if (*wd->period != DBL_MAX && 
                        (task.slot + 1) * task.T > *wd->period * *wd->slots) {
                    pthread_sem_up(1, wd->sem_worker_available);
                    pthread_mutex_unlock(wd->task_mutex);
                    continue;
                }
                *wd->slots = task.slot + 1;
            } else
                *wd->energy = energy;
            
            memcpy(wd->params, &task.pc, sizeof(protocol_params_t));
            *wd->period = task.T;
        }

        pthread_sem_up(1, wd->sem_worker_available);
        pthread_mutex_unlock(wd->task_mutex);
    }
    printf("[%d] offline\n", thread_id);
    return NULL;
}


unsigned long time_delta(struct timeval *start, struct timeval *end)
{
    double t1, t2;

    t1 = (double) start->tv_sec * 1000 + (double) start->tv_usec / 1000;
    t2 = (double) end->tv_sec * 1000 + (double) end->tv_usec / 1000;
    return t2 - t1;
}


double get_latency_params(double latency, double probability, double *period, 
        protocol_params_t *params)
{
    struct timeval start, end;
    struct timezone tz;
    unsigned long total_states = 0, states_completed = 0;
    unsigned long elapsed, estimated;

    int i, j, max_slots, max_samples;
    double min_energy = DBL_MAX;
    double lambda;
    pthread_t *threads;
    int thread_num = sysconf(_SC_NPROCESSORS_ONLN);

    int finish = 0;
    
    pthread_sem_t sem_worker_available;
    pthread_sem_t sem_new_task;
    pthread_sem_t sem_task_buffered;
    pthread_mutex_t task_mutex = PTHREAD_MUTEX_INITIALIZER;

    struct worker_task task;
    struct worker_data worker_data = {
        .probability = probability,
        .task = &task,
        
        .finish = &finish,
        
        .energy = &min_energy,
        .slots = NULL,
        .params = params,
        .period = period,

        .sem_new_task = &sem_new_task,
        .sem_task_buffered = &sem_task_buffered,
        .sem_worker_available = &sem_worker_available,
        .task_mutex = &task_mutex,
    };

    assert(period != NULL);
    assert(params != NULL);
    
    pthread_sem_init(0, &sem_worker_available);
    pthread_sem_init(0, &sem_new_task);
    pthread_sem_init(0, &sem_task_buffered);

    latency *= 100;
    max_slots = latency / 2 / (2 * MINttx + trx);

    printf("running on %d threads\n", thread_num);
    threads = malloc(thread_num * sizeof(pthread_t));
    for (i = 0; i < thread_num; i++)
        pthread_create(threads + i, NULL, worker_thread, &worker_data);


    for (i = 0; i < max_slots; i++) {
        task.slot = i;
        task.T = latency / (i + 1);
        lambda = get_lambda(task.T);
        task.lb = 2 * lambda;

        if (2 * M_PI * MINttx / task.T > task.lb)
            task.lb = 2 * M_PI * MINttx / task.T;

        total_states += (M_PI - lambda) / task.lb - 1;
    }


    gettimeofday(&start, &tz);
    pthread_mutex_lock(&task_mutex);
    pthread_sem_down(1, &sem_worker_available, &task_mutex);
    
    for (i = 0; i < max_slots; i++) {
        task.slot = i;
        task.T = latency / (i + 1);
        lambda = get_lambda(task.T);
        task.lb = 2 * lambda;

        if (2 * M_PI * MINttx / task.T > task.lb)
            task.lb = 2 * M_PI * MINttx / task.T;

        max_samples = (M_PI - lambda) / task.lb - 1;

        if (energy_per_time(task.T, task.lb, 1) > min_energy) {
            printf("stopping at %d periods, as min(I)=%.2f mA from now\n", 
                    i + 1, energy_per_time(task.T, task.lb, 1) / 100);
            break;
        }

        for (j = 1; j <= max_samples; j++) {
            task.ub = (M_PI - lambda) / (j + 1);
            task.pc.lambda = lambda;
            task.pc.samples = j;

            if (energy_per_time(task.T, task.lb, j) > min_energy) {
                states_completed += max_samples - j + 1;
                printf("stopping samples at %d, as min(I)=%.2f mA from now\n",
                        j, energy_per_time(task.T, task.lb, j) / 100);
                break;
            }

            pthread_sem_up(1, &sem_new_task);
            pthread_sem_down(1, &sem_task_buffered, &task_mutex);
            pthread_sem_down(1, &sem_worker_available, &task_mutex);
            
            states_completed++;
            gettimeofday(&end, &tz);
            elapsed = time_delta(&start, &end);
            estimated = elapsed * total_states / states_completed;
            printf("exploring at %6.2f%% remaining %lds\n", 
                    states_completed * 100. / total_states,
                    (estimated - elapsed) / 1000);
        }
    }
    pthread_mutex_unlock(&task_mutex);
    gettimeofday(&end, &tz);
    elapsed = time_delta(&start, &end);

    finish = 1;
    pthread_mutex_lock(&task_mutex);
    pthread_sem_up(thread_num, &sem_new_task);
    pthread_mutex_unlock(&task_mutex);

    printf("waiting for all workers\n");
    for (i = 0; i < thread_num; i++)
        pthread_join(threads[i], NULL);
    
    printf("\nexplored a total of %ld states in %ld.%lds\n\n", total_states,
            elapsed / 1000, elapsed % 1000);
    
    free(threads);
    return min_energy;
}


static double try_latency(int thread_num, double latency, double max_energy, 
        double probability, struct worker_data *wd)
{
    int i, j, max_slots, max_samples;
    double lambda;
    struct worker_task *task = wd->task;

    printf("trying latency %.2f ms\n", latency / 100);
    max_slots = latency / 2 / (2 * MINttx + trx);
    
    *wd->period = DBL_MAX;
    *wd->slots = 0;

    pthread_mutex_lock(wd->task_mutex);
    pthread_sem_down(1, wd->sem_worker_available, wd->task_mutex);
    
    for (i = 0; i < max_slots; i++) {
        task->slot = i;
        task->T = latency / (i + 1);
        lambda = get_lambda(task->T);
        task->lb = 2 * lambda;

        if (2 * M_PI * MINttx / task->T > task->lb)
            task->lb = 2 * M_PI * MINttx / task->T;

        max_samples = (M_PI - lambda) / task->lb - 1;

        if (energy_per_time(task->T, task->lb, 1) > max_energy) {
            printf("stopping at %d periods, as min(I)=%.2f mA from now\n",
                    i + 1, energy_per_time(task->T, task->lb, 1) / 100);
            break;
        }

        if (*wd->slots > 0) 
            break;

        for (j = 1; j <= max_samples; j++) {
            task->ub = (M_PI - lambda) / (j + 1);
            task->pc.lambda = lambda;
            task->pc.samples = j;

            if (energy_per_time(task->T, task->lb, j) > max_energy) {
                printf("stopping samples at %d, as min(I)=%.2f from now\n", j, 
                        energy_per_time(task->T, task->lb, j));
                break;
            }
            if (*wd->slots > 0) 
                break;

            pthread_sem_up(1, wd->sem_new_task);
            pthread_sem_down(1, wd->sem_task_buffered, wd->task_mutex);
            pthread_sem_down(1, wd->sem_worker_available, wd->task_mutex);
        }
    }
    pthread_sem_down(thread_num - 1, wd->sem_worker_available, wd->task_mutex);
    pthread_sem_up(thread_num, wd->sem_worker_available);
    pthread_mutex_unlock(wd->task_mutex);

    return *wd->slots * *wd->period;
}


double get_lifetime_params(double lifetime, double probability, double *period,
        protocol_params_t *params)
{
    double lb, ub, middle;
    double last_latency, actual_latency;
    unsigned long calls;
    double max_energy = BATTERY / lifetime;
    int thread_num = sysconf(_SC_NPROCESSORS_ONLN);

    pthread_t *threads;
    pthread_sem_t sem_worker_available;
    pthread_sem_t sem_new_task;
    pthread_sem_t sem_task_buffered;
    pthread_mutex_t task_mutex = PTHREAD_MUTEX_INITIALIZER;
    int finish = 0;
    int slots = 0;
    int i;

    struct worker_task task;
    struct worker_data worker_data = {
        .probability = probability,
        .task = &task,
        
        .finish = &finish,
        
        .energy = &max_energy,
        .slots = &slots,
        .params = params,
        .period = period,

        .sem_new_task = &sem_new_task,
        .sem_task_buffered = &sem_task_buffered,
        .sem_worker_available = &sem_worker_available,
        .task_mutex = &task_mutex,
    };

    assert(period != NULL);
    assert(params != NULL);

    pthread_sem_init(0, &sem_worker_available);
    pthread_sem_init(0, &sem_new_task);
    pthread_sem_init(0, &sem_task_buffered);

    threads = malloc(thread_num * sizeof(pthread_t));
    for (i = 0; i < thread_num; i++)
        pthread_create(threads + i, NULL, worker_thread, &worker_data);

    lb = 4 * MINttx;
    ub = MAXLATENCY;
    middle = (ub - lb) / 2 + lb;
 
    last_latency = ub;
    actual_latency = try_latency(thread_num, last_latency, max_energy, 
            probability, &worker_data);

    if (actual_latency == 0) {
        actual_latency = DBL_MAX;
        goto lifetime_terminate;
    }

    for (calls = 0; calls < MAX_CALLS * 100; calls++) {
        actual_latency = try_latency(thread_num, middle, max_energy, 
                probability, &worker_data);

        if (actual_latency != 0) {
            double delta;

            delta = fabs(middle - last_latency);
            if (delta / last_latency < TOL_REL) {
                break;
            }
            last_latency = ub = middle;
        } else 
            lb = middle;

        middle = (ub - lb) / 2 + lb;
    }

lifetime_terminate:
    finish = 1;
    pthread_mutex_lock(&task_mutex);
    pthread_sem_up(thread_num, &sem_new_task);
    pthread_mutex_unlock(&task_mutex);

    printf("waiting for all workers\n");
    for (i = 0; i < thread_num; i++)
        pthread_join(threads[i], NULL);
    
    free(threads);

    return actual_latency;
}

