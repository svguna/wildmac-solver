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
#include <pthread.h>
#include <assert.h>
#include "pthread_sem.h"

int pthread_sem_down(int val, pthread_sem_t *sem, pthread_mutex_t *mutex)
{
    if (sem->val - val >= 0) {
        sem->val -= val;
        return 0;
    }

    do {
        pthread_cond_wait(&sem->cond, mutex);
    } while (sem->val - val < 0);
    
    sem->val -= val;
    return 0;
}


int pthread_sem_up(int val, pthread_sem_t *sem)
{
    int i;
    sem->val += val;
    for (i = 0; i < val; i++)
        pthread_cond_signal(&sem->cond);
    return 0;
}


int pthread_sem_init(int val, pthread_sem_t *sem)
{
    assert(sem != NULL);

    sem->val = val;
    pthread_cond_init(&sem->cond, NULL);
    return 0;
}


int pthread_sem_destroy(pthread_sem_t *sem)
{
    pthread_cond_destroy(&sem->cond);
    return 0;
}

