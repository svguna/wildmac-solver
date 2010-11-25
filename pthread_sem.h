/*
 * wildmac-solver - returns the proper configuration of the wildmac protocol,
 * given a desired detection latency and probability.
 * Copyright (C) 2010  Stefan Guna, http://disi.unitn.it/~guna
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
#ifndef __PTHREAD_SEM_H
#define __PTHREAD_SEM_H

#include <pthread.h>

struct pthread_sem {
    int val;
    pthread_cond_t cond;
};
typedef struct pthread_sem pthread_sem_t;


int pthread_sem_down(int val, pthread_sem_t *sem, pthread_mutex_t *mutex);
int pthread_sem_up(int val, pthread_sem_t *sem);
int pthread_sem_init(int val, pthread_sem_t *sem);

#endif
