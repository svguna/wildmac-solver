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
