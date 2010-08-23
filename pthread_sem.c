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
    sem->val += val;
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

