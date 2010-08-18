#ifndef __WILDMAC_H
#define __WILDMAC_H

struct protocol_params {
    double tau;
    double lambda;
    int samples;
    double on; // solely to improve performance
};
typedef struct protocol_params protocol_params_t;

#define ACTIVE(p) (2 * M_PI - (p)->on)
#define SET_ON(p) (p)->on = ((p)->samples + 1) * (p)->tau + (p)->lambda

#define HASH_SIZE 16

#endif

