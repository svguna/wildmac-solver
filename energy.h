#ifndef __ENERGY_H
#define __ENERGY_H

#include "wildmac.h"

// Currents are given in tens of uA.
// Time is given in tens of us.

#define Itx 1740
#define Irx 1970
#define Ioff 2
#define trx 25
#define MINttx 1000
#define MAX_EVAL 100

#define MAX_CALLS 1000
#define TOL_REL 1e-2

double get_protocol_parameters(double latency, double probability, 
        double *period, protocol_params_t *params);



#endif

