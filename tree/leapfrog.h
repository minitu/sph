#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "state.h"
#include "params.h"

void leapfrog_start(sim_state_t* s, sim_param_t* param, double dt);
void leapfrog_step(sim_state_t* s, sim_param_t* param, double dt);

#endif /* LEAPFROG_H */
