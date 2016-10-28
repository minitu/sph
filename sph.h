#include <iostream>
#include <cmath>

#define VERSION_TAG "SPHView01"

typedef struct sim_param_t_ {
  char *fname;                // file name
  int nframes;                // number of frames
  int npframe;                // steps per frame
  float h;                    // particle size
  float dt;                   // time step
  float rho0;                 // reference density
  float k;                    // bulk modulus
  float mu;                   // viscosity
  float g;                    // gravity strength
} sim_param_t;

int get_params(int argc, char **argv, sim_param_t *params);

typedef struct sim_state_t_ {
  // FIXME move n and mass to param?
  int n;                      // number of particles
  float mass;                 // particle mass
  float *__restrict__ rho;    // densities
  float *__restrict__ x;      // positions
  float *__restrict__ vh;     // velocities (half step)
  float *__restrict__ v;      // velocities (full step)
  float *__restrict__ a;      // acceleration
} sim_state_t;

sim_state_t *alloc_state(int n); // TODO
void free_state(sim_state_t *state); // TODO

uint32_t htonf(void *data) {
  return htonl(*(uint32_t*)data);
}
