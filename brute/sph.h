#include <cstdint>
#include <string>
#include <cassert>

#define VERSION_TAG "SPHView01"

using namespace std;

/* Parameters */
typedef struct sim_param_t_ {
  string fname;                // file name
  int nframes;                // number of frames
  int npframe;                // steps per frame
  float h;                    // particle size
  float dt;                   // time step
  float rho0;                 // reference density
  float k;                    // bulk modulus
  float mu;                   // viscosity
  float g;                    // gravity strength
} sim_param_t;

static void default_params(sim_param_t *params);
int get_params(int argc, char **argv, sim_param_t *params);
static void print_usage();

/* State */
typedef struct sim_state_t_ {
  int n;                      // number of particles
  float mass;                 // particle mass
  float *__restrict__ rho;    // densities
  float *__restrict__ x;      // positions
  float *__restrict__ vh;     // velocities (half step)
  float *__restrict__ v;      // velocities (full step)
  float *__restrict__ a;      // acceleration
} sim_state_t;

sim_state_t *alloc_state(int n);
void free_state(sim_state_t *state);

/* Computations */
void compute_density(sim_state_t *state, sim_param_t *params);
void compute_accel(sim_state_t *state, sim_param_t *params);

/* Leapfrog integration */
void leapfrog_step(sim_state_t *state, double dt);
void leapfrog_start(sim_state_t *state, double dt);

/* Reflection boundary conditions */
static void damp_reflect(int which, float barrier, float *x, float *v, float *vh);
static void reflect_bc(sim_state_t *state);

/* Initialization */
typedef int (*domain_fun_t)(float, float);
int box_indicator(float x, float y);
int circ_indicator(float x, float y);
sim_state_t *place_particles(sim_param_t *param, domain_fun_t indicatef);
void normalize_mass(sim_state_t *state, sim_param_t *param);
sim_state_t *init_particles(sim_param_t *param);

/* Output */
void write_header(FILE *fp, int n);
void write_frame_data(FILE *fp, int n, float *x, int *c);
uint32_t htonf(void *data);
