#include <iostream>
#include <getopt.h>
#include <cmath>
#include <cstring>
#include <arpa/inet.h>
#include "sph.h"

static void default_params(sim_param_t *params) {
  params->fname = "run.out";
  params->nframes = 400;
  params->npframe = 100;
  params->dt = 1e-4;
  params->h = 5e-2;
  params->rho0 = 1000;
  params->k = 1e3;
  params->mu = 0.1;
  params->g = 9.8;
}

int get_params(int argc, char **argv, sim_param_t *params) {
  extern char *optarg;
  const char *optstring = "h:o:F:f:t:s:d:k:v:g:";
  int c;

  default_params(params);
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
      case 'h':
        print_usage();
        return -1;
      case 'o':
        params->fname = optarg; break;
      case 'F':
        params->nframes = atoi(optarg); break;
      case 'f':
        params->npframe = atoi(optarg); break;
      case 't':
        params->dt = (float)atof(optarg); break;
      case 's':
        params->h = (float)atof(optarg); break;
      case 'd':
        params->rho0 = (float)atof(optarg); break;
      case 'k':
        params->k = (float)atof(optarg); break;
      case 'v':
        params->mu = (float)atof(optarg); break;
      case 'g':
        params->g = (float)atof(optarg); break;
      default:
        fprintf(stderr, "Unknown option\n");
        return -1;
    }
  }
  return 0;
}

static void print_usage() {
  sim_param_t param;
  default_params(&param);
  fprintf(stderr,
      "SPH\n"
      "\t-h: print this message\n"
      "\t-o: output file name (%s)\n"
      "\t-F: number of frames (%d)\n"
      "\t-f: steps per frame (%d)\n"
      "\t-t: time step (%e)\n"
      "\t-s: particle size (%e)\n"
      "\t-d: reference density (%g)\n"
      "\t-k: bulk modulus (%g)\n"
      "\t-v: dynamic viscosity (%g)\n"
      "\t-g: gravitational strength (%g)\n",
      param.fname.c_str(), param.nframes, param.npframe,
      param.dt, param.h, param.rho0,
      param.k, param.mu, param.g);
}

sim_state_t *alloc_state(int n) {
  sim_state_t *state = (sim_state_t *)malloc(sizeof(sim_state_t));
  state->n = n;
  state->rho = (float *)malloc(sizeof(float) * n);
  state->x = (float *)malloc(sizeof(float) * n * 2);
  state->vh = (float *)malloc(sizeof(float) * n * 2);
  state->v = (float *)malloc(sizeof(float) * n * 2);
  state->a = (float *)malloc(sizeof(float) * n * 2);

  return state;
}

void free_state(sim_state_t *state) {
  assert(state != NULL);
  free(state);
}

void compute_density(sim_state_t *state, sim_param_t *params) {
  int n = state->n;
  float *__restrict__ rho = state->rho;
  const float *__restrict__ x = state->x;
  
  float h = params->h;
  float h2 = h * h;
  float h8 = (h2 * h2) * (h2 * h2);
  float C = 4 * state->mass / M_PI / h8;

  memset(rho, 0, n*sizeof(float));
  // FIXME brute force
  for (int i = 0; i < n; i++) {
    rho[i] += 4 * state->mass / M_PI / h2; // FIXME move out of the loop
    for (int j = i+1; j < n; j++) {
      float dx = x[2*i+0] - x[2*j+0];
      float dy = x[2*i+1] - x[2*j+1];
      float r2 = dx * dx + dy * dy;
      float z = h2 - r2;
      if (z > 0) {
        float rho_ij = C * z * z * z;
        rho[i] += rho_ij;
        rho[j] += rho_ij;
      }
    }
  }
}

void compute_accel(sim_state_t *state, sim_param_t *params) {
  // unpack basic parameters
  const float h = params->h;
  const float rho0 = params->rho0;
  const float k = params->k;
  const float mu = params->mu;
  const float g = params->g;
  const float mass = state->mass;
  const float h2 = h * h;

  // unpack system state
  const float *__restrict__ rho = state->rho;
  const float *__restrict__ x = state->x;
  const float *__restrict__ v = state->v;
  float *__restrict__ a = state->a;
  int n = state->n;

  // compute density and color
  compute_density(state, params);

  // start with gravity and surface forces
  for (int i = 0; i < n; i++) {
    a[2*i+0] = 0;
    a[2*i+1] = -g;
  }

  // constants for interaction term
  float C0 = mass / M_PI / (h2 * h2);
  float Cp = 15 * k;
  float Cv = -40 * mu;

  // compute interation forces
  // FIXME brute force
  for (int i = 0; i < n; i++) {
    const float rhoi = rho[i];
    for (int j = i + 1; j < n; j++) {
      float dx = x[2*i+0] - x[2*j+0];
      float dy = x[2*i+1] - x[2*j+1];
      float r2 = dx * dx + dy * dy;
      if (r2 < h2) {
        const float rhoj = rho[j];
        float q = sqrt(r2) / h;
        float u = 1 - q;
        float w0 = C0 * u / rhoi / rhoj;
        float wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
        float wv = w0 * Cv;
        float dvx = v[2*i+0] - v[2*j+0];
        float dvy = v[2*i+1] - v[2*j+1];
        a[2*i+0] += (wp * dx + wv * dvx);
        a[2*i+1] += (wp * dy + wv * dvy);
        a[2*j+0] -= (wp * dx + wv * dvx);
        a[2*j+1] -= (wp * dy + wv * dvy);
      }
    }
  }
}

void leapfrog_step(sim_state_t *state, double dt) {
  const float *__restrict__ a = state->a;
  float *__restrict__ vh = state->vh;
  float *__restrict__ v = state->v;
  float *__restrict__ x = state->x;
  int n = state->n;

  for (int i = 0; i < 2 * n; i++) vh[i] += a[i] * dt;
  for (int i = 0; i < 2 * n; i++) v[i] = vh[i] + a[i] * dt / 2;
  for (int i = 0; i < 2 * n; i++) x[i] += vh[i] * dt;

  reflect_bc(state);
}

void leapfrog_start(sim_state_t *state, double dt) {
  const float *__restrict__ a = state->a;
  float *__restrict__ vh = state->vh;
  float *__restrict__ v = state->v;
  float *__restrict__ x = state->x;
  int n = state->n;
  
  for (int i = 0; i < 2 * n; i++) vh[i] = v[i] + a[i] * dt / 2;
  for (int i = 0; i < 2 * n; i++) v[i] += a[i] * dt;
  for (int i = 0; i < 2 * n; i++) x[i] += vh[i] * dt;

  reflect_bc(state);
}

static void damp_reflect(int which, float barrier,
    float *x, float *v, float *vh) {
  // coefficent of restitution
  const float DAMP = 0.75;

  // ignore degenerate cases
  if (v[which] == 0)
    return;

  // scale back the distance traveled based on time from collision
  float tbounce = (x[which] - barrier) / v[which];
  x[0] -= v[0] * (1 - DAMP) * tbounce;
  x[1] -= v[1] * (1 - DAMP) * tbounce;

  // reflect the position and velocity
  x[which] = 2 * barrier - x[which];
  v[which] = -v[which];
  vh[which] = -vh[which];

  // damp the velocities
  v[0] *= DAMP; vh[0] *= DAMP;
  v[1] *= DAMP; vh[1] *= DAMP;
}

static void reflect_bc(sim_state_t *state) {
  // boundaries of the computational domain
  const float XMIN = 0.0;
  const float XMAX = 1.0;
  const float YMIN = 0.0;
  const float YMAX = 1.0;

  float *__restrict__ vh = state->vh;
  float *__restrict__ v = state->v;
  float *__restrict__ x = state->x;
  int n = state->n;

  for (int i = 0; i < n; i++, x +=2, v += 2, vh += 2) {
    if (x[0] < XMIN) damp_reflect(0, XMIN, x, v, vh);
    if (x[0] > XMAX) damp_reflect(0, XMAX, x, v, vh);
    if (x[1] < YMIN) damp_reflect(1, YMIN, x, v, vh);
    if (x[1] > YMAX) damp_reflect(1, YMAX, x, v, vh);
  }
}

int box_indicator(float x, float y) {
  return (x < 0.5) && (y < 0.5);
}

int circ_indicator(float x, float y) {
  float dx = x - 0.5;
  float dy = y - 0.3;
  float r2 = dx * dx + dy * dy;
  return (r2 < 0.25 * 0.25);
}

sim_state_t *place_particles(sim_param_t *param,
    domain_fun_t indicatef) {
  float h = param->h;
  float hh = h / 1.3;

  // count mesh points that fall in indicated region
  int count = 0;
  for (float x = 0; x < 1; x += hh)
    for (float y = 0; y < 1; y += hh)
      count += indicatef(x, y);

  // populate the particle data structure
  sim_state_t *state = alloc_state(count);
  int p = 0;
  for (float x = 0; x < 1; x += hh) {
    for (float y = 0; y < 1; y+= hh) {
      if (indicatef(x, y)) {
        state->x[2*p+0] = x;
        state->x[2*p+1] = y;
        state->v[2*p+0] = 0;
        state->v[2*p+1] = 0;
        p++;
      }
    }
  }

  return state;
}

void normalize_mass(sim_state_t *state, sim_param_t *param) {
  state->mass = 1;
  compute_density(state, param);
  float rho0 = param->rho0;
  float rho2s = 0;
  float rhos = 0;
  for (int i = 0; i < state->n; i++) {
    rho2s += (state->rho[i]) * (state->rho[i]);
    rhos += state->rho[i];
  }
  state->mass *= (rho0 * rhos / rho2s);
}

sim_state_t *init_particles(sim_param_t *param) {
  sim_state_t *state = place_particles(param, box_indicator);
  normalize_mass(state, param);
  return state;
}

void write_header(FILE *fp, int n) {
  float scale = 1.0;
  uint32_t nn = htonl((uint32_t)n);
  uint32_t nscale = htonf(&scale);
  fwrite(&nn, sizeof(nn), 1, fp);
  fwrite(&nscale, sizeof(nscale), 1, fp);
}

void write_frame_data(FILE *fp, int n, float *x, int *c) {
  for (int i = 0; i < n; i++) {
    uint32_t xi = htonf(x++);
    uint32_t yi = htonf(x++);
    fwrite(&xi, sizeof(xi), 1, fp);
    fwrite(&yi, sizeof(yi), 1, fp);
    uint32_t ci0 = c ? *c++ : 0;
    uint32_t ci = htonl(ci0);
    fwrite(&ci, sizeof(ci), 1, fp);
  }
}

uint32_t htonf(void *data) {
  return htonl(*(uint32_t*)data);
}
