#include <cstdio>
#include "sph.h"

void check_state(sim_state_t *state) {
  for (int i = 0; i < state->n; i++) {
    float xi = state->x[2*i+0];
    float yi = state->x[2*i+1];
    assert(xi >= 0 || xi <= 1);
    assert(yi >= 0 || yi <= 1);
  }
}

int main(int argc, char **argv) {
  sim_param_t params;
  if (get_params(argc, argv, &params) != 0)
    exit(-1);
  
  sim_state_t *state = init_particles(&params);
  FILE *fp = fopen(params.fname.c_str(), "w");
  int nframes = params.nframes;
  int npframe = params.npframe;
  float dt = params.dt;
  int n = state->n;
  clock_t t = clock();

  write_header(fp, n);
  write_frame_data(fp, n, state->x, NULL);
  compute_accel(state, &params);
  leapfrog_start(state, dt);
  check_state(state);
  for (int frame = 1; frame < nframes; frame++) {
    for (int i = 0; i < npframe; i++) {
      compute_accel(state, &params);
      leapfrog_step(state, dt);
      check_state(state);
    }
    write_frame_data(fp, n, state->x, NULL);
  }
  printf("Ran in %f seconds\n", ((float)(clock() - t) / CLOCKS_PER_SEC));

  fclose(fp);
  free_state(state);
}
