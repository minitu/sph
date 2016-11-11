#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "io.h"
#include "params.h"
#include "state.h"
#include "interact.h"
#include "leapfrog.h"
#include "timing.h"

extern float density_time;
extern float accel_time;

/*@q
 * ====================================================================
 */

/*@T
 * \section{Initialization}
 *
 * We've hard coded the computational domain to a unit box, but we'd prefer
 * to do something more flexible for the initial distribution of fluid.
 * In particular, we define the initial geometry of the fluid in terms of an
 * {\em indicator function} that is one for points in the domain occupied
 * by fluid and zero elsewhere.  A [[domain_fun_t]] is a pointer to an
 * indicator for a domain, which is a function that takes two floats and
 * returns 0 or 1.  Two examples of indicator functions are a little box
 * of fluid in the corner of the domain and a circular drop.
 *@c*/
typedef int (*domain_fun_t)(float, float, float);

int box_indicator(float ww, float x, float y)
{
    return (x < ww/2) && (y < ww/2);
}

int circ_indicator(float ww, float x, float y)
{
    float dx = (x-ww/2);
    float dy = (y-ww/3);
    float r2 = dx*dx + dy*dy;
    return (r2 < ww/4*ww/4);
}

/*@T
 *
 * The [[place_particles]] routine fills a region (indicated by the
 * [[indicatef]] argument) with fluid particles.  The fluid particles
 * are placed at points inside the domain that lie on a regular mesh
 * with cell sizes of $h/1.3$.  This is close enough to allow the
 * particles to overlap somewhat, but not too much.
 *@c*/
sim_state_t* place_particles(sim_param_t* param, 
                             domain_fun_t indicatef)
{
    float h  = param->h;
    float ww = param->ww;
    float hh = h/1.3;

    // Count mesh points that fall in indicated region.
    int count = 0;
    for (float x = 0; x < ww; x += hh)
        for (float y = 0; y < ww; y += hh)
            count += indicatef(ww, x,y);

    // Populate the particle data structure
    sim_state_t* s = alloc_state(count);
    int p = 0;
    for (float x = 0; x < ww; x += hh) {
        for (float y = 0; y < ww; y += hh) {
            if (indicatef(ww,x,y)) {
                s->x[2*p+0] = x;
                s->x[2*p+1] = y;
                s->v[2*p+0] = 0;
                s->v[2*p+1] = 0;
                ++p;
            }
        }
    }
    return s;    
}

/*@T
 *
 * The [[place_particle]] routine determines the initial particle
 * placement, but not the desired mass.  We want the fluid in the
 * initial configuration to exist roughly at the reference density.
 * One way to do this is to take the volume in the indicated body of
 * fluid, multiply by the mass density, and divide by the number of
 * particles; but that requires that we be able to compute the volume
 * of the fluid region.  Alternately, we can simply compute the
 * average mass density assuming each particle has mass one, then use
 * that to compute the particle mass necessary in order to achieve the
 * desired reference density.  We do this with [[normalize_mass]].
 * @c*/
void normalize_mass(sim_state_t* s, sim_param_t* param)
{
    s->mass = 1;
    compute_density(s, param);
    float rho0 = param->rho0;
    float rho2s = 0;
    float rhos  = 0;
    for (int i = 0; i < s->n; ++i) {
        rho2s += (s->rho[i])*(s->rho[i]);
        rhos  += s->rho[i];
    }
    s->mass *= ( rho0*rhos / rho2s );
}

sim_state_t* init_particles(sim_param_t* param)
{
    sim_state_t* s = place_particles(param, circ_indicator);
    normalize_mass(s, param);
    return s;
}

/*@T
 * \section{The [[main]] event}
 *
 * The [[main]] routine actually runs the time step loop, writing
 * out files for visualization every few steps.  For debugging
 * convenience, we use [[check_state]] before writing out frames,
 * just so that we don't spend a lot of time on a simulation that
 * has gone berserk.
 *@c*/

void check_state(sim_state_t* s, sim_param_t* param)
{
    for (int i = 0; i < s->n; ++i) {
        float xi = s->x[2*i+0];
        float yi = s->x[2*i+1];
        assert( xi >= 0 || xi <= param->ww );
        assert( yi >= 0 || yi <= param->ww );
    }
}

int main(int argc, char** argv)
{
    sim_param_t params;
    if (get_params(argc, argv, &params) != 0)
        exit(-1);
    sim_state_t* state = init_particles(&params);
    FILE* fp    = fopen(params.fname, "w");
    int nframes = params.nframes;
    int npframe = params.npframe;
    float dt    = params.dt;
    int n       = state->n;

    tic(0);
    write_header(fp, n, params.ww);
    write_frame_data(fp, n, state->x, NULL);
    compute_accel(state, &params);
    leapfrog_start(state, &params, dt);
    check_state(state, &params);
    for (int frame = 1; frame < nframes; ++frame) {
        for (int i = 0; i < npframe; ++i) {
            compute_accel(state, &params);
            leapfrog_step(state, &params, dt);
            check_state(state, &params);
        }
        write_frame_data(fp, n, state->x, NULL);
    }
    printf("Ran in %g seconds\n\tdensity_time %g\n\taccel_time %g\n", toc(0), density_time, accel_time);

    fclose(fp);
    free_state(state);
}
