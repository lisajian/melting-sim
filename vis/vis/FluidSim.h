#ifndef FLUID_SIM
#define FLUID_SIM

#include "Particles.h"

#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

struct FluidSim {
    // Default settings
    FluidSim() {
        dt = 0.1;
        rho = 100;
        h = 2;
        k = 0.1;
        del_q = 0.1 * h;
        c = 0.01;
        n = 4;
        eps = 1000;
        eps_vort = 0.01;
        all_particles = Particles(h);
        solverIterations = 1;

        // dt = 0.1;
        // rho = 100;
        // h = 2;
        // k = 0.1;
        // del_q = 0.1 * h;
        // c = 0.01;
        // n = 4;
        // eps = 1000;
        // eps_vort = 0.01;
        // all_particles = Particles(h);
        // solverIterations = 1;
    }

    // Size of time step
    double dt;
    // Default rest density
    double rho;
    // Used in all smoothing calculations, h
    // := maximum radius between "neighboring" particles
    double h;
    // Used in tensile instability, k;
    // usually k = 0.1 works well
    double k;
    // Used in tensile instability, del_q;
    // usually del_q = 0.1h to 0.3h works well
    double del_q;
    // Used in tensile instability, n
    double n;
    // Used in determining lambda_i's, relaxation term
    double eps;
    // Used in vorticity calculation
    double eps_vort;
    // Used in large simulation loop
    int solverIterations;
    // Used in XSPH viscocity calculation
    double c;
    Particles all_particles = Particles(1);

    void reset();
    void step();
    void render();
};

#endif FLUID_SIM
