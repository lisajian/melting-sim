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
        rho = 16;
        h = 0.15;
        k = 0.001;
        del_q = 0.01 * h;
        n = 4;
        eps = 1000;
        all_particles = Particles();
        solverIterations = 1;
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
    // Used in large simulation loop
    int solverIterations;
    Particles all_particles;

    void reset();
    void step();
    void render();
};

#endif FLUID_SIM