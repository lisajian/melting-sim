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
        dt = 0.5;
        mass = 1;
        default_forces = {glm::dvec3(0, -9.8, 0)};
        rho = 1;
        h = 3;
        k = 0.1;
        del_q = 0.1 * h;
        n = 4;
        eps = 0.05;
        all_particles = Particles();
    }

    // Size of time step
    double dt;
    // Default masses for each particle
    double mass;
    // Default external forces
    std::vector<glm::dvec3> default_forces;
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
    Particles all_particles;

    void reset();
    void step();
};

#endif FLUID_SIM