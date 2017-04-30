#ifndef FLUID_SIM
#define FLUID_SIM

#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

struct FluidSim {
    // Size of time step
    double dt = 0.5;
    // Default masses for each particle
    double mass = 1;
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

    // Should just call particles.reset()
    void reset();
    // Calls step particles.step() with appropriate params
    void step();
};

#endif FLUID_SIM