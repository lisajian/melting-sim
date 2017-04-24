// Single particle struct

#ifndef PARTICLE_H
#define PARTICLE_H

#include <glm/glm.hpp>
#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

struct Particle {
    Particle(glm::dvec3 init_pos)
        : p(init_pos) {}

    Particle(glm::dvec3 init_pos, glm::dvec3 init_velo)
        : p(init_pos), v(init_velo) {}

    glm::dvec3 p;
    glm::dvec3 v;
    glm::dvec3 forces;
    float mass = 800; // TODO: Masses shouldn't be this heavy...
};

#endif /* PARTICLE_H */