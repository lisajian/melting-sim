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
        : last_pos(init_pos), curr_pos(init_pos), vdt(glm::dvec3(0, 0, 0)) {}

    // Particle(glm::dvec3 init_pos, glm::dvec3 init_velo)
    //     : last_pos(init_pos), curr_pos(init_pos), vdt(init_velo) {}

    glm::dvec3 last_pos;
    glm::dvec3 curr_pos;
    glm::dvec3 vdt;
    glm::dvec3 forces;
    float mass = 800; // TODO: Masses shouldn't be this heavy...
};

#endif /* PARTICLE_H */