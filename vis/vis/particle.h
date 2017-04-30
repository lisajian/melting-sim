// Single particle struct

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <glm/glm.hpp>
#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

struct Particle {
    // Particle(glm::dvec3 init_pos)
    //     : last_pos(init_pos), curr_pos(init_pos) {}

    Particle(glm::dvec3 init_pos)
        : curr_pos(init_pos) {
            last_pos = (init_pos); // TODO: Need a better way of getting initial vel
            neighbors = std::vector<Particle>();
    }

    glm::dvec3 last_pos;
    glm::dvec3 curr_pos;
    glm::dvec3 del_p;
    glm::dvec3 vdt;
    glm::dvec3 new_vdt;
    glm::dvec3 forces;
    glm::dvec3 adjustment_vec;
    std::vector<Particle> neighbors;
    float mass = 100000; // TODO: Masses shouldn't be this heavy...
    double C = 0.0;
    double lambda = 0.0;

    // TODO: need a predicted position, pos*

};

#endif /* PARTICLE_H */