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
            neighbors = std::vector<Particle>();
    }

    glm::dvec3 curr_pos;
    glm::dvec3 pred_pos;
    glm::dvec3 del_p;
    glm::dvec3 vdt;
    glm::dvec3 new_vdt;
    glm::dvec3 forces;
    glm::dvec3 adjustment_vec;
    glm::dvec3 wall_collide;
    glm::dvec3 wall_collide_f;
    glm::dvec3 particle_collide;
    std::vector<Particle> neighbors;
    float mass;
    double C;
    double lambda;
    int id; // TODO: remove this

};

#endif /* PARTICLE_H */