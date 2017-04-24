// An axis aligned bounding box

#ifndef BBOX_H
#define BBOX_H

// #include <iostream>
#include "particle.h"
// #include <glm/glm.hpp>
// #include <vector>
// #if defined(__APPLE_CC__)
// #include <GLUT/glut.h>
// #else
// #include <GL/glut.h>
// #include <math.h>
// #endif

struct BBox {
    // Constructor for an invalid bounding box
    public:
        BBox() {
            double inf = std::numeric_limits<double>::infinity();
            b_min = glm::dvec3(inf, inf, inf);
            b_max = glm::dvec3(-inf, -inf, -inf);
        }

        // Constructor for a bounding box with specified min and max points
        BBox(glm::dvec3 init_min, glm::dvec3 init_max)
            : b_min(init_min), b_max(init_max) {}

        glm::dvec3 b_min;
        glm::dvec3 b_max;

        // Returns true if p will collide with box in given
        // time step and adjusts velocity vector accordingly
        bool collides(Particle &p) {
            glm::dvec3 next_step = p.v + p.p;
            bool c = false;

            // Collided with bounding box
            // TODO: Dampening doesn't look right
            // TODO: Implement self collisions
            // TODO: Bounce's angle seems a bit strange..
            if (next_step[0] > b_max[0] || next_step[0] < b_min[0]) {
                p.v[0] *= -1;
                // p.p[0] += 0.01;
                c = true;
            }
            if (next_step[1] > b_max[1] || next_step[1] < b_min[1]) {
                p.v[1] *= -1;
                // p.p[1] += 0.01;
                c = true;
            }
            if (next_step[2] > b_max[2] || next_step[2] < b_min[2]) {
                p.v[2] *= -1;
                // p.p[2] += 0.01;
                c = true;
            }
            return c;
        }
};

#endif /* BBOX_H */











