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
            // TODO: Have a separate method for verlett integration
            double d = 0.05;
            glm::dvec3 next_step = p.curr_pos + (1 - d) * p.vdt + (p.forces / (double) p.mass);
            bool c = false;

            // Collided with bounding box; make adjustment in position
            // TODO: Dampening doesn't look right
            // TODO: Implement self collisions
            // TODO: Bounce's angle seems a bit strange..
            if (next_step.x > b_max.x || next_step.x < b_min.x) {
                p.vdt.x *= -1;
                // p.last_pos = p.curr_pos;
                // p.curr_pos.x *= -1;
                c = true;
            }
            if (next_step.y > b_max.y || next_step.y < b_min.y) {
                p.vdt.y *= -1;
                // p.last_pos = p.curr_pos;
                // p.curr_pos.y *= -1;
                c = true;
            }
            if (next_step.z > b_max.z || next_step.z < b_min.z) {
                p.vdt.z *= -1;
                // p.last_pos = p.curr_pos;
                // p.curr_pos.z *= -1;
                c = true;
            }
            return c;
        }
};

#endif /* BBOX_H */











