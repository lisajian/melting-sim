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
        bool collides(Particle &p, double dt) {
            // TODO: Have a separate method for verlett integration
            double d = 0.15;
            // glm::dvec3 next_step = p.pred_pos + p.del_p;
            glm::dvec3 next_step = p.pred_pos;
            bool c = false;
            double clamp = 0.1;
            p.wall_collide = glm::dvec3(1, 1, 1);
            p.wall_collide_f = glm::dvec3(1, 1, 1);

            // Collided with bounding box; make adjustment in position
            // TODO: Implement self collisions
            if (next_step.x > b_max.x || next_step.x < b_min.x) {
                double adj = fmin(std::abs(next_step.x - b_min.x), std::abs(next_step.x - b_max.x));
                p.pred_pos += adj * (p.curr_pos - next_step);
                p.wall_collide = glm::dvec3(-(1 - d), 1 - d, 1 - d);
                c = true;
                std::cout << "out of x" << std::endl;
                std::cout << p.curr_pos.x << ", " << p.curr_pos.y << ", " << p.curr_pos.z << std::endl;
            }
            if (next_step.y > b_max.y || next_step.y < b_min.y) {
                double adj = fmin(std::abs(next_step.y - b_min.y), std::abs(next_step.y - b_max.y));
                p.pred_pos += adj * (p.curr_pos - next_step);
                p.wall_collide = glm::dvec3(1 - d, -(1 - d), 1 - d);
                c = true;
                std::cout << "out of y" << std::endl;
                std::cout << p.curr_pos.x << ", " << p.curr_pos.y << ", " << p.curr_pos.z << std::endl;
            }
            if (next_step.z > b_max.z || next_step.z < b_min.z) {
                double adj = fmin(std::abs(next_step.z - b_min.z), std::abs(next_step.z - b_max.z));
                p.pred_pos += adj * (p.curr_pos - next_step);
                p.wall_collide = glm::dvec3(1 - d, 1 - d, -(1 - d));
                c = true;
                std::cout << "out of z" << std::endl;
                std::cout << p.curr_pos.x << ", " << p.curr_pos.y << ", " << p.curr_pos.z << std::endl;
            }

            double min_dist = 0.05;
            // Stop the particle if the velocity is small and it's close to
            // a face of the bounding box
            // TODO: Do this for all faces
            if (next_step.y - b_min.y < 0) {
                std::cout << "===================" << std::endl; 
                std::cout << "p.id: " << p.id << std::endl;
                std::cout << "glm::length((1.0 / dt) * (p.pred_pos + p.del_p - p.curr_pos)): " << glm::length((1.0 / dt) * (p.pred_pos + p.del_p - p.curr_pos)) << std::endl;
                std::cout << "glm::length(p.vdt): " << glm::length(p.vdt) << std::endl;
                std::cout << "std::abs(next_step.y - b_min.y): " << std::abs(next_step.y - b_min.y) << std::endl;
                // std::cout << "std::abs(next_step.y - b_max.y): " << std::abs(next_step.y - b_max.y) << std::endl;
                // std::cout << "std::abs(next_step.x - b_min.x): " << std::abs(next_step.x - b_min.x) << std::endl;
                // std::cout << "std::abs(next_step.x - b_max.x): " << std::abs(next_step.x - b_max.x) << std::endl;
                // std::cout << "std::abs(next_step.z - b_min.z): " << std::abs(next_step.z - b_min.z) << std::endl;
                // std::cout << "std::abs(next_step.z - b_max.z): " << std::abs(next_step.z - b_max.z) << std::endl;
            }
            // std::cout << "===================" << std::endl; 
            // std::cout << "p.id: " << p.id << std::endl;
            // std::cout << "glm::length((1.0 / dt) * (p.pred_pos + p.del_p - p.curr_pos)): " << glm::length((1.0 / dt) * (p.pred_pos + p.del_p - p.curr_pos)) << std::endl;
            // std::cout << "glm::length(p.vdt): " << glm::length(p.vdt) << std::endl;
            // std::cout << "std::abs(next_step.y - b_min.y): " << std::abs(next_step.y - b_min.y) << std::endl;
            // std::cout << "std::abs(next_step.y - b_max.y): " << std::abs(next_step.y - b_max.y) << std::endl;
            // std::cout << "std::abs(next_step.x - b_min.x): " << std::abs(next_step.x - b_min.x) << std::endl;
            // std::cout << "std::abs(next_step.x - b_max.x): " << std::abs(next_step.x - b_max.x) << std::endl;
            // std::cout << "std::abs(next_step.z - b_min.z): " << std::abs(next_step.z - b_min.z) << std::endl;
            // std::cout << "std::abs(next_step.z - b_max.z): " << std::abs(next_step.z - b_max.z) << std::endl;
            if (glm::length((1.0 / dt) * (p.pred_pos + p.del_p - p.curr_pos)) < clamp && next_step.y - b_min.y < 0) {
                p.wall_collide = glm::dvec3(1-d, 0, 1-d);
                // p.wall_collide_f = glm::dvec3(1, 0, 1);
                p.pred_pos.y = b_min.y;
            }
            return c;
        }
};

#endif /* BBOX_H */











