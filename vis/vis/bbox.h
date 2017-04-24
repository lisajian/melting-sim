// An axis aligned bounding box

#ifndef BBOX_H
#define BBOX_H

#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

struct BBox {
    // Constructor for an invalid bounding box
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

    bool collides(Particle p);
};

#endif /* BBOX_H */
