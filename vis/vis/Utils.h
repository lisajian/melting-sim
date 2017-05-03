#ifndef UTILS
#define UTILS

#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

// poly6_const := (315.0 / (64.0 * M_PI * pow(h, 9.0)))
double poly6_kernel(glm::dvec3 r, double poly6_const, double h) {
    double r2 = glm::dot(r, r);
    double h2 = h * h;
    if (r2 <= h2 && h >= 0) {
        return poly6_const * std::pow(h * h - glm::dot(r, r), 3.0);
    }
    return 0;
}

// spiky_const := (-45.0 / (M_PI * pow(h, 6.0)))
glm::dvec3 spiky_kernel_grad(glm::dvec3 r, double spiky_const, double h) {
    double rsq = glm::length(r);
    if (rsq <= h) {
        return std::pow(h - rsq, 2.0) * spiky_const * r / rsq;
    }
    return glm::dvec3(0, 0, 0);
}

#endif UTILS