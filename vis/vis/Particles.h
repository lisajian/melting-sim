/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Particles.h
 * Author: swl
 *
 * Created on April 15, 2016, 12:16 PM
 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <iostream>
#include "bbox.h"
#include "particle.h"
// #include <glm/glm.hpp>
#include <vector>
// #if defined(__APPLE_CC__)
// #include <GLUT/glut.h>
// #else
// #include <GL/glut.h>
// #include <math.h>
// #endif


class Particles {
    public:
        Particles();
        // Particles(BBox b) : bbox(b) {} // TODO: Particles don't need to know about bounding box...

        void render() const;
        void reset();
        void step(double dt, double h, double rho, double eps, double k, \
                  double const_n, double del_q, int solverIterations); // simulate one frame
        void find_neighboring(double h, Particle &p);
    private:
        std::vector<Particle> particles;
        BBox bbox;
        double default_mass;
        std::vector<glm::dvec3> default_forces;
};

#endif /* PARTICLES_H */

