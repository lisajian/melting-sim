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
#include <unordered_map>
// #if defined(__APPLE_CC__)
// #include <GLUT/glut.h>
// #else
// #include <GL/glut.h>
// #include <math.h>
// #endif


class Particles {
    public:
        Particles(double const_h);
        // Particles(BBox b) : bbox(b) {} // TODO: Particles don't need to know about bounding box...

        void render() const;
        void reset();
        void step(double dt, double h, double rho, double eps, double k, \
                  double const_n, double del_q, int solverIterations, double c, double eps_vort); // simulate one frame
        void find_neighboring(double h, Particle &p);
        void build_spatial_map();
        float hash_position(glm::dvec3 pos);
        void particle_collisions(Particle &p);
        void calc_spiky(Particle &p, double spiky_const, double h);
        void calc_poly6(Particle &p, double poly6_const, double h);
    private:
        std::vector<Particle> particles;
        BBox bbox;
        std::unordered_map<float, std::vector<Particle *> *> map;
        double default_mass;
        std::vector<glm::dvec3> default_forces;
        int nx;
        int ny;
        int nz;
        double w;
        double h;
        double t;
        double poly6_const;
        double spiky_const;
};

#endif /* PARTICLES_H */
