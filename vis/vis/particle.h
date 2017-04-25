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
    // Particle(glm::dvec3 init_pos)
    //     : last_pos(init_pos), curr_pos(init_pos) {}

    Particle(glm::dvec3 init_pos)
        : curr_pos(init_pos) {
        	last_pos = (init_pos); // TODO: Need a better way of getting initial vel
        }

    // Particle(glm::dvec3 init_pos, glm::dvec3 init_velo)
    //     : last_pos(init_pos), curr_pos(init_pos), vdt(init_velo) {}

    glm::dvec3 last_pos;
    glm::dvec3 curr_pos;
    glm::dvec3 vdt;
    glm::dvec3 new_vdt;
    glm::dvec3 forces;
    glm::dvec3 adjustment_vec;
    std::vector<Particle> neighbors;
    float mass = 1000; // TODO: Masses shouldn't be this heavy...

    // Collides with other particle, assumes other particle has same mass
    // Stores new velocity in self.new_vdt
    // Source: https://en.wikipedia.org/wiki/Elastic_collision
    bool particle_collide(Particle other) {
    	double min_dist = 0.1; // TODO: Have a better way of keeping track of min allotted dist between particles
    	if (glm::length(curr_pos - other.curr_pos) < min_dist) {
    		new_vdt = other.vdt;
    		return true;
    	}
    	return false;
    }

    // Apply correction vector if particles are within min_dist of each other
    bool particle_adjust(Particle other) {
    	double min_dist = 0.1; // TODO: Have a better way of keeping track of min allotted dist between particles
    	double dist = glm::length(curr_pos - other.curr_pos);
    	if (dist < min_dist) {
    		glm::dvec3 correction = last_pos - curr_pos;
    		// TODO: maybe just adjust by vdt?
    		adjustment_vec += correction * ((min_dist - dist) / glm::length(correction)); 
    		return true;
    	}
    	return false;
    }
};

#endif /* PARTICLE_H */