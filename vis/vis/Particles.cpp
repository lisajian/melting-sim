/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Particles.cpp
 * Author(s): swl
 * 
 * Created on April 15, 2016, 12:16 PM
 */

#include "Particles.h"
#include <limits>

Particles::Particles() {
    bbox = BBox(glm::dvec3(-2, -2, -2), glm::dvec3(2, 2, 2));
    default_forces = {glm::dvec3(0, -9.8, 0)};
    default_mass = 10;
    nx = 5;
    ny = 1;
    nz = 1;
    reset();
}

// Resets all particles back to default state.
void Particles::reset() {
    // Number of particles in each dimension
    particles.clear();
    float y_offset = 1;
    float z_offset = -2;
    float d = 0.1;
    for(int x=0; x<nx; x++)
    {
        for(int y=0; y<ny; y++)
        {
            for(int z=0; z<nz; z++)
            {
                Particle par(glm::dvec3((x+0.5-nx*0.5)*d + z_offset, (y+0.5)*d-1.0 + y_offset, (z+0.5-nz*0.5)*d));
                par.forces = glm::dvec3(0.0, 0, 0.0); // TODO: Handle external forces better. This is gravity
                for (glm::dvec3 f : default_forces) {
                    par.forces += f;
                }
                par.mass = default_mass;
                particles.push_back(par);
            }
        }
    }
}

// Single update step for all particles
void Particles::step(double dt, double h, double rho, double eps, double k, \
                     double const_n, double del_q, int solverIterations) {
    // Generic time step update
    for (auto &p : particles) {
        p.vdt = p.vdt + dt * p.forces / (double) p.mass;
        p.pred_pos = p.curr_pos + dt * p.vdt;
        // p.curr_pos = p.pred_pos; // TODO: comment out
        // bbox.collides(p, dt); // TODO: comment out
    }

    // return; // TODO: remove this

    // Get neighbors. Distance between point centers is hardcoded
    // double inf = std::numeric_limits<double>::infinity();
    for (auto &p : particles) {
        find_neighboring(h, p);
    }

    for (int i = 0; i < solverIterations; i++) {
        // Determine lambda_i
        for (auto &p : particles) {
            // Use poly6 kernel for rho (pressure)
            double rho_i = 0.0;
            for (auto &n : p.neighbors) {
                glm::dvec3 diff = p.curr_pos - n.curr_pos;
                rho_i += n.mass * std::pow(h * h - glm::dot(diff, diff), 3.0);
            }
            rho_i *= 315.0 / (64.0 * M_PI * std::pow(h, 9.0));
            p.C = rho_i / rho - 1.0;

            // Use spiky kernel for gradients
            glm::dvec3 grad_i_W = glm::dvec3(0, 0, 0);
            double grad_j_W = 0.0;
            for (auto &n : p.neighbors) {
                glm::dvec3 diff = p.curr_pos - n.curr_pos;
                double len = glm::length(diff);
                glm::dvec3 intermed = (pow(h - len, 2.0) / len) * diff;
                intermed *= - 45.0 / (M_PI * std::pow(h, 6.0));
                
                grad_i_W += intermed;
                grad_j_W += glm::dot(intermed, intermed);
            }
            p.lambda = - p.C / (grad_j_W + glm::dot(grad_i_W, grad_i_W) + eps);
        }
        // Determine delta p
        for (auto &p : particles) {
            glm::dvec3 new_del_p = glm::dvec3(0, 0, 0);
            for (auto &n : p.neighbors) {
                glm::dvec3 diff = p.curr_pos - n.curr_pos;

                // Get the gradient
                double len = glm::length(diff);
                glm::dvec3 intermed = (pow(h - len, 2.0) / len) * diff;
                intermed *= - 45.0 / (M_PI * std::pow(h, 6.0));

                // Calculate s_corr using poly6 kernel
                double s_corr = std::pow(h * h - glm::dot(diff, diff), 3.0);
                s_corr /= std::pow(h * h - (del_q * del_q), 3.0);
                s_corr = std::pow(s_corr, const_n) * -1 * k;

                new_del_p += (p.lambda + n.lambda + s_corr) * intermed;
            }
            p.del_p = new_del_p / rho;
            // TODO: Collision detection?
            bbox.collides(p, dt);
        }
        // Update pred_pos
        for (auto &p : particles) {
            p.pred_pos += p.del_p;
            // bbox.collides(p, dt);
        }
    }

    // Update velocity, position, vorticity, XSPH
    for (auto &p : particles) {
        p.vdt = (1.0 / dt) * (p.pred_pos - p.curr_pos);
        p.curr_pos = p.pred_pos;
    }
}

// Find all neighbors n of p such that the distance between the
// centers of n and p is at most h
void Particles::find_neighboring(double h, Particle &p) {
    p.neighbors.clear();
    // // Naive way of finding neighbors
    // int j = 0;
    // for (auto &other : particles) {
    //     if (&p != &other && glm::length(p.curr_pos - other.curr_pos) <= h) {
    //         p.neighbors.push_back(other);
    //     }
    // }
    float hash = hash_position(p.curr_pos);
    if (map.count(hash) == 1) {
        std::vector<Particle *> *matches = map.at(hash);
        for (Particle *m : *matches) {
            if (m != &p) {
              glm::dvec3 p_pos = p.curr_pos;
              glm::dvec3 m_pos = m->curr_pos;
              double dist = glm::length(p_pos - m_pos);
              if (dist <= h) {
                p.neighbors.push_back(*m);
              }
            }
        }
    }
}

// bin each particle's position
void Particles::build_spatial_map() {
    map.clear();
    for (auto &p : particles) {
        float hash = hash_position(p.curr_pos);
        if (map.count(hash) == 1) { // already in map
            std::vector<Particle *> *p_vec = map[hash];
            p_vec->push_back(&p);
        } else {
            std::vector<Particle *> *p_vec = new std::vector<Particle *>();
            p_vec->push_back(&p);
            std::pair<float, std::vector<Particle*> *> record(hash, p_vec);
            map.insert(record);
        }
    }
}

float Particles::hash_position(glm::dvec3 pos) {
    // find width and height
    // TODO: make radius, d, nx, ny, nz member variables
    float radius = 0.04;
    float d = 0.1;
    double width = (2 * radius * nx) + (d * (nx - 1));
    double height = (2 * radius * ny) + (d * (ny - 1));

    // partition space into 3d boxes with dimensions w * h * t
    double w = 3 * width / nx;
    double h = 3 * height / ny;
    double t = fmax(w, h);

    glm::dvec3 truncatedPos = glm::dvec3(floor(pos.x / w), floor(pos.y / h), floor(pos.z / t));
    return sqrt(2) * truncatedPos.x + sqrt(3) * truncatedPos.y + sqrt(5) * truncatedPos.z;
}

void Particles::render() const
{
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 }; // Larger number == lest shiny
    GLfloat light_position[] = { 10.0, 10.0, 10.0, 0.0 };
    glShadeModel (GL_SMOOTH);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glColor3f(0.2, 0.5, 0.8);
    glColorMaterial(GL_FRONT, GL_SPECULAR);
    glColor3f(0.9, 0.9, 0.9);
    glColorMaterial(GL_FRONT, GL_AMBIENT);
    glColor3f(0.2, 0.5, 0.8);
    
    int i = 0 ;
    for(const Particle &par : particles)
    {    
        // Source for push/pop: http://www.swiftless.com/tutorials/opengl/pop_and_push_matrices.html

        glPushMatrix(); // Set where to start the current object transformations
        glTranslatef(par.curr_pos.x, par.curr_pos.y, par.curr_pos.z); // Multiply current matrix by a translation matrix
        glutSolidSphere(0.04, 10, 10); // Render a solid sphere
        glPopMatrix(); // End the current object transformations
        i++;
    }
    
    glPopAttrib();
}

