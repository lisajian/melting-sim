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
#include "Utils.h"
#include <limits>
#include <omp.h>
#include <valarray>

const double SQ2 = sqrt(2);
const double SQ3 = sqrt(3);
const double SQ5 = sqrt(5);

Particles::Particles(double const_h) {
    bbox = BBox(glm::dvec3(-2, -2, -2), glm::dvec3(2, 2, 2));
    default_forces = {glm::dvec3(0, -9.8, 0)};
    default_mass = 1;
    nx = 5;
    ny = 5;
    nz = 5;
    /////////////////////////////////////////////////////////////
    float radius = 0.04;
    float d = 0.1;
    double width = (2 * radius * nx) + (d * (nx - 1));
    double height = (2 * radius * ny) + (d * (ny - 1));

    // partition space into 3d boxes with dimensions w * h * t
    w = 1.5 * width / nx;
    h = 1.5 * height / ny;
    t = fmax(w, h);
    /////////////////////////////////////////////////////////////

    poly6_const = 315.0 / (64.0 * M_PI * std::pow(const_h, 9.0));
    spiky_const = - 45.0 / (M_PI * std::pow(const_h, 6.0));

    reset();
}

// Resets all particles back to default state.
void Particles::reset() {
    // Number of particles in each dimension
    particles.clear();
    float x_offset = 0;
    float y_offset = 1;
    float z_offset = 0;
    float d = 0.1;
    int i = 0;
    for(int x=0; x<nx; x++)
    {
        for(int y=0; y<ny; y++)
        {
            for(int z=0; z<nz; z++)
            {
                Particle par(glm::dvec3((x+0.5-nx*0.5)*d + x_offset, (y+0.5)*d-1.0 + y_offset, (z+0.5-nz*0.5)*d + z_offset));
                par.forces = glm::dvec3(0.0, 0, 0.0);
                for (glm::dvec3 f : default_forces) {
                    par.forces += f;
                }
                par.vdt = glm::dvec3();
                par.mass = default_mass;
                par.id = i;
                i++;
                particles.push_back(par);
            }
        }
    }
}

// Single update step for all particles
void Particles::step(double dt, double h, double rho, double eps, double k, \
                     double const_n, double del_q, int solverIterations, double c, double eps_vort) {

    // Generic time step update
    // #pragma omp parallel
    for (auto &p : particles) {
        // std::cout << p.curr_pos.x << ", " << p.curr_pos.y << ", " << p.curr_pos.z << std::endl;
        p.vdt = 0.85 * (p.vdt + dt * p.forces / (double) p.mass);
        p.pred_pos = p.curr_pos + dt * p.vdt;
        p.num_collisions = 0;
        p.adjustment_vec = glm::dvec3(0, 0, 0);
        // p.curr_pos = p.pred_pos; // TODO: comment out
        // bbox.collides(p, dt); // TODO: comment out
    }

    // Get neighbors. Distance between point centers is hardcoded
    build_spatial_map();
    // #pragma omp parallel
    for (auto &p : particles) {
        find_neighboring(h, p);
        p.spiky = std::vector<glm::dvec3>(p.neighbors.size());
    }

    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Particle &p = particles.at(i);
        calc_spiky(p, spiky_const, h);
        // calc_poly6(p, poly6_const, h);
    }

    // constants for kernels
    // double poly6_const = 315.0 / (64.0 * M_PI * std::pow(h, 9.0));
    // double spiky_const = - 45.0 / (M_PI * std::pow(h, 6.0));

    // double poly6_del_q = poly6_kernel(del_q * glm::dvec3(1, 0, 0), poly6_const, h);
    double s_corr_denom = pow(h - del_q, 3.0);
    glm::dvec3 zero = glm::dvec3(0, 0, 0);

    // #pragma omp parallel
    for (int i = 0; i < solverIterations; i++) {
        // Determine lambda_i
        // for (auto &p : particles) {
        #pragma omp parallel for
        for (int k = 0; k < particles.size(); k++) {
            Particle &p = particles.at(k);
            // Use poly6 kernel for rho (pressure)
            double rho_i = 0.0;
            // Use spiky kernel for gradients
            glm::dvec3 grad_i_W = glm::dvec3(0, 0, 0);
            double grad_j_W = 0.0;
            // for (auto &n : p.neighbors) {
            for (int j = 0; j < p.neighbors.size(); j++) {
                Particle &n = p.neighbors.at(j);
                glm::dvec3 diff = p.curr_pos - n.curr_pos;
                rho_i += pow(h * h - glm::dot(diff, diff), 3.0);

                glm::dvec3 intermed = p.spiky[j];

                grad_i_W += intermed;
                grad_j_W += glm::dot(intermed, intermed);
            }

            rho_i *= default_mass * poly6_const;
            p.C = rho_i / rho - 1.0;

            p.lambda = - p.C / (grad_j_W + glm::dot(grad_i_W, grad_i_W) + eps);
        }


        // Determine delta p
        #pragma omp parallel for
        for (int k = 0; k < particles.size(); k++) {
            Particle &p = particles.at(k);
            glm::dvec3 new_del_p = glm::dvec3(0, 0, 0);
            // for (auto &n : p.neighbors) {
            for (int j = 0; j < p.neighbors.size(); j++) {
                Particle &n = p.neighbors.at(j);

                glm::dvec3 intermed = p.spiky[j];

                /////////////////////// to undo: ///////////////////////
                // double s_corr = p.poly6[j] / poly6_del_q;
                // s_corr = std::pow(s_corr, const_n) * -1 * k;

                ///////////// TRY CALC S_CORR W/ SPIKY /////////////////
                double s_corr = (h - glm::length(p.curr_pos - n.curr_pos), 3.0);
                s_corr /= s_corr_denom;
                s_corr = -k * pow(s_corr, const_n);
                ////////////////////////////////////////////////////////

                new_del_p += (p.lambda + n.lambda + s_corr) * intermed;
            }
            p.del_p = new_del_p / rho;
            ///////////////////////////////////////////////////////////
            /////////////// Update pred_pos here instead //////////////
            p.pred_pos += p.del_p;
            bbox.collides(p, dt);
            ///////////////////////////////////////////////////////////
            // TODO: Collision detection?
            // bbox.collides(p, dt);
            // handle particle collisions with each other
            // particle_collisions(p);
        }

        ///////////////////////////////////////////////////////////////
        /////////// This was moved into the prev block ////////////////
        // Update pred_pos
        // for (auto &p : particles) {
        // #pragma omp parallel for
        // for (int k = 0; k < particles.size(); k++) {
        //     Particle &p = particles.at(k);
        //     p.pred_pos += p.del_p;
        //     if (p.num_collisions != 0) {
        //         p.pred_pos += (1.0 / (double) p.num_collisions) * p.adjustment_vec;
        //     }
        //     bbox.collides(p, dt);
        //     // // handle particle collisions with each other
        //     // particle_collisions(p);
        // }
        /////////////////////////////////////////////////////////////////
    }

    // Update velocity, position, vorticity, XSPH
    // #pragma omp parallel
    #pragma omp parallel for
    for (int k = 0; k < particles.size(); k++) {
        Particle &p = particles.at(k);
        p.vdt = (1.0 / dt) * (p.pred_pos - p.curr_pos);

        // bbox.collides(p, dt);

        // TODO: Do this after?
        p.vdt.x *= p.wall_collide.x;
        p.vdt.y *= p.wall_collide.y;
        p.vdt.z *= p.wall_collide.z;

        // Apply XSPH viscocity
        // Determine eta for vorticity calculation
        glm::dvec3 visc = glm::dvec3(0, 0, 0);
        glm::dvec3 omega = glm::dvec3(0, 0, 0);
        glm::dvec3 eta = glm::dvec3(0, 0, 0);
        // for (auto &n : p.neighbors) {
        for (int j = 0; j < p.neighbors.size(); j++) {
            Particle &n = p.neighbors.at(j);

            glm::dvec3 vij = n.vdt - p.vdt;
            omega += glm::cross(vij, -p.spiky[j]);

            double r = glm::length(p.curr_pos - n.curr_pos);
            double weight = - pow(r, 3.0) / (2 * pow(h, 3.0)) + pow(r, 2.0) / pow(h, 2.0);
            weight += h / (2 * r) - 1;
            visc += weight * vij;
        }
        if (!glm::all(glm::equal(omega, zero))) {
            eta = glm::normalize(omega);
        }

        glm::dvec3 vorticity = eps_vort * glm::cross(eta, omega);

        // Update with vorticity
        p.forces += vorticity;

        // Update velocity with viscocity
        p.vdt += c * visc;

        p.curr_pos = p.pred_pos;
    }
}

void Particles::calc_spiky(Particle &p, double spiky_const, double h) {
    for (int i = 0; i < p.neighbors.size(); i++) {
        Particle &n = p.neighbors.at(i);
        p.spiky[i] = spiky_kernel_grad(p.curr_pos - n.curr_pos, spiky_const, h);
    }
}

void Particles::calc_poly6(Particle &p, double poly6_const, double h) {
    for (int i = 0; i < p.neighbors.size(); i++) {
        Particle &n = p.neighbors.at(i);
        p.poly6[i] = poly6_kernel(p.curr_pos - n.curr_pos, poly6_const, h);
    }
}

void Particles::particle_collisions(Particle &p) {
    float radius = 0.1;
    // std::cout << "=======================" << std::endl;
    // std::cout << "p.id: " << p.id << std::endl;
    for (auto &n : p.neighbors) {
        // std::cout << "neighbors with: " << n.id << std::endl;
        double diff = glm::length(p.pred_pos - n.pred_pos);
        p.particle_collide = glm::dvec3(1, 1, 1);
        if (diff < 2 * radius) {
            // std::cout << "... and collided" << n.id << std::endl;
            // std::cout << "p.pred_pos before: (" << p.pred_pos.x << ", " << p.pred_pos.y << ", " << p.pred_pos.z << ")" << std::endl;
            // std::cout << "n.pred_pos before: (" << n.pred_pos.x << ", " << n.pred_pos.y << ", " << n.pred_pos.z << ")" << std::endl;
            // apply correction vector to both
            glm::dvec3 p_corr = p.pred_pos - p.curr_pos;
            glm::dvec3 n_corr = n.pred_pos - n.curr_pos;
            double p_corr_len = glm::length(p_corr) - (diff / 2);
            double n_corr_len = glm::length(n_corr) - (diff / 2);
            p_corr = glm::normalize(p_corr);
            n_corr = glm::normalize(n_corr);
            // p.pred_pos = p.curr_pos + ((double) p_corr_len * p_corr);
            // n.pred_pos = n.curr_pos + ((double) n_corr_len * n_corr);
            p.adjustment_vec += p_corr_len * p_corr;
            n.adjustment_vec += n_corr_len * n_corr;
            p.num_collisions++;
            n.num_collisions++;
            // std::cout << "p.pred_pos after: (" << p.pred_pos.x << ", " << p.pred_pos.y << ", " << p.pred_pos.z << ")" << std::endl;
            // std::cout << "n.pred_pos after: (" << n.pred_pos.x << ", " << n.pred_pos.y << ", " << n.pred_pos.z << ")" << std::endl;
        }
    }
}


// Find all neighbors n of p such that the distance between the
// centers of n and p is at most h
void Particles::find_neighboring(double h, Particle &p) {
    p.neighbors.clear();
    // // Naive way of finding neighbors
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
              double dist = glm::length(p.curr_pos - m->curr_pos);
              if (dist <= 0.11) {
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
    // float radius = 0.04;
    // float d = 0.1;
    // double width = (2 * radius * nx) + (d * (nx - 1));
    // double height = (2 * radius * ny) + (d * (ny - 1));
    //
    // // partition space into 3d boxes with dimensions w * h * t
    // double w = 1.5 * width / nx;
    // double h = 1.5 * height / ny;
    // double t = fmax(w, h);

    // glm::dvec3 truncatedPos = glm::dvec3(floor((pos.x + 2)/ w), floor((pos.y + 2) / h), floor((pos.z + 2) / t));
    // std::cout << "pos x: " << truncatedPos.x << ", pos y: " << truncatedPos.y << ", pos z: " << truncatedPos.z << std::endl;
    float hash = SQ2 * floor((pos.x + 2)/ w) + SQ3 * floor((pos.y + 2) / h) + SQ5 * floor((pos.z + 2) / t);
    // std::cout << "hash: " << hash << std::endl;
    return hash;
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
