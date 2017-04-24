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

Particles::Particles() {
    // std::cout << "MADE" << std::endl;
    bbox = BBox(glm::dvec3(-2, -2, -2), glm::dvec3(2, 2, 2)); // TODO: remove janky default
    reset();
}

// Resets all particles back to default state.
void Particles::reset() {
    // Number of particles in each dimension
    // std::cout << "MADE" << std::endl;
    int nx = 1;
    int ny = 1;
    int nz = 1;
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
                par.forces = glm::dvec3(0.0, -9.8, 0.0); // TODO: Handle external forces better. This is gravity
                particles.push_back(par);
            }
        }
    }
}

// Single update step for all particles
void Particles::step() {
    for (auto &p : particles) {
        // Use Verlett integration
        double d = 0.05;
        p.vdt = p.curr_pos - p.last_pos;
        std::cout << "Before adjustment " << p.vdt.y << std::endl;
        bool c = bbox.collides(p);
        glm::dvec3 new_pos = p.curr_pos + (1 - d) * p.vdt + (p.forces / (double) p.mass);
        p.last_pos = p.curr_pos;
        p.curr_pos = new_pos;
        if (c && p.curr_pos.y < -2) {
            std::cout << "Not adjusted" << std::endl;
            std::cout << p.curr_pos.y << std::endl;
        }
        // p.last_pos = p.curr_pos;
        // p.curr_pos += p.v;
        // p.v += p.forces / (double) p.mass;
        // bbox.collides(p);
    }
}

// 
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

