// Axis aligned bounding box

#include "bbox.h"

// Returns true if p will collide with box in given
// time step and adjusts velocity vector accordingly
bool BBox::collides(Particle p) {
    glm::dvec3 next_step = p.v + p.p;
    bool c = false;

    // Collided with bounding box
    if (next_step[0] > b_max[0] || next_step[0] < b_min[0] || 
        next_step[1] > b_max[1] || next_step[1] < b_min[1] ||
        next_step[2] > b_max[2] || next_step[2] < b_min[2]) {

        c = true;
        p.v *= -1; // Send the ball going the opposite direction (idt this is correct)
        p.v *= 0.9; // TODO: Keep track of dampening term better.
    }
    return c;
}