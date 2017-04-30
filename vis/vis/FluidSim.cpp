#include "FluidSim.h"

// Should just call particles.reset()
void FluidSim::reset() {
    all_particles.reset();
}

// Calls step particles.step() with appropriate params
void FluidSim::step() {
    all_particles.step(dt, h, rho, eps, k, n, del_q, solverIterations);
}

void FluidSim::render() {
    all_particles.render();
}