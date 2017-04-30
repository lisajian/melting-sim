/***   NOTES   ***/

// assign sphere radius to a variable called radius

// put in particles.h:
// unordered_map<float, vector<PointMass *> *> map;

// put in particle.h:
// vector<Particle *> neighbors;

// need a constant p0 for rest density

// should go after first 'for' loop
void Particles::step() { 
  // build map to be used for spatial hashing
  build_spatial_map();
  // for each particle: 
  for (auto &p : particles) {
    find_neighboring(d, p); // d = support

    // for # iterations n:
      double p_i = 0;
      double w_p = 0;
      double w_s = 0;
      // double grad_c = 0;
      // for each of neighboring particles
      for (auto &n : p.neighbors) {
        // calculate p_i (given by SPH density estimator)
        double r = p.curr_pos - n.curr_pos;
        double r2 = (r.x * r.x + r.y * r.y + r.z * r.z);
        // if (0 <= sqrt(r2) && sqrt(r2) <= d) {
        if (r2 <= d * d) {
          // calculate poly6 kernel
          w_p = pow((d * d - r2), 3);
          w_p *= 315 / (64 * PI * pow(d, 9)); // This can be moved out of the for loop
          // calculate spiky kernel 
          w_s = r / (sqrt(r2)) * pow((d - sqrt(r2)), 2);
          w_s *= -45 / (PI * pow(d, 6));
        } else {
          w_p = 0;
          w_s = 0;
        }
        p_i += n.mass * w;

        // gradient of constraint function
        // grad_c += 1 / p0 * w_s;
        // TODO: Gradient of constraint function will be a vector
      }
      // constraint function C_i(p1,...pn):
      double c = (p_i / p0) - 1;

      // calculate lambda_i
      double lambda_i = -c / () // TODO: how to calculate grad C??
  }
}

void Particles::find_neighboring(double h, Particle &p) {
  p.neighbors.clear();
  float hash = hash_position(p.curr_pos);
  vector<Particle *> *matches = map.at(hash);
  for (Particle *m : *matches) {
    if (m != &p) {
      glm::dvec3 p_pos = p.curr_pos;
      glm::dvec3 m_pos = m.curr_pos;
      double dist = glm::length(p_pos - m_pos);
      if (dist <= 1.5 * h) {
        p.neighbors.push_back(m);
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
      vector<Particle *> *p_vec = map[hash];
      p_vec->push_back(&p);
    } else {
      vector<Particle *> *p_vec = new vector<Particle *>();
      p_vec->push_back(&p);
      std::pair<float, vector<Particle*> *> record(hash, p_vec);
      map.insert(record);
    }
  }
}

/**  find neighboring particles **/
float Particles::hash_position(glm::dvec3 pos) {
  // find width and height
  double width = (2 * radius * nx) + (d * (nx - 1));
  double height = (2 * radius * ny) + (d * (ny - 1));

  // partition space into 3d boxes with dimensions w * h * t
  double w = 3 * width / nx;
  double h = 3 * height / ny;
  double t = max(w, h);

  glm::dvec3 truncatedPos = Vector3D(floor(pos.x / w), floor(pos.y / h), floor(pos.z / t));
}