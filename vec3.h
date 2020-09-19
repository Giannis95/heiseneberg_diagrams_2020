//
// Created by ioannis charalampidis on 25/02/2020.
//

#ifndef HEISENBERG_VEC3_H
#define HEISENBERG_VEC3_H

#include <array>
#include <cmath>
#include <random>
#include <iostream>

#include "constants.h"
#include "global.h"

using Vec3 = std::array<double, 3>;

inline double dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double length(const Vec3& a) {
    return sqrt(dot(a,a));
}

inline void normalize(Vec3& a) {
    const double l = length(a);
    for (auto i : {0,1,2}) {
        a[i] = a[i] / l;
    }
}

inline Vec3 random_unit_vector(){
  std::uniform_real_distribution<> dist;
  double v1, v2, s;
  do {
    v1 = -1.0 + 2.0 * dist(global::generator);
    v2 = -1.0 + 2.0 * dist(global::generator);
    s = (v1 * v1) + (v2 * v2);
  } while (s > 1.0);

  auto ss = sqrt(1.0 - s);

  return {2.0 * v1 * ss, 2.0 * v2 * ss, 1.0 - 2.0 * s};
}
//TODO : Test the mz_perpendicular inline function in the code-runner
inline double mz_perpendicular(const Vec3& mag){
  return (sqrt((mag[0]*mag[0])+(mag[1]*mag[1])))/num_spins;
}
#endif //HEISENBERG_VEC3_H
