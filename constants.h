#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace metad {
    constexpr double kBoltzmann = 1.38064852e-23;   // J K^-1  || NIST 2014 CODATA
    const double kPi = 3.14159265358979323846264338327950288;
    const double kBohrMagneton = 9.274009994E-24;  // J T^-1  || NIST 2014 CODATA

}
//const int iterations = 50000000; // Monte Carlo Steps
const int L_x = 10;
const int L_y = 10;
const int L_z = 10;
const double J = 3e-21; // Exchange Constant Joules of Fe
const double Uniaxial_Anisotropy_k = 1.4e-23; // Uniaxial anisotropy constant order of 2 or so smaller than the exchange energy;
const int num_spins = L_x * L_y * L_z; //  Total Number of Spins

#include "functions.h"
//METADYNAMICS
const double gaussian_amplitude = 1e-24; // Height of the gaussian 1.0e-24 same as the paper
const double gaussian_width = 0.03; // Width of the gaussian 1.4e-2 same as the paper
const std::vector<double> sample_points = linear_space(-1.0, 1.0,0.01); // Used to store the predefined gaussian possition in the [-1,1] range.
const std::vector<double> sample_points_m_perpendicular = linear_space (0,1,0.01); //bug fixed.

//TEMPERED PARAMETERS
const double temp_bias = 200;//\Delta T ,

//file_names

#endif