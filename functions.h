#include <vector>
#include <array>
#include "vec3.h"
#include <fstream>
#include <iostream>
#include <functional>
#include "global.h"
#include "vec3.h"
#include <vector>
#include <random>
#include <iostream>



using NeighbourList = std::vector<std::vector<int>>;

NeighbourList generate_neighbour_list(const int Nx, const int Ny,const int Nz);

double gaussian(const double&x, const double& center, const double& amplitude, const double& width);

double calculate_energy_difference(const std::vector<double>& potential);

double one_spin_energy(const int index, const std::vector<Vec3> &spins, const std::vector<std::vector<int>>& nbr_list, const Vec3 External_Field);

void display_vector (std::vector<Vec3>& vy);

void uniform_spin_initialisation(std::vector<Vec3>&spins_uniform);

void random_unit_vector_initialisation(std::vector<Vec3>& spins);
Vec3 small_angle_trial_move(const Vec3& spin, const double& sigma);

bool metropolis_algorithm(double &delta_energy, const double &T);

// M E T A D Y N A M I C S //

Vec3 total_magnetisation(const std::vector<Vec3> &spins);

double  m_z_difference_after_trial_move (std::vector<Vec3> &m_z_trial, double &m_z_current, int spin_selected, Vec3 &initial_spin );

bool floats_are_equal(const double& x, const double& y, double epsilon);

double linear_interpolation(const double& x,const double& x_lower, const double& x_upper,const double& y_lower, const double& y_upper);

double interpolated_potential(const std::vector<double>& sample_points, const std::vector<double>& discrete_potential, const double& value);

void insert_gaussian(const double& center, const double& amplitude, const double& width,const std::vector<double>& sample_points, std::vector<double>& discrete_potential);

std::vector<double> linear_space(const double& min, const double& max, const double& step);

double amplitude_tempering(const std::vector<double>& discrete_potential,const double& center,const std::vector<double>& sample_points);

void physical_region_indices (const std::vector<double>& points, double &lower_limit, double &upper_limit);

void sigma_value(double &sigma, double &sim_temp,const std::vector<std::vector<int>>& neighbour_list, const Vec3 B_External_Field);

// 2D Gaussian Metadynamics //
double gaussian_2D(const double&x, const double&y, const double& x0, const double &y0, const double& amplitude);

void insert_2D_gaussian(const double& x, const double& y , const double& amplitude,
						const std::vector<double>& sample_points, const std::vector<double>& sample_points_m_perpendicular, std::vector<std::vector<double>>& potential_2D);

double interpolated_2D_potential(const std::vector<std::vector<double>>& discrete_potential, const double& m,const double m_p);

double linear_interpolation_2(const double& x,const double& x_lower, const double& x_upper,const double& Q1, const double& Q2);

double energy_difference_2d (const std::vector<std::vector<double>>&potential);

double amplitude_tempering_2d(const std::vector<std::vector<double>>& potential, const double m, const double m_p);

// Print Functions //
void potential_2d_file_print(const std::vector<std::vector<double>> &potential_2D_values,std::ofstream &potential_2d_file);
void metadynamics_values_print(const std::vector<std::vector<double>> &potential_2d_values,const int &iteration,const double &calculated_gaussian_amplitude,
							   std::ofstream & energy_barrier_output_file,std::ofstream &gaussian_height_output_file);

void plain_heisenberg_values_print(const int &iteration, const Vec3 &magnetisation, std::ofstream &magnetisation_values_file);

void monte_carlo_step(const bool flag_metadynamics,std::vector<Vec3> &spin_vector, Vec3 magnetisation, const double sigma_value, const std::vector<std::vector<int>> &neighbour_list,
					  const Vec3 external_field_b, const double sim_temperature, double moves_accepted,const std::vector<std::vector<double>> &potential_2d);

void simulation_progress_status(int iteration, int total_iterations, double temp,bool flag_metad);