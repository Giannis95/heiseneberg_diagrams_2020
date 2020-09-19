#include <iostream>
#include "constants.h"
#include "functions.h"
#include "global.h"
#include "vec3.h"
#include <vector>
#include <array>
#include <functional>
#include <algorithm>
#include <string>

using namespace std;
using Vec2 = std::array<double, 2>;

int main(int argc, char *argv[]) {
  std::random_device r;
  global::generator.seed(r());
  double Temperature = 10.0;
  Vec3 magnetisation = {0.0, 0.0, 0.0};
  Vec3 B_External_Field = {0.0, 0.0, 0.0};
  vector<Vec3> spins(num_spins, Vec3{0.0, 0.0, 1.0}); // Vec 3 {x , y , z }

  double amplitude = gaussian_amplitude;
  //This 4 double variables are used for the extra region in the boundary condition for 1D gaussians.
  double lower_vector_index =
	  0.0; //lower index of the physical region [-1,1] of the sample_points and discrete_potential vector (Index of element -1)
  double upper_vector_index = 0.0;//Upper Index . (Index of the element 1)
  double lower_vector_index_m_perpendicular =
	  0.0; //lower index of the physical region [0,1] of the sample_points and discrete_potential vector (Index of element 0)
  double upper_vector_index_m_perpendicular = 0.0;//Upper Index . (Index of the element 1)
  physical_region_indices(sample_points, lower_vector_index = -1, upper_vector_index = 1);
  physical_region_indices(sample_points_m_perpendicular,lower_vector_index_m_perpendicular = 0,
						  upper_vector_index_m_perpendicular =1); // TODO: there is a bug in the linear_space function (goes up to 0.99 instead of 1 need to fix then change the upper limit)

  vector<vector<double>> potential_2D(sample_points.size(), vector<double>(sample_points_m_perpendicular.size(), 0));
  const auto neighbour_list = generate_neighbour_list(L_x, L_y, L_z);
  //random_unit_vector_initialisation(spins);

  for (int temp_arguments = 2; temp_arguments < argc; ++temp_arguments) {
	bool flag_meta = false;
	string flag_check("meta");

	if (argv[1] == flag_check) {
	  flag_meta = true;
	}

	Temperature = stod(argv[temp_arguments]);
	uniform_spin_initialisation(spins);
	magnetisation = total_magnetisation(spins);
	double sigma = 0.0;
	sigma_value(sigma, Temperature, neighbour_list, B_External_Field);
	int integer_temperature = (int)Temperature;
	string filename1 = "potential_landscape_" + to_string(integer_temperature) + ".csv";
	string filename2 = "sigma_value_" + to_string(integer_temperature) + ".csv";
	string filename3 = "energy_barrier_" + to_string(integer_temperature) + ".csv";
	string filename4 = "gaussian_height_evolution_" + to_string(integer_temperature) + ".csv";
	string filename5_heisenberg = "magnetisation_values_" + to_string(integer_temperature) + ".csv";
	ofstream potential_file;
	ofstream info_file;
	ofstream energy_barrier_file;
	ofstream gaussian_height_file;
	ofstream magnetisation_file_heisenberg;

	if (!flag_meta) {
	  magnetisation_file_heisenberg.open(filename5_heisenberg);
	} else {
	  potential_file.open(filename1);
	  energy_barrier_file.open(filename3);
	  gaussian_height_file.open(filename4);
	}

//	info_file.open(filename2);
//	info_file << "Temperature: " << Temperature << " Sigma Value: " << sigma << endl;

	auto Random_Spin_Index = std::bind(uniform_int_distribution<int>(0, num_spins - 1), global::generator);
	double moves_accepted = 0.0;
	double m_perpendicular = 0.0;
	int iterations = 50000000; // Monte Carlo Steps for metadynamics (but they are too many for a heisenberg model);

	if (!flag_meta) {
	  iterations = 10000;
	}
	for (int range=0; range <100;++range) {
	  Temperature+=10;
	  uniform_spin_initialisation(spins);
	  //random_unit_vector_initialisation(spins);
	  magnetisation = total_magnetisation(spins);
	  //double sigma = 0.0;
	  sigma_value(sigma, Temperature, neighbour_list, B_External_Field);

	  for (auto i = 0; i < iterations; ++i) {
		simulation_progress_status(i, iterations, Temperature, flag_meta);
		magnetisation = total_magnetisation(spins);
		if (i > 1000 && i % 1000 == 0 && flag_meta) {
		  m_perpendicular = mz_perpendicular(magnetisation);
		  amplitude = amplitude_tempering_2d(potential_2D, magnetisation[2] / num_spins, m_perpendicular);
		  insert_2D_gaussian(magnetisation[2] / num_spins, m_perpendicular, amplitude, sample_points,
							 sample_points_m_perpendicular, potential_2D); //every t MCS add a 2D gaussian
		}
		monte_carlo_step(flag_meta, spins, magnetisation, sigma, neighbour_list, B_External_Field,
						 Temperature, moves_accepted, potential_2D);

		/* if (i % 100000 == 0) {
		   if (flag_meta) {
			 metadynamics_values_print(potential_2D, i, amplitude, energy_barrier_file, gaussian_height_file);
		   } else {
			 plain_heisenberg_values_print(i, magnetisation, magnetisation_file_heisenberg);
		   }
		 }*/
	  }
	  cout << magnetisation[2] / num_spins << " " << magnetisation[1] / num_spins << " " << magnetisation[0] / num_spins
		   << endl;
	  magnetisation_file_heisenberg<< Temperature<<"," <<magnetisation[0] / num_spins << "," << magnetisation[1] / num_spins << "," << magnetisation[2] / num_spins<< endl;
	  cout << range << " out of : " << "63" <<endl;
	}
	if (flag_meta) {
	  potential_2d_file_print(potential_2D, potential_file);
	}

  }
  return EXIT_SUCCESS;
}
