#include <functional>
#include "constants.h"
#include "functions.h"
#include "global.h"
#include "vec3.h"
#include <vector>
#include <random>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;
void print_vector (const vector<Vec3>& vy)
{
    for (auto i = 0; i < vy.size(); i++){// loops through each row of vy
        cout << endl;
      for (auto j = 0; j < vy[i].size(); j++) // loops through each element of each row
            cout << " " << vy[i][j] << " ";// prints the jth element of the ith row
            cout << endl;

    }
}

void uniform_spin_initialisation(vector<Vec3>&spins_uniform)
{
    for (auto i = 0; i <spins_uniform.size(); ++i)
    {
        spins_uniform[i] = {0.0, 0.0 ,1};

    }

}

void random_unit_vector_initialisation(vector<Vec3> &spins){
    for (auto i = 0; i < num_spins; ++i){
        spins[i]=random_unit_vector();
    }
}

NeighbourList generate_neighbour_list(const int Nx, const int Ny, const int Nz ) {
    vector<vector<vector<int>>> lattice(Nx, vector<vector<int>>(Ny, vector<int>(Nz)));
    int counter = 0;// number the lattice
    for (auto k = 0; k < Nz; ++k) {
        for (auto i = 0; i < Nx; ++i) {
            for (auto j = 0; j < Ny; ++j) {
                lattice[i][j][k] = counter;
                counter++;
            }
        }
    }

    NeighbourList neighbour_list(Nx * Ny * Nz); // create a double vector to pass the neighbour list

    counter = 0;
    for (auto k = 0; k < Nz; ++k) {
        for (auto i = 0; i < Nx; ++i) {
            for (auto j = 0; j < Ny; ++j) {
                std::vector<int> nbrs;

                //            assert(((i+1) + Nx)%Nx < Nx && ((i+1) + Nx)%Nx > 0);
                //            assert(((i-1) + Nx)%Nx < Nx && ((i-1) + Nx)%Nx > 0);
                //            assert(((j+1) + Ny)%Ny < Ny && ((j+1) + Ny)%Ny > 0);
                //            assert(((j-1) + Ny)%Ny < Ny && ((j-1) + Ny)%Ny > 0);

                nbrs.push_back(lattice[((i + 1) + Nx) % Nx][j][k]); // Down
                nbrs.push_back(lattice[((i - 1) + Nx) % Nx][j][k]); // Up
                nbrs.push_back(lattice[i][((j + 1) + Ny) % Ny][k]); // Right
                nbrs.push_back(lattice[i][((j - 1) + Ny) % Ny][k]); //  Left
                nbrs.push_back(lattice[i][j][((k + 1) + Nz) % Nz]); // Front
                nbrs.push_back(lattice[i][j][((k - 1) + Nz) % Nz]); // Back

                neighbour_list[counter] = nbrs;
                counter++;
            }
        }
    }

    return neighbour_list;
    /*for (int i=0; i<neighbour_list.size(); ++i) {
       cout << "Spin " <<i << ": ";
       for (int j = 0; j < neighbour_list[i].size(); ++j) {
           cout << neighbour_list[i][j] << " ";
       }
       cout << endl;
} */
}

double one_spin_energy(const int i, const std::vector<Vec3> &spins, const std::vector<std::vector<int>>& nbr_list, const Vec3 External_Field) { // calculate the energy of a single spin
   using namespace metad;

   assert(i >=0 && i<num_spins);

   const Vec3 anisotropy_axis = {0, 0, 1};

   double exchange_energy =0.0;
   for (const auto &j : nbr_list[i]) {
       assert(j >=0 && j < num_spins);
     exchange_energy += -2 * J * dot(spins[i], spins[j]);
   }

   double zeeman_energy = -kBohrMagneton * dot(External_Field, spins[i]);
   //double  uniaxial_anisotropy_energy = - Uniaxial_Anisotropy_k * (dot(spins[i], anisotropy_axis) * dot(spins[i], anisotropy_axis));
    double  uniaxial_anisotropy_energy = - Uniaxial_Anisotropy_k * pow(dot(spins[i], anisotropy_axis), 2) ;

   return exchange_energy + zeeman_energy + uniaxial_anisotropy_energy;
}

Vec3 small_angle_trial_move(const Vec3& spin, const double& sigma) {
  // make a small deviation of the spin by adding a fraction of random vector
  // (uniform on sphere) and normalising the result
  auto trial_spin = random_unit_vector();
  for (auto n : {0, 1, 2}) {
    trial_spin[n] = spin[n] + sigma * trial_spin[n];
  }
  normalize(trial_spin);
  return trial_spin;
}


// return true if move was accepted
bool metropolis_algorithm(double &delta_energy, const double &T) {
    using namespace metad;

    if (delta_energy < 0.0) {
      return true;
    }
   uniform_real_distribution<double> distribution(0,1);

   double alpha = distribution(global::generator);
   if (alpha > exp(-delta_energy / (kBoltzmann * T))) {
       return false;
   }
   return true;
}

Vec3 total_magnetisation(const vector<Vec3> &spins)
{
    Vec3 m={0, 0, 0};
    for (auto i =0; i <spins.size(); ++i) {
        for (auto n=0; n<3; ++n){
            m[n] +=spins[i][n];
        }
    }
    return m;
}

double  m_z_difference_after_trial_move (vector<Vec3> &m_z_trial, double &m_z_current, int spin_selected, Vec3 &initial_spin)
{
  cout << ((m_z_current* num_spins ) - (m_z_trial[spin_selected][2]-initial_spin[2] ))/num_spins;
    return ((m_z_current* num_spins ) - (m_z_trial[spin_selected][2]-initial_spin[2] ))/num_spins;
}

//////////////////////////////////////////////////////////
//   M E T A D Y N A M I C S   F U N C T I O N S        //
//////////////////////////////////////////////////////////

bool floats_are_equal(const double& x, const double& y, const double epsilon = 1e-8) {
  return abs(x - y) < epsilon;
}

// creates a vector of points between min and max with given step size
vector<double> linear_space(const double& min, const double& max, const double& step) {
  assert(min < max);
  vector<double> space;
  double value = min;
  do {
	space.push_back(value);
	value += step;
  } while (value < max+step);

  return space;
}


double linear_interpolation(const double& x,
							const double& x_lower, const double& x_upper,
							const double& y_lower, const double& y_upper) {
  assert(x_lower < x_upper);
  assert(x > x_lower || floats_are_equal(x, x_lower));
  assert(x < x_upper || floats_are_equal(x, x_upper));
auto a =y_lower + (x - x_lower)*(y_upper - y_lower) / (x_upper - x_lower);

  return y_lower + (x - x_lower)*(y_upper - y_lower) / (x_upper - x_lower);
}

double linear_interpolation_2(const double& x,
							  const double& x_lower, const double& x_upper,
							  const double& Q1, const double& Q2) {
  cout<< "linear space x_lower: " <<x_lower<<" x_upper: "<<x_upper<<endl;
  assert(x_lower < x_upper);
  assert(x > x_lower || floats_are_equal(x, x_lower));
//		assert(x < x_upper || floats_are_equal(x, x_upper));
  auto a = (x_lower-x)/(x_upper-x_lower);
  auto b =(x-x_lower)/(x_upper-x_lower);
  return a*Q1 + b*Q2;
}

double interpolated_potential(const vector<double>& sample_points, const vector<double>& discrete_potential, const double& value) {
  assert(is_sorted(begin(sample_points), end(sample_points)));
  assert(value > sample_points.front() || floats_are_equal(sample_points.front(), value));
  assert(value < sample_points.back() || floats_are_equal(sample_points.back(), value));
  assert(sample_points.size() == discrete_potential.size());
  // TODO: Test if this gives the correct points for the interpolation
  auto lower = floor((value - sample_points[0]) / (sample_points[1] - sample_points[0]));
  auto upper = lower+1;
  assert(lower < upper);
  cout << "Indices Lower:" << lower <<endl;
  return linear_interpolation(value, sample_points[lower], sample_points[upper],
							  discrete_potential[lower], discrete_potential[upper]);
}

double amplitude_tempering(const vector<double>& discrete_potential,const double& center,const vector<double>& sample_points)
{
  return gaussian_amplitude*exp(-interpolated_potential(sample_points,discrete_potential,center)/(metad::kBoltzmann*temp_bias));
}

double gaussian(const double&x, const double& center, const double& amplitude, const double& width) {
  return amplitude*exp(- (x-center)*(x-center ) / (2.0 * width*width));
}

void insert_gaussian(const double& center, const double& amplitude, const double& width,
					 const vector<double>& sample_points, vector<double>& discrete_potential) {

  assert(sample_points.size() == discrete_potential.size());
  for (auto i = 0; i < sample_points.size(); ++i) {
	  discrete_potential[i] += gaussian(sample_points[i], center, amplitude, width);
  }

  // calculate the center position for a gaussian according to mirror boundary conditions
  double mirrored_center;
  if (center >=0) {
    mirrored_center = 2 - center;
  } else {
    mirrored_center = -2 - center;
  }
  assert(mirrored_center >= -2 && mirrored_center <= 2);

  // insert the mirrored gaussian
  for (auto i = 0; i < sample_points.size(); ++i) {
    discrete_potential[i] += gaussian(sample_points[i], mirrored_center, amplitude, width);
  }
}

void physical_region_indices (const vector<double>& points, double &lower_limit, double &upper_limit){
  //lower_limit= -1.00;
  //upper_limit= 1.00;
  for (double i = 0; i < points.size(); ++i){
	if (points[i] > lower_limit) {
	  cout << endl<<"Lower Index: " << i<<" Element at index " << points[i]<<endl;
	  lower_limit= i;
	  assert(points[i]<=lower_limit);
	  break;
	}}
  for (double ii=lower_limit; ii<points.size();++ii){
	if (floats_are_equal(upper_limit,points[ii])) {
	  cout <<  " Upper Index: " << ii<<" Element at index " << points[ii] <<endl;
	  upper_limit = ii;
	  assert(points[ii]<=upper_limit);
	  break;
	}}

}


void sigma_value(double &sigma,  double &sim_temp,const std::vector<std::vector<int>>& neighbour_list, const Vec3 B_External_Field){
  auto Random_Spin_Index = bind(uniform_int_distribution<int>(0, num_spins - 1), global::generator);
  double moves_accepted_fun = 0;
  //sigma = 0.0; //since I have an increasing temperature I used the previous sigma value to start itterating, to increase the efficiency.

  	while (moves_accepted_fun/ double( double(num_spins))<=50 || moves_accepted_fun/ double(double(num_spins))>60){
    sigma +=0.01;
  	Vec3 magnetisation_fun = {0.0, 0.0, 0.0};
	vector<Vec3> spins_fun(num_spins, Vec3{0.0, 0.0, 1.0}); // Vec 3 {x , y , z } const auto neighbour_list = generate_neighbour_list(L_x, L_y, L_z);
	magnetisation_fun = total_magnetisation(spins_fun);
	moves_accepted_fun=0.0;

  for (auto i = 0; i < 100; ++i) {
	magnetisation_fun = total_magnetisation(spins_fun);
	for (auto n = 0; n < num_spins; ++n)
	{
	  auto const random_index = Random_Spin_Index();   // Choose a random spin index for the trial move
	  auto initial_spin = spins_fun[random_index];        // Save the initial spin before replacing with the trial spin  (Vec3)
    auto initial_energy = one_spin_energy(random_index, spins_fun, neighbour_list, B_External_Field);
   auto trial_spin = small_angle_trial_move(spins_fun[random_index], sigma);
   spins_fun[random_index] = trial_spin;

auto trial_energy = one_spin_energy(random_index, spins_fun, neighbour_list, B_External_Field);

Vec3 trial_magnetisaion;
for (auto n : {0, 1, 2}) {
trial_magnetisaion[n] = magnetisation_fun[n] - initial_spin[n] + trial_spin[n];
}
auto energy_difference_fun = (trial_energy  - initial_energy);
bool move_accepted_fun = metropolis_algorithm(energy_difference_fun,sim_temp);

if (move_accepted_fun) {
magnetisation_fun = trial_magnetisaion;
moves_accepted_fun++;
} else {
spins_fun[random_index] = initial_spin;
}
}
	}
	//  cout << "acceptance rate: " << moves_accepted_fun / double(double(num_spins)) << " sigma: "<< sigma <<endl;
  }
}

// Calculates the energy difference in the potential using the central 50%
// The left and right 25% are fake data for the mirrored boundaries
double calculate_energy_difference(const vector<double> &potential) {
  const auto margin = potential.size()/4;
  const double max = *max_element(potential.begin()+margin, potential.end()-margin);
  const double min = *min_element(potential.begin()+margin, potential.end()-margin);
  return max - min;
}

//*************************** 2D Gaussian Functions ************************************
double gaussian_2D(const double&x, const double&y, const double& x0, const double &y0, const double& amplitude) //x = mz, y=mz_perpendicular, x0: sample points for mz, yO:sample points for mz_p [0,1]
{
  return amplitude*exp(- (((x-x0)*(x-x0))+((y - y0)*(y - y0)))/(2.0 * gaussian_width*gaussian_width));//xo is on X
}

void insert_2D_gaussian(const double& x, const double& y , const double& amplitude,
						const vector<double>& sample_points, const vector<double>& sample_points_m_perpendicular, vector<vector<double>>& potential_2D) {
  for (int i = 0; i < sample_points.size(); ++i) { //TODO: simplify the input parameter, use only the potential 2D vector.
	for (int ii = 0; ii < sample_points_m_perpendicular.size(); ++ii) { //sample_points_m_perpendicular.size()
	  potential_2D[i][ii] += gaussian_2D(sample_points[i], sample_points_m_perpendicular[ii], x, y, amplitude);
	 // cout <<potential_2D[i][ii]<<"1";
	}}
}

//--------------------------------Bilinear Interpolation----------------------------------------------//
double interpolated_2D_potential(const vector<vector<double>>& discrete_potential, const double& m,const double m_p)
{
  auto lower_m = floor((m - sample_points[0]) / (sample_points[1] - sample_points[0])); //lower_y index for the discrete_potential
  auto upper_m = lower_m+1;

  auto lower_m_p = floor((m_p - sample_points_m_perpendicular[0]) / (sample_points_m_perpendicular[1] - sample_points_m_perpendicular[0]));//lower_x index for the discrete_potential
  auto upper_m_p = lower_m_p+1;

  assert(lower_m < upper_m);
  assert(lower_m_p < upper_m_p);
  //(corner potential values around the  point of interest) : f(x1,y1)=Q(11) , f(x1,y2)=Q(12), f(x2,y1), f(x2,y2)
  double Q11=discrete_potential[lower_m][lower_m_p];
  //cout <<"Q11 = " <<Q11<<" ,";
  double Q12=discrete_potential[lower_m][upper_m_p];
  //cout <<"Q12 = " <<Q12<<" ,";
  double Q21=discrete_potential[upper_m][lower_m_p];
 // cout <<"Q21 = " <<Q21<<" ,";
  double Q22=discrete_potential[upper_m][upper_m_p];
  //cout <<"Q22 = " <<Q22<<endl;
  auto R1 = linear_interpolation(m_p,sample_points_m_perpendicular[lower_m_p],sample_points_m_perpendicular[upper_m_p],Q11,Q21);
  auto R2 = linear_interpolation(m_p,sample_points_m_perpendicular[lower_m_p],sample_points_m_perpendicular[upper_m_p],Q12,Q22);
 // cout<<"R1,R2 : " << R1<< ", "<<R2 <<endl;
 // cout << linear_interpolation(m,sample_points[lower_m],sample_points[upper_m],R1,R2) << endl;
  return linear_interpolation(m,sample_points[lower_m],sample_points[upper_m],R1,R2); //interpolate between R1 and R2
}

double energy_difference_2d (const vector<vector<double>>&potential) {
  auto min=1.0;
  auto max=0.0;
  for (auto i=0; i<potential.size(); ++i){
    for (auto ii=0; ii<potential[i].size();++ii){
      if (potential[i][ii]>max){
        max = potential[i][ii];
      }
      if (potential[i][ii]<min){
        min =potential[i][ii];
      }

      assert(min < max || floats_are_equal(min,max));
    }
  }
  return max-min;
}

double amplitude_tempering_2d(const vector<vector<double>>& potential, const double m, const double m_p)
{
  return gaussian_amplitude*exp(-interpolated_2D_potential(potential,m,m_p)/(metad::kBoltzmann*temp_bias));
}

//####### PRINT FUNCTIONS ##########
void potential_2d_file_print(const vector<vector<double>> &potential_2D_values,ofstream &potential_2d_file){
  for (auto i=0; i<sample_points.size();++i){
	for(auto ii=0; ii<sample_points_m_perpendicular.size();++ii){
	  potential_2d_file << sample_points_m_perpendicular[ii] << "," << sample_points[i] << "," << potential_2D_values[i][ii] << endl;
	}}
}

void metadynamics_values_print(const vector<vector<double>> &potential_2d_values,const int &iteration,const double &calculated_gaussian_amplitude,
	ofstream & energy_barrier_output_file,ofstream &gaussian_height_output_file){
  auto energy_barrier_2d = energy_difference_2d(potential_2d_values);
  energy_barrier_output_file << iteration << "," << energy_barrier_2d << endl;
  gaussian_height_output_file << iteration << "," << calculated_gaussian_amplitude << endl;
}

void plain_heisenberg_values_print(const int &iteration, const Vec3 &magnetisation, ofstream &magnetisation_values_file){
  magnetisation_values_file << iteration<<","<<magnetisation[0]<<","<<magnetisation[1]<<","<<magnetisation[2]<<endl;
}

// Monte Carlo Step function (both for the plain heisenberg model and Metadynamics)
void monte_carlo_step(const bool flag_metadynamics,vector<Vec3> &spin_vector, Vec3 magnetisation, const double sigma_value, const vector<vector<int>> &neighbour_list,
					 const Vec3 external_field_b, const double sim_temperature, double moves_accepted,const vector<vector<double>> &potential_2d){
  std::random_device r;
  global::generator.seed(r());
  auto energy_difference = 0.0;
  auto initial_potential = 0.0;
  auto Random_Spin_Index = std::bind(uniform_int_distribution<int>(0, num_spins - 1), global::generator);
  for (auto n = 0; n < num_spins; ++n) //monte carlo step
  {
	auto const random_index = Random_Spin_Index();   // Choose a random spin index for the trial move
	auto initial_spin =spin_vector[random_index];      // Save the initial spin before replacing with the trial spin  (Vec3)
	auto initial_energy = one_spin_energy(random_index, spin_vector, neighbour_list, external_field_b);
	if (flag_metadynamics) {
	   auto m_perpendicular = mz_perpendicular(magnetisation); //TODO: separate for 1 and 2 cv's.
	   initial_potential = interpolated_2D_potential(potential_2d, magnetisation[2] / num_spins, m_perpendicular);
	}

	// Replace the initial spin with a random trial spin
	auto trial_spin = small_angle_trial_move(spin_vector[random_index], sigma_value);
	spin_vector[random_index] = trial_spin;
	//calculates the trial spin energy and potential.
	auto trial_energy = one_spin_energy(random_index, spin_vector, neighbour_list, external_field_b);
	Vec3 trial_magnetisaion;
	for (auto n : {0, 1, 2}) {
	  trial_magnetisaion[n] = magnetisation[n] - initial_spin[n] + trial_spin[n];
	}
	if (flag_metadynamics) {
	  auto trial_m_perpendicular = mz_perpendicular(trial_magnetisaion);
	  auto trial_potential = interpolated_2D_potential(potential_2d,trial_magnetisaion[2] / num_spins,trial_m_perpendicular); //TODO: initialised in the function different from previous version
	   energy_difference =(trial_energy + trial_potential) - (initial_energy + initial_potential); //potential only for metadynamics
	}
	if (!flag_metadynamics) {
	 energy_difference = (trial_energy) - (initial_energy);
  }
	bool move_accepted = metropolis_algorithm(energy_difference,sim_temperature);
	if (move_accepted) {
	  magnetisation = trial_magnetisaion;
	  moves_accepted++;
	} else {
	  spin_vector[random_index] = initial_spin;
	}
  }
}

void simulation_progress_status(int iteration, const int total_iterations, double temp, bool flag_metad){
  string simulation_parameter_output;
  auto simulation_steps_scale =0;
  if (iteration % 1000 == 0) {
	if (flag_metad) {
	  simulation_parameter_output = "metadynamics";
	  simulation_steps_scale =100000;
	} else {
	  simulation_parameter_output = "Plain Heisenberg";
	  simulation_steps_scale =10000;
	}
	cout <<iteration/simulation_steps_scale  << "out of: " << total_iterations / simulation_steps_scale << " Iterations for Temperature: " << temp << " -"
		 << simulation_parameter_output << endl;
  } else {
	return;
  }}
