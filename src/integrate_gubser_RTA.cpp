#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "coord_functions.cpp"
#include "soln_functions.cpp"
#include "Parameter.h"
#include "FileIO.cpp"

int main(void)
{
  printf("This code is designed to find the Boltzmann distribution function f(tau, x, y; |p|, phi_p, v_z)... \n");
  printf("...given by the solution to the RTA Boltzmann EQN subject to the Gubser symmetry. \n");
  printf("It requires as input the hatted Temperature profile as a function of Gubser time rho : ^T(rho) \n");

  //read in the parameters
  parameters params;
  std::cout << "reading in parameters" << "\n";
  readInParameters(params);

  int n_rho = params.n_rho; //grid size for rho
  double c = params.c; // 5 eta / s constant
  double rho0 = params.rho0;
  double T_hat0 = params.T_hat0;
  int n_grid_rho = params.n_grid_rho;
  double rho_min = params.rho_min;
  double delta_rho = params.delta_rho;
  int n_grid_p_omega_hat = params.n_grid_p_omega_hat;
  double p_omega_hat_min = params.p_omega_hat_min;
  double delta_p_omega_hat = params.delta_p_omega_hat;
  int n_grid_p_eta_hat = params.n_grid_p_eta_hat;
  double p_eta_hat_min = params.p_eta_hat_min;
  double delta_p_eta_hat = params.delta_p_eta_hat;

  std::cout << "allocating space "<< "\n";
  double rho[n_rho]; //gubser time rho
  double e_rho[n_rho]; //real part of hatted energy density
  double T_rho[n_rho];
  double dummy;

  std::cout << "Reading ^e profile from input/energy.dat "<< "\n";
  std::ostringstream edat_stream;
  edat_stream << "input/energy.dat";
  std::ifstream edat(edat_stream.str().c_str());
  for (int i = 0; i < n_rho; i++)
  {
    // file format is :  rho Re{e} Im{e}
    edat >> rho[i] >> e_rho[i] >> dummy;
    //std::cout << rho[i] << "\t" << e_rho[i] << "\n";
  }
  //now find the hatted Temperature
  for (int i = 0; i < n_rho; i++) T_rho[i] = pow(M_PI * M_PI * e_rho[i] / 3., 0.25);

  //now create a gsl interp of T_hat(rho)
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *T_spline = gsl_spline_alloc(gsl_interp_cspline, n_rho);
  gsl_spline_init(T_spline, rho, T_rho, n_rho);

  //calculate the distribution function for various values of rho, p_omega_hat and p_eta_hat
  double f_soln_gubser[n_grid_rho][n_grid_p_omega_hat][n_grid_p_eta_hat];
  #pragma omp parallel for collapse(3)
  for (int irho = 0; irho < n_grid_rho; irho++)
  {
    double rho = rho_min + (double)(irho) * delta_rho;
    std::cout << "rho = " << rho << "\n";
    for (int ip_omega_hat = 0; ip_omega_hat < n_grid_p_omega_hat; ip_omega_hat++)
    {
      double p_omega_hat = p_omega_hat_min + (double)(ip_omega_hat) * delta_p_omega_hat;
      for (int ip_eta_hat = 0; ip_eta_hat < n_grid_p_eta_hat; ip_eta_hat++)
      {
        double p_eta_hat = p_eta_hat_min + (double)(ip_eta_hat) * delta_p_eta_hat;
        f_soln_gubser[irho][ip_omega_hat][ip_eta_hat] = f_solution(rho, p_omega_hat, p_eta_hat, T_spline, acc, params);
      } // for (int ip_eta_hat = 0; ...
    } // for (int ip_omega_hat = 0; ...
  } // for (int irho = 0; ...

  //now dump the solution to file
  std::ofstream myfile;
  myfile.open("output/f_rho_pomegahat_petahat.dat");
  for (int irho = 0; irho < n_grid_rho; irho++)
  {
    double rho = rho_min + (double)(irho) * delta_rho;
    for (int ip_omega_hat = 0; ip_omega_hat < n_grid_p_omega_hat; ip_omega_hat++)
    {
      double p_omega_hat = p_omega_hat_min + (double)(ip_omega_hat) * delta_p_omega_hat;
      for (int ip_eta_hat = 0; ip_eta_hat < n_grid_p_eta_hat; ip_eta_hat++)
      {
        double p_eta_hat = p_eta_hat_min + (double)(ip_eta_hat) * delta_p_eta_hat;
        myfile << rho << " " << p_omega_hat << " " << p_eta_hat << " " << f_soln_gubser[irho][ip_omega_hat][ip_eta_hat] << "\n";
      } // for (int ip_eta_hat = 0; ...
    } // for (int ip_omega_hat = 0; ...
  } // for (int irho = 0; ...
  myfile.close();

  std::cout << "Done generating f(rho, p_omega_hat, p_eta_hat) "<< "\n";

  gsl_spline_free(T_spline);
  gsl_interp_accel_free(acc);

}
