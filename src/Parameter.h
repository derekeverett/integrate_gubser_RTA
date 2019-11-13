#pragma once
struct parameters
{
  //parameters for constructing solution in gubser coordinates
  int n_rho; //grid size for rho
  double c; // 5 eta / s constant
  double rho0; //initial gubse time
  double T_hat0; //initial ^T
  int n_grid_rho; //number of points to compute in rho
  double rho_min; //min value of rho
  double delta_rho; //rho spacing
  int n_grid_p_omega_hat; // number of points in p_omega_hat
  double p_omega_hat_min; //min value of p_omega_hat
  double delta_p_omega_hat; //p_omega_hat spacing
  int n_grid_p_eta_hat;
  double p_eta_hat_min;
  double delta_p_eta_hat;
  //parameters for constructing solution in milne coords
  int n_r;
  double delta_r;
  int n_p_mag;
  double delta_p_mag;
  int n_phip;
  int n_vz;
  double tau0;
  double q0;
};
