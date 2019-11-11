#pragma once

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

#ifdef _OPENMP
#include <omp.h>
#endif

double rho_func(double tau, double r, double q)
{
    double num = 1. - (q*q) * (tau*tau) + (q*q) * (r*r);
    double den = 2. * q * tau + 1e-10;
    return -asinh(num / den);
}

double theta_func(double tau, double r, double q)
{
  double num = 2. * q * r;
  double den = 1. + (q*q) * (tau*tau) - (q*q) * (r*r) + 1e-10;
  return atan(num / den);
}

double kappa_func(double tau, double r, double q)
{
  double num = 2. * q*q * tau * r;
  double den = 1. + q*q * ( tau*tau + r*r) + 1e-10;
  return atanh(num / den);
}

double p_theta_hat_func(double p_tau, double p_r, double q)
{
  double num = 2. * q * p_r;
  double den = 1. + (q*q) * (p_tau*p_tau) - (q*q) * (p_r*p_r) + 1e-10;
  return atan(num / den);
}

double tau_func(double rho, double theta, double q)
{
  double num = 1. / cosh(rho);
  double den = cos(theta) - tanh(rho) + 1e-10;
  return num / den / q;
}

double r_func(double rho, double theta, double q)
{
  double num = sin(theta);
  double den = cos(theta) - tanh(rho) + 1e-10;
  return num / den / q;
}

double p_omega_hat_func(double p_theta_hat, double p_phi_hat, double theta)
{
  double c = (p_phi_hat / (sin(theta) + 1e-10) );
  return p_theta_hat*p_theta_hat + c*c;
}

double p_theta_func(double rho, double theta, double q, double p_tau, double p_r)
{
  double del_tau_del_theta_num = - sin(theta) / cosh(rho);
  double c = (cos(theta) - tanh(rho));
  double den = c*c + 1e-10;
  double del_tau_del_theta = del_tau_del_theta_num / den / q;

  double del_r_del_theta_num = cos(theta)*(cos(theta) - tanh(rho)) + sin(theta)*sin(theta);
  double del_r_del_theta = del_r_del_theta_num / den / q;

  return del_tau_del_theta * p_tau + del_r_del_theta * p_r;
}
