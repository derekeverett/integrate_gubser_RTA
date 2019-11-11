#pragma once
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "coord_functions.cpp"
#include "Parameter.h"

double f_eq(double p_omega_hat, double p_eta_hat, double rho, double T_hat)
{
  double p_rho_hat = sqrt( (p_omega_hat*p_omega_hat / (cosh(rho)*cosh(rho))) + p_eta_hat*p_eta_hat );
  return exp(-p_rho_hat / T_hat);
}

double DampingFunction(double rho2, double rho1, gsl_spline *T_spline, gsl_interp_accel *acc, parameters &params)
{
  //define an integration step size
  double drho = params.delta_rho;
  double c = params.c;
  double result = 0.;
  double rho = rho1;
  while (rho < rho2)
  {
    result += gsl_spline_eval(T_spline, rho, acc);
    rho += drho;
  }
  result *= drho;
  return exp(-1.0 * result / c);
}

double f_solution(double rho, double p_omega_hat, double p_eta_hat, gsl_spline *T_spline, gsl_interp_accel *acc, parameters &params)
{
  double drho = params.delta_rho;
  double rho0 = params.rho0;
  double T_hat0 = params.T_hat0;

  //the first term
  double I1 = DampingFunction(rho, rho0, T_spline, acc, params) * f_eq(p_omega_hat, p_eta_hat, rho0, T_hat0);
  //the second term
  double I2 = 0.;
  double rho_p = rho0;
  while (rho_p < rho)
  {
    double T_hat_rho_p = gsl_spline_eval(T_spline, rho_p, acc);
    I2 += DampingFunction(rho, rho_p, T_spline, acc, params) * T_hat_rho_p * f_eq(p_omega_hat, p_eta_hat, rho_p, T_hat_rho_p);
    rho_p += drho;
  }
  I2 *= drho;
  return I1 + I2;
}
