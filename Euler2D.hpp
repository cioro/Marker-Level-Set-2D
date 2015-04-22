#pragma once
#ifndef EULER2D_H
#define EULER2D_H
/*
namespace {
static const double GAMMA = 0.15;
}*/

enum VARIABLE {RHO, V, U, PRESSURE, MOMENT_U,MOMENT_V,ENERGY};
class Euler{

public:
  //Data members
  const double gamma = 1.4;
  
  struct W_state{
    double rho;
    double u;
    double v;
    double P;
    W_state();
    W_state(double rho, double u,double v, double P);
    double get_var(VARIABLE v);
    void print()const; 
 };
  
  struct U_state{
    double rho;
    double moment_u;
    double moment_v;
    double energy;
    U_state();
    U_state(double rho, double moment_u,double moment_v, double energy);
    double get_var(VARIABLE v);
    void print()const;  
};

  //Constructors; 
  Euler();

  //Do not think I need explicit constructors other than the empty one;
 // Euler(W_state& w);
 // Euler(U_state& u);

  //Member functions
  double a();//Calculate speed of sound a=sqrt(gamma*P/rho);
  // double a(U_state& u); //could not find formula to calculate speed of sound from conserved var
  double a(const W_state& w);

  double int_energy(const W_state& w);
  
  U_state flux(const U_state& u);
  U_state flux_y(const U_state& u);

  U_state CfromP(const W_state& w);
  W_state PfromC(const U_state& u);

  
};

#endif
