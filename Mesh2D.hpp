//Mesh class:
//Initialises a 1d vector of primitive variables,
//Calculates the boundary conditions
//Calculates the adequate dt for each time step
//HLLC flux calculator
//WAF flux calculator

///// 

#ifndef MESH_H
#define MESH_H

#include <cstdio>
#include <vector>
#include <string>
#include "Euler2D.hpp"
#include <blitz/array.h>
#include "MLS_2D.hpp"

class Mesh{
public:
  //Pointer to Euler class
  Euler * ptr_euler;
  
  //Parameters
  //std::vector<double> axis;
  int ncells; //The domain is a ncellsxncells matrix

  double x_min;
  double x_max;
  double dx;
  double y_min;
  double y_max;
  double dy;
  double time;
  double cfl;
  std::string BC;
  //Boundary functions
  Euler::W_state (*boundary1)(Euler::W_state  w);
  Euler::W_state (*boundary2)(Euler::W_state  w);
  int nGhost;

   //Data
  blitz::Array<Euler::U_state,2> Bdata;
  blitz::Array<MLS,2> MLS_data;
  blitz::Array<double,1> xaxis;
  blitz::Array<double,1> yaxis;

  //Constructor
  Mesh();
  Mesh(int ncells, double x_min, double x_max,double y_min,double y_max,double cfl, Euler::U_state (*f)(double x,double y),Euler::W_state (*b1)(Euler::W_state w), Euler::W_state (*b2)(Euler::W_state w), int nGhost,std::string arg_BC, MLS (*level_set)(double x, double y));
  ~Mesh();
  //print to screen and print to file
  void print()const;
  void save_u_state(std::string name)const;
  void save_w_state(std::string name)const;
  void slice_x_axis(std::string name)const;
  void slice_y_axis(std::string name)const;
 
  //Apply BCs
  void applyBC();
  //Apply Ghost_Fluid_conditions
  void applyGhost(double obj_speed_x,double obj_speed_y);
  // void update_level_set(double dt,double time, double obj_speed_x, double obj_speed_y,double obj_x_c,double obj_y_c, double (*level_set)(double x, double y, double x_c, double y_c));
  void advect_level_set(double dt, double time, double obj_speed_x,double obj_speed_y);
  //Calculate dt 
  double Calculate_dt();
  
  void reset(Euler::U_state (*f)(double x,double y));
  void reset(Euler::U_state (*f)(double x,double y),MLS (*level_set)(double x, double y));
  void reset_BC(Euler::W_state (*b1)(Euler::W_state w),Euler::W_state (*b2)(Euler::W_state w),std::string arg_BC); 
  
};

void flux_and_update(Mesh &m,double dt,std::string sweep_order);


blitz::Array<Euler::U_state,1> HLLC_U_state(Euler::U_state U_state_L, Euler::U_state U_state_R);

blitz::Array<Euler::U_state,1> WAF_1D(blitz::Array<Euler::U_state,1> input_data, double dt, double ds, double ncells, double nGhost,std::string limiter,std::string sweep);

double minmod(double r, double c);
double superbee(double r, double c);
int sign(double c);
#endif
