#ifndef MLS2D_H
#define MLS2D_H

#include<vector>
#include<initializer_list>
#include<functional>
#include<unordered_map>
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
#include"MLS_cell.hpp"
#include"MLS_particle.hpp"

namespace std {
  template <> 
  struct hash <std::pair<int,int>> {
   size_t operator () (const std::pair<int,int> & p )const {
      return p.first+1000*p.second;
    }
  };

}//end of std namespace


namespace MLS{
  class Mesh{

  public:

    int ncells;
    int nGhost;
    double cfl;
    double dx,dy;
    double x_min,x_max,y_min,y_max;
    double dt;
    double time;
    double T_max;
    int iter_counter;
    int numMarkers;
 
    double Z_area;

    blitz::Array<Cell,2> MLS_data;
    blitz::Array<double,1> xaxis;
    blitz::Array<double,1> yaxis;
    blitz::Array<double,1> x_cell_axis;
    blitz::Array<double,1> y_cell_axis;
    std::vector<Particle> MLS_markers;

    std::unordered_map<std::pair<int,int>,std::vector<std::pair<Particle,Particle>>> interface;
    std::unordered_map<std::pair<int,int>,std::vector<std::pair<Particle,Particle>>> plus_dx_LS;
    std::unordered_map<std::pair<int,int>,std::vector<std::pair<Particle,Particle>>> minus_dx_LS;

    double (*speed_x)(double x, double y, double t,double T);
    double (*speed_y)(double x, double y, double t,double T);
  
    //constructor
    Mesh();
    Mesh(int numMarkers, double T, int ncells,int nGhost, double x_min, double x_max,double y_min,double y_max,double cfl, double (*speed_x)(double x, double y, double t, double T), double (*speed_y)(double x, double y, double t,double T), Cell (*level_set)(double x, double y,double time,double T));

    ~Mesh();
    void Calculate_dt();
    void applyBC();
    void applyBC_periodic();
    void applyBC(blitz::Array<double,2> & input);
    void applyBC_periodic(blitz::Array<double,2> & input);
    void advect_level_set();
    void save_to_file(std::string name)const;
    void vtk_output(std::string filename)const;
    void vtk_output_marker(std::string filename)const;

    double D_minus(int i, int j, std::string dir,const blitz::Array<double,2> & input);
    double D_plus(int i, int j, std::string dir,const blitz::Array<double,2> & input);
    double WENO(int i, int j, std::string stencil, std::string dir,const blitz::Array<double,2> & input);
    blitz::Array<double,2> spatial_WENO(blitz::Array<double,2> input);
    blitz::Array<double,2> spatial_WENO_X(blitz::Array<double,2> input);
    blitz::Array<double,2> spatial_WENO_X_periodic(blitz::Array<double,2> input);
    blitz::Array<double,2> spatial_WENO_Y(blitz::Array<double,2> input);
    blitz::Array<double,2> spatial_WENO_Y_periodic(blitz::Array<double,2> input);
    blitz::Array<double,2> spatial_first(blitz::Array<double,2> input);
    blitz::Array<double,2> spatial_first_X(blitz::Array<double,2> input);
    blitz::Array<double,2> spatial_first_Y(blitz::Array<double,2> input);
    void advect_RK();
    void advect_RK_WENO();
    void advect_RK_WENO_periodic();
    void advect_markers();
    void RK(Particle & p);
    void correction1();
    void marchingSquares();
    int Case(double threshold, double v1, double v2, double v3, double v4);
    std::vector<std::pair<Particle,Particle>> squareReconstruction(int,int,int,double,double,double,double,double);
    double dist(Particle,Particle);
    double distance_squared(Particle,Particle);
    double dot_prod(Particle,Particle);
    double minimum_distance(Particle,Particle,Particle);
  };
}
#endif
