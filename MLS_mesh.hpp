#ifndef MLS2D_H
#define MLS2D_H

#include<vector>
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
#include"MLS_cell.hpp"

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
    int iter_counter;
 

    blitz::Array<Cell,2> MLS_data;
    blitz::Array<double,1> xaxis;
    blitz::Array<double,1> yaxis;

    Cell (*speed_x)(double x, double y, double t);
    Cell (*speed_y)(double x, double y, double t);
  
    //constructor
    Mesh();
    Mesh(int ncells,int nGhost, double x_min, double x_max,double y_min,double y_max,double cfl, Cell (*speed_x)(double x, double y, double t), Cell (*speed_y)(double x, double y, double t), Cell (*level_set)(double x, double y,double Atime));

    ~Mesh();
    void Calculate_dt();
    void applyBC();
    void advect_level_set();
    void save_to_file(std::string name)const;


  };
}
#endif
