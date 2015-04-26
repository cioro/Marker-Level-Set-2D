#ifndef MLS_CELL_H
#define MLS_CELL_H

#include<vector>
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
namespace MLS{
  class Cell{

  public:
    double phi;
    double phi_u,phi_v;

    //constructor
    Cell();
    Cell(double arg_phi, double phi_u, double phi_v);
 
    ~Cell();
    void print()const;

  };
}
#endif
