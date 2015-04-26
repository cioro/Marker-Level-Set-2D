#include<vector>
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
#include"MLS_mesh.hpp"
#include"MLS_cell.hpp"
namespace MLS{
  //constructor
  Cell::Cell(){};
  Cell::Cell(double arg_phi, double arg_phi_u, double arg_phi_v):phi(arg_phi), phi_u(arg_phi_u),phi_v(arg_phi_v){};
 
  Cell::~Cell(){};
  void Cell::print()const{
    std::cout << phi << "\t" <<phi_u << "\t" << phi_v << "\n";
 
  };

}
