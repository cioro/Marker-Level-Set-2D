#ifndef MLS2D_H
#define MLS2D_H

#include<vector>
#include"Mesh2D.hpp"
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
#include"Euler2D.hpp"

class MLS{

public:
  double phi;

  //constructor
  MLS();
  MLS(double arg_phi);
 
};

#endif
