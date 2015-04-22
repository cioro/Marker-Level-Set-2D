#include<vector>
#include"Mesh2D.hpp"
#include"MLS_2D.hpp"
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
#include"Euler2D.hpp"

MLS::MLS() : phi (0) {}
MLS::MLS(double arg_phi): phi(arg_phi) {}


