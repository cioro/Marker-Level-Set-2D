#define _USE_MATH_DEFINES 

#include"MLS_mesh.hpp"
#include"MLS_cell.hpp"
#include<vector>
#include<cmath>
#include<fstream>
#include<iostream>
#include<assert.h>
#include<string>
#include<blitz/array.h>

//Speed functions

double spiral_speed_x(double x, double y, double t,double T){
return  (-2*M_PI*sin(M_PI*x)*sin(2*M_PI*y)*cos(M_PI*t/T));
}
double spiral_speed_y(double x, double y, double t, double T){
return  (2*M_PI*sin(M_PI*y)*sin(2*M_PI*x)*cos(M_PI*t/T));
}

//Inialising level set fcn
MLS::Cell level_set_circle(double x, double y,double t,double T){
  
MLS::Cell cell;
double phi = 0.0;
double r = 0.2;
double x_0 = 0.6;
double y_0 = 0.5;

//Inside the circle-solid
if( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) <= (r*r) ){
    
phi = r - sqrt((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0));
    
}else if( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) > (r*r) ){

phi = r - sqrt((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0));
  
}else{
std::cout << "CIRCLE LEVEL_SET. error in setting phi" << "\n";
}
cell.phi = phi;
cell.phi_u = (-2*M_PI*sin(M_PI*x)*sin(2*M_PI*y)*cos(M_PI*t/T));
cell.phi_v = (2*M_PI*sin(M_PI*y)*sin(2*M_PI*x)*cos(M_PI*t/T));

return cell;
}

int main(){

//Set parameters cfl, x_min,x_max
int ncells=1000;
int nGhost=1;
double x_min=0.0;
double x_max=2.0;
double y_min=0.0;
double y_max=2.0;
double cfl=0.9;
double T_max = 2.0;
//Construct Level set mesh
MLS::Mesh m(T_max,ncells, nGhost, x_min,x_max, y_min, y_max, cfl,spiral_speed_x, spiral_speed_y, level_set_circle);
m.applyBC();
m.Calculate_dt();

std::string Snap = "Snap_";

//applyBC

//Evolution loop
for(double t=0; t<T_max;t+=m.dt){

std::cout <<m.time << "\n";
m.advect_level_set();    
m.Calculate_dt();
m.applyBC();
m.iter_counter++;
m.save_to_file(Snap);
m.time += m.dt;
}

}
