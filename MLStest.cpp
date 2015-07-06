#define _USE_MATH_DEFINES 

#include"MLS_mesh.hpp"
#include"MLS_cell.hpp"
#include<vector>
#include <unordered_map>
#include<functional>
#include<cmath>
#include<fstream>
#include<iostream>
#include<assert.h>
#include<string>
#include<blitz/array.h>

//Speed functions
double Z_speed_x(double x, double y, double t, double T){
  return ((M_PI/314)*(0.5-y));
}

double Z_speed_y(double x, double y, double t, double T){
  return ((M_PI/314)*(x-0.5));
}

double spiral_speed_x(double x, double y, double t,double T){
  return 0.0;//(-2*M_PI*sin(M_PI*x)*sin(2*M_PI*y)*cos(M_PI*t/T));
}
double spiral_speed_y(double x, double y, double t, double T){
  return -1.0;//(2*M_PI*sin(M_PI*y)*sin(2*M_PI*x)*cos(M_PI*t/T));
}

MLS::Cell Zalesak_disk(double x, double y,double time, double T_max){
  
  MLS::Cell cell;
  double phi_circle = 0.0;
  double r = (3.0/20.0);
  double x_c = 0.5;
  double y_c = 0.75;

  //Inside the circle-solid
  if( ((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c)) <= (r*r) ){
    
    phi_circle = r - sqrt((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c));
    
  }else if( ((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c)) > (r*r) ){

    phi_circle = r - sqrt((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c));
  
  }else{
    std::cout << "CIRCLE LEVEL_SET. error in setting phi" << "\n";
  }

  double phi_rectangle = 0.0;
  double square_diam_x = 0.025;
  double square_diam_y = 0.125;
  
  double x_0 = 0.5;
  double y_0 = 0.725;
  
  double x_dis;
  double y_dis;
  
  double x_prime = std::abs(x-x_0);
  double y_prime = std::abs(y-y_0);

  x_dis = x_prime - square_diam_x;
  y_dis = y_prime - square_diam_y;

  if(x_dis <= 0  && y_dis <= 0){ // case both inside
    if(x_dis > y_dis){
      phi_rectangle = -x_dis;
    }else{
      phi_rectangle = -y_dis;
    }
  }else if(x_dis <= 0 && y_dis > 0){ // case y outside, x inside
    phi_rectangle = -y_dis;
  }else if(y_dis <= 0 && x_dis > 0){ // case x outside, y inside
    phi_rectangle = -x_dis;
  }else{ // case both outside
    phi_rectangle = -std::sqrt(y_dis*y_dis + x_dis*x_dis);
  }
  
  phi_rectangle = -1*phi_rectangle;
  double phi;
  phi = std::min(phi_circle, phi_rectangle);



  cell.phi = phi;
  cell.phi_u = (M_PI/314)*(0.5-y);
  cell.phi_v = (M_PI/314)*(x-0.5);
  return cell;

}


MLS::Cell level_set_diagonal(double x, double y, double time, double T_max){

  MLS::Cell cell;
  double phi = 0.0;

  phi = x + y - 1;
  cell.phi = phi;
  cell.phi_u = 0;
  cell.phi_v = 0;
  
  return cell;

}


MLS::Cell level_set_wall(double x, double y,double time, double T_max){
  
  MLS::Cell cell;
  double phi = 0.0;
  const double x_0 = 0.5;
  
  //Solid
  if (y<=x_0){
    phi = (x_0-y);   
    
  }else if (y>x_0){
    phi = (x_0-y);
    
  }else{
    std::cout<<"Something went wrong inside the level set function in main.cpp"<<std::endl;
    
  }  
  cell.phi = phi;
  cell.phi_u = 0.0;
  cell.phi_v = 0.0;
  
  return cell;

}





MLS::Cell zero(double x, double y, double t, double T){
  MLS::Cell cell;
 
  cell.phi = 2;
  cell.phi_u = (M_PI/314)*(0.5-y);
  cell.phi_v = (M_PI/314)*(x-0.5);
   
  return cell;
}

//Inialising level set fcn
MLS::Cell level_set_circle(double x, double y,double t,double T){
  
  MLS::Cell cell;
  double phi = 0.0;
  double r = (3.0/20.0);
  double x_0 = 0.5;
  double y_0 = 0.75;

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

int main(int argc, char* argv[]){

  //Set parameters cfl, x_min,x_max
  int ncells=atof(argv[1]);
  int nGhost=5;
  double x_min=-0.5;
  double x_max=1.5;
  double y_min=-0.5;
  double y_max=1.5;
  double cfl=0.6;
  double T_max = 628.3185;
  int numMarkers = ncells*4;
  //Construct Level set mesh
  MLS::Mesh m(numMarkers,T_max,ncells, nGhost, x_min,x_max, y_min, y_max, cfl, Z_speed_x, Z_speed_y,Zalesak_disk);
  m.applyBC();
  
  m.marchingSquares();
  //m.correction1();
  /*
  std::string Snap = "Snap_";
  //m.save_to_file(Snap);
  m.vtk_output(Snap);
  m.vtk_output_marker(Snap);
  
  //Evolution loop
  
   
 for(double t = 0; t < T_max; t += m.dt){

    m.Calculate_dt();
    if(m.time+m.dt > T_max){
      m.dt = T_max-m.time;
    }
    //std::cout <<m.time << "\n";
    //m.advect_level_set();    
    // m.advect_WENO();
    // m.advect_RK_WENO_periodic();
    // m.advect_RK();
    m.advect_RK_WENO();
    m.applyBC();
    m.advect_markers();
    m.correction1();
    if((m.iter_counter%10) == 0){
      //m.save_to_file(Snap);
      m.vtk_output(Snap);
      m.vtk_output_marker(Snap);
    }
    
    m.time += m.dt;
  
    m.iter_counter++;
    }
 //std::cout<< m.time << "\n";
 // m.save_to_file(Snap);
 m.vtk_output(Snap);
 m.vtk_output_marker(Snap);
  */
 return 0;
  
}
