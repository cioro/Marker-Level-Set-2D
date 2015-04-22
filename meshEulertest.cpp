#include"Mesh2D.hpp"
#include"Euler2D.hpp"
#include"MLS_2D.hpp"
#include<vector>
#include<fstream>
#include<iostream>
#include<assert.h>
#include<string>
#include<sstream>
#include<blitz/array.h>

MLS level_set(double x, double y){
  double arg_phi = -1.0;
  MLS level_set;
  level_set.phi = arg_phi;
  return level_set;
}

MLS level_set_circle(double x, double y){
  
  MLS level_set;
  double arg_phi = 0.0;
  double r = 0.2;
  double x_0 = 0.6;
  double y_0 = 0.5;

  //Inside the circle-solid
  if( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) <= (r*r) ){
    
    arg_phi = r - sqrt((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0));
    
  }else if( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) > (r*r) ){

    arg_phi = r - sqrt((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0));
  
  }else{
    std::cout << "CIRCLE LEVEL_SET. error in setting phi" << "\n";
  }

  level_set.phi = arg_phi;
}

MLS level_set_square_update(double x, double y, double x_c, double y_c){

  MLS level_set;//level set default constructor
  double square_diam = 0.4;
  
  double x_0 = x_c;
  double y_0 = y_c;

  double x_dis;
  double y_dis;

  double x_prime = std::abs(x-x_0);
  double y_prime = std::abs(y-y_0);

  x_dis = x_prime - 0.5*square_diam;
  y_dis = y_prime - 0.5*square_diam;

  double arg_phi;

  if(x_dis <= 0  && y_dis <= 0){// case both inside
    if(x_dis > y_dis){
      arg_phi  = -x_dis;
    }else{
      arg_phi = -y_dis;
    }
  }else if(x_dis <= 0 && y_dis > 0){ // case y outside, x inside
    arg_phi = -y_dis;
  }else if(y_dis <= 0 && x_dis > 0){ // case x outside, y inside
    arg_phi = -x_dis;
  }else{// case both outside
    arg_phi = -std::sqrt(y_dis*y_dis + x_dis*x_dis);
  }
  level_set.phi = arg_phi;
  return level_set;
  
}
MLS level_set_square_moving(double x, double y){
  
  MLS level_set;
  double square_diam = 0.4;
  
  double x_0 = 1.5;
  double y_0 = 0.5;
  
  double x_dis;
  double y_dis;
  
  double x_prime = std::abs(x-x_0);
  double y_prime = std::abs(y-y_0);

  x_dis = x_prime - 0.5*square_diam;
  y_dis = y_prime - 0.5*square_diam;

  double arg_phi;

  if(x_dis <= 0  && y_dis <= 0){// case both inside
    if(x_dis > y_dis){
      arg_phi = -x_dis;
    }else{
      arg_phi = -y_dis;
    }
  }else if(x_dis <= 0 && y_dis > 0){ // case y outside, x inside
    arg_phi = -y_dis;
  }else if(y_dis <= 0 && x_dis > 0){ // case x outside, y inside
    arg_phi = -x_dis;
  }else{ // case both outside
    arg_phi = -std::sqrt(y_dis*y_dis + x_dis*x_dis);
  }
  level_set.phi= arg_phi;
  return level_set;
  
}


MLS level_set_square(double x, double y){

  MLS level_set;
  double square_diam = 0.4;
  
  double x_0 = 0.5;
  double y_0 = 0.5;
  
  double x_dis;
  double y_dis;
  
  double x_prime = std::abs(x-x_0);
  double y_prime = std::abs(y-y_0);

  x_dis = x_prime - 0.5*square_diam;
  y_dis = y_prime - 0.5*square_diam;

  double arg_phi;

  if(x_dis <= 0  && y_dis <= 0){ // case both inside
    if(x_dis > y_dis){
      arg_phi = -x_dis;
    }else{
      arg_phi = -y_dis;
    }
  }else if(x_dis <= 0 && y_dis > 0){ // case y outside, x inside
    arg_phi = -y_dis;
  }else if(y_dis <= 0 && x_dis > 0){ // case x outside, y inside
    arg_phi = -x_dis;
  }else{ // case both outside
    arg_phi = -std::sqrt(y_dis*y_dis + x_dis*x_dis);
  }
  level_set.phi=arg_phi;
  return level_set;
  
}

MLS level_set_diamond(double x, double y){

  MLS level_set;
  double arg_phi;
  double square_diam = 0.4;
  double x_0 = 0.5;
  double y_0 = 0.5;
  
  double x_dis = std::abs(x-x_0);
  double y_dis = std::abs(y-y_0);
  
  arg_phi= (0.5*square_diam - x_dis - y_dis) * std::sqrt(2) * 0.5;
  level_set.phi =arg_phi;

  return level_set;
}
MLS level_set_wall(double x, double y){
  double arg_phi = 0.0;
  const double x_0 = 0.5;
  
  //Solid
  if (x<=x_0){
    arg_phi = (x_0-x);   
    
  }else if (x>x_0){
    arg_phi = (x_0-x);
    
  }else{
    std::cout<<"Something went wrong inside the level set function in main.cpp"<<std::endl;
   
  }  

  MLS level_set;
  level_set.phi = arg_phi;

  return level_set;

}



//Reflective boundary function
Euler::W_state Refl_Left_Bound(Euler::W_state w){
     
  Euler::W_state w_result;
  
  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.v   =   w.v;
  w_result.P   =   w.P;

  return w_result;
}

Euler::W_state Refl_Right_Bound(Euler::W_state w ){
 
 Euler::W_state w_result;  

  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.v   =   w.v;
  w_result.P   =   w.P;

  return w_result;
  }

//Transmissive boundary function
Euler::W_state Trans_Left_Bound(Euler::W_state w){
 
 Euler::W_state w_result;
   
  w_result.rho =  w.rho;
  w_result.u   =  w.u;
  w_result.v   =  w.v;
  w_result.P   =  w.P;
 
  return w_result;
}


Euler::W_state Trans_Right_Bound(Euler::W_state w){
 
  
  Euler::W_state w_result;
  
  w_result.rho =   w.rho;
  w_result.u   =   w.u;
  w_result.v   =   w.v;
  w_result.P   =   w.P;
 
  return w_result;
}
//Initial condition function
//Right going wave

Euler::U_state Ghost_moving_flow(double x, double y){

  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  const Euler::W_state wL(1.0,1.0,0.0,1.0);
  Euler::U_state uL = e.CfromP(wL);
  return uL;

}

Euler::U_state Ghost_moving_square(double x, double y){

  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  const Euler::W_state wL(1.0,0.0,0.0,1.0);
  Euler::U_state uL = e.CfromP(wL);
  return uL;

}


Euler::U_state Ghost_square(double x, double y){
 
  const double x_0 = 0.2;
  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  
  if (x<=x_0){
    
    const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,0.0,0.0,1.0);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital ghost_square function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}

Euler::U_state Ghost_diamond(double x, double y){
 
  const double x_0 = 0.2;
  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  
  if (x<=x_0){
    
    const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,0.0,0.0,1.0);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital square_circle function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state Ghost_circle(double x, double y){
 
  const double x_0 = 0.2;
  
  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  
  if (x<=x_0){
    
    const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,0.0,0.0,1.0);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital Ghost_circle function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state Sod_xdir(double x, double y){

  const double x_0 = 0.5;
  //Should do this assert(x_0 < x_max && x_0 > x_min); to avoid errors
  //but I would have to add to more args to the fcn.

  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  
  if (x<=x_0){
    
    const Euler::W_state wL(0.125,0.0,0.0,0.1);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,0.75,0.0,1.0);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }
}

Euler::U_state Sod_ydir(double x,double y){
  
  const double y_0 = 0.3;
  
  Euler e;

 if (y<=y_0){
    
   const Euler::W_state wL(1.0,0.0,0.75,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (y>y_0){
    const Euler::W_state wR(0.125,0.0,0.0,0.1);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }
}

Euler::U_state cylinder(double x,double y){
  
  const double x_0 = 1.0;
  const double y_0 = 1.0;
  const double r  = 0.4;
 
  Euler e;

  if ( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) <= (r*r)){
    
    const Euler::W_state wL(1.0,0.0,0.0,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }else  if ( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) > (r*r)  ){
    const Euler::W_state wR(0.125,0.0,0.0,0.1);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside initial cylindrical explosion fcn"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state diagonal(double x,double y){
  
  const double x_0 = 0.4;
  
  Euler e;

  if (y < (1-x)){
    
    const Euler::W_state wL(1.0,0.0,0.0,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }else if (y >= (1-x)){
    const Euler::W_state wR(0.125,0.0,0.0,0.1);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}

Euler::U_state diagonal_downwards(double x,double y){
  
  const double x_0 = 0.4;
  
  Euler e;

  if (y < (1-x)){
    
    const Euler::W_state wL(0.125,0.0,0.0,0.1);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }else if (y >= (1-x)){

    const Euler::W_state wR(1.0,0.0,0.0,1.0);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state collision_perpendicular_to_x_axis(double x,double y){
  
  const double x_0 = 0.3;
  const double x_1 = 0.7;
  
  Euler e;

 if (x<=x_0){
    
   const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
   Euler::U_state uL = e.CfromP(wL);
   return uL;
 }else if(x_0 < x && x < x_1){
   const Euler::W_state wM(1.0,0.0,0.0,1.0);
   Euler::U_state uM = e.CfromP(wM);
   return uM;
 }else if (x >= x_1){
   const Euler::W_state wR(1.3764,-0.394,0.0,1.5698);
   Euler::U_state uR = e.CfromP(wR);
   return uR;
 }else{
   std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
   Euler::U_state u;
   return u;
 }

}

Euler::U_state collision_perpendicular_to_y_axis(double x,double y){
  
  const double y_0 = 0.3;
  const double y_1 = 0.7;
  
  Euler e;

 if (y<=y_0){
    
   const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
   Euler::U_state uL = e.CfromP(wL);
   return uL;
 }else if(y_0 < y && y < y_1){
   const Euler::W_state wM(1.0,0.0,0.0,1.0);
   Euler::U_state uM = e.CfromP(wM);
   return uM;
 }else if (y >= y_1){
   const Euler::W_state wR(1.3764,-0.394,0.0,1.5698);
   Euler::U_state uR = e.CfromP(wR);
   return uR;
 }else{
   std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
   Euler::U_state u;
   return u;
 }

}

double xTransform(double v_n, double v_t, double a)
{
	double  norm = 1/std::sqrt(1 + a*a);
		
	return v_t * norm - v_n * a * norm;
}

double yTransform(double v_n, double v_t, double a)
{
	double  norm = 1/std::sqrt(1 + a*a);

	return v_t * a * norm + v_n * norm;
}

Euler::U_state collision_diagonal(double x,double y){
  
  const double x_0 = 0.3;
  const double x_1 = 0.7;
 
  Euler e;
  

  if(y <= (0.3-x) ){
    std::cout << " y : " << y << " 0.3-x " << x << "\n";
    const Euler::W_state wL(1.3764,0.2786,0.2786,1.5698);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }else if(((0.3-x) < y) && (y < (1.7-x))){
    const Euler::W_state wM(1.0,0.0,0.0,1.0);
    Euler::U_state uM = e.CfromP(wM);
    return uM;
  }else if(y >= (1.7-x)){
    std::cout << " y : " << y << " 0.7-x " << x << "\n";
    const Euler::W_state wR(1.3764,-0.2786,-0.2786,1.5698);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
 }

}

Euler::U_state collision_angle(double x,double y){
  
  const double x_0 = 0.3;
  const double x_1 = 0.7;
  
  Euler e;

  if (y <= (0.3-1.5*x)  ){
    
    const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
   Euler::U_state uL = e.CfromP(wL);
   return uL;
  }else if((0.3-1.5*x) > y && y < (0.7-1.5*x)){
    const Euler::W_state wM(1.0,0.0,0.0,1.0);
   Euler::U_state uM = e.CfromP(wM);
   return uM;
  }else if (y >= (0.7-1.5*x)){
    const Euler::W_state wR(1.3764,-0.394,0.0,1.5698);
   Euler::U_state uR = e.CfromP(wR);
   return uR;
 }else{
   std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
   Euler::U_state u;
   return u;
 }

}

//Exact riemann solver function

Euler::U_state Exact_solver(double x){

  Euler::U_state u_empty;

  return u_empty;
}


int main(int argc, char* argv[]){


  //Parameters of the problem
  double x_min = 0.0, x_max = 1.0; //domain length
  double y_min = 0.0, y_max = 1.0;
  if(argv[1]==std::string("cylinder")){
      x_max = 2.0;
      y_max = 2.0;
    }

  //Object Level SET DATA
  double obj_vel_x = 0.0;
  double obj_vel_y = 0.0;
  double obj_x_c = 0.5;
  double obj_y_c = 0.5;
 

  if(argv[1]==std::string("ghost_moving_square")){
      x_max = 2.0;
      obj_vel_x = -1.0;
      obj_vel_y = 0.0;
      obj_x_c = 1.5;
      obj_y_c = 0.5;
  }

 if(argv[1]==std::string("ghost_moving_flow")){
      x_max = 2.0;
    }
 

  double cfl = 0.9;
  int ncells;
  if(argc > 2){
    ncells = atof(argv[2]);
  }else{
    ncells = 100;
  }
  double dt;
  int nGhost = 2;  
  double T_max;
  std::string BC;
  std::string Snap;
  std::string slice_x_axis;
  std::string slice_y_axis;
  std::string file_init_w;
  std::string file_output;
  Mesh m(ncells, x_min, x_max,y_min,y_max,cfl, Sod_xdir, Trans_Left_Bound, Trans_Right_Bound, nGhost,BC, level_set);

  if(argc == 4){
    if(argv[3]==std::string("transmissive")){
      BC = "Transmissive";
      m.reset_BC(Trans_Left_Bound, Trans_Right_Bound,BC);
    }else if(argv[3]==std::string("reflective")){
      BC = "Reflective";
      m.reset_BC(Refl_Left_Bound, Refl_Right_Bound,BC);
    }else{
      std::cout << "No BC specified. Default Reflective BC";
      BC = "Reflective";
      m.reset_BC(Refl_Left_Bound, Refl_Right_Bound,BC);
    }
 }
  m.applyBC();
  
  if(argv[1]==std::string("Sod_xdir")){
    T_max = 0.20;
    m.reset(Sod_xdir);
    Snap = "Sod_xdir_snap_";
    slice_x_axis = "slice_x_axis_Sod_xdir_";
    slice_y_axis = "slice_y_axis_Sod_xdir_";
    file_init_w = "Initial_Sod_xdir_";
    file_output = "Sod_xdir_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing Sod_xdir test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("Sod_ydir")){ 
    T_max = 0.25;
    m.reset(Sod_ydir);
    Snap = "Sod_ydir_snap_";
    slice_x_axis = "slice_x_axis_Sod_ydir_";
    slice_y_axis = "slice_y_axis_Sod_ydir_";
    file_init_w = "Initial_Sod_ydir_";
    file_output = "Sod_ydir_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing Sod_ydir test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("ghost_circle")){ 
    T_max = 0.4;
    m.reset(Ghost_circle,level_set_circle);
    Snap = "ghost_circle_snap_";
    slice_x_axis = "slice_x_axis_ghost_circle_";
    slice_y_axis = "slice_y_axis_ghost_circle_";
    file_init_w = "Initial_ghost_circle_";
    file_output = "ghost_circle_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing ghost_circle test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
    }else if(argv[1]== std::string("ghost_square")){ 
    T_max = 0.4;
    m.reset(Ghost_square,level_set_square);
    Snap = "ghost_square_snap_";
    slice_x_axis = "slice_x_axis_ghost_square_";
    slice_y_axis = "slice_y_axis_ghost_square_";
    file_init_w = "Initial_ghost_square_";
    file_output = "ghost_square_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing ghost_square test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
    }else if(argv[1]== std::string("ghost_diamond")){ 
    T_max = 0.4;
    m.reset(Ghost_diamond,level_set_diamond);
    Snap = "ghost_diamond_snap_";
    slice_x_axis = "slice_x_axis_ghost_diamond_";
    slice_y_axis = "slice_y_axis_ghost_diamond_";
    file_init_w = "Initial_ghost_diamond_";
    file_output = "ghost_diamond_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing ghost_diamond test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";

    }else if(argv[1]== std::string("ghost_moving_square")){ 
    T_max = 1.0;
    m.reset(Ghost_moving_square,level_set_square_moving);
    Snap = "ghost_moving_square_snap_";
    slice_x_axis = "slice_x_axis_ghost_moving_square_";
    slice_y_axis = "slice_y_axis_ghost_moving_square_";
    file_init_w = "Initial_ghost_moving_square_";
    file_output = "ghost_moving_square_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing ghost_moving_square test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("ghost_moving_flow")){ 
    T_max = 1.0;
    m.reset(Ghost_moving_flow,level_set_square);
    Snap = "ghost_moving_flow_snap_";
    slice_x_axis = "slice_x_axis_ghost_moving_flow_";
    slice_y_axis = "slice_y_axis_ghost_moving_flow_";
    file_init_w = "Initial_ghost_moving_flow_";
    file_output = "ghost_moving_flow_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing ghost_moving_flow test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("cylinder")){ 
    T_max = 0.25;
    m.reset(cylinder);
    Snap = "cylinder_snap_";
    slice_x_axis = "slice_x_axis_cylinder_";
    slice_y_axis = "slice_y_axis_cylinder_";
    file_init_w = "Initial_cylinder_";
    file_output = "cylinder_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing cylinder test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("diagonal")){ 
    T_max = 0.25;
    m.reset(diagonal);
    Snap = "diagonal_snap_";
    slice_x_axis = "slice_x_axis_diagonal_";
    slice_y_axis = "slice_y_axis_diagonal_";
    file_init_w = "Initial_diagonal_";
    file_output = "diagonal_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing diagonal test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("diagonal_downwards")){ 
    T_max = 0.25;
    m.reset(diagonal);
    Snap = "diagonal_downwards_snap_";
    slice_x_axis = "slice_x_axis_diagonal_downwards_";
    slice_y_axis = "slice_y_axis_diagonal_downwards_";
    file_init_w = "Initial_diagonal_downwards_";
    file_output = "diagonal_output_downwards_";
    
  }else if(argv[1]== std::string("collision_perpendicular_to_x_axis")){ 
    T_max = 0.4;
    m.reset(collision_perpendicular_to_x_axis);
    Snap = "collision_perpendicular_to_x_axis_snap_";
    slice_x_axis = "slice_x_axis_collision_perpendicular_to_x_axis_";
    slice_y_axis = "slice_y_axis_collision_perpendicular_to_x_axis_";
    file_init_w = "Initial_collision_perpendicular_to_x_axis_";
    file_output = "collision_perpendicular_to_x_axis_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing collsion perpedicular to x axis test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("collision_perpendicular_to_y_axis")){ 
    T_max = 0.4;
    m.reset(collision_perpendicular_to_y_axis);
    Snap = "collision_perpendicular_to_y_axis_snap_";
    slice_x_axis = "slice_x_axis_collision_perpendicular_to_y_axis_";
    slice_y_axis = "slice_y_axis_collision_perpendicular_to_y_axis_";
    file_init_w = "Initial_collision_perpendicular_to_y_axis_";
    file_output = "collision_perpendicular_to_y_axis_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing collision_perpendicular_to_y_axis test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("collision_diagonal")){ 
    T_max = 0.4;
    m.reset(collision_diagonal);
    Snap = "collision_diagonal_snap_";
    slice_x_axis = "slice_x_axis_collision_diagonal_";
    slice_y_axis = "slice_y_axis_collision_diagonal_";
    file_init_w = "Initial_collision_diagonal_";
    file_output = "collision_diagonal_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing collision_diagonal test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else if(argv[1]== std::string("collision_angle")){ 
    T_max = 0.4;
    m.reset(collision_angle);
    Snap = "collision_angle_snap_";
    slice_x_axis = "slice_x_axis_collision_angle_";
    slice_y_axis = "slice_y_axis_collision_angle_";
    file_init_w = "Initial_collision_angle_";
    file_output = "collision_angle_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Performing collision_angle test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }else{
    T_max = 0.25;
    m.reset(cylinder);
    Snap = "cylinder_snap_";
    slice_x_axis = "slice_x_axis_cylinder_";
    slice_y_axis = "slice_y_axis_cylinder_";
    file_init_w = "Initial_cylinder_";
    file_output = "cylinder_output_";
    m.save_w_state(file_init_w);
    std::cout <<"Default case :Performing cylinder test with ncells: "<<ncells <<" up to time T_max: " << T_max << "\n";
  }
 
  m.applyBC();
  m.applyGhost(obj_vel_x,obj_vel_y);
  
  for(double t = 0; t<T_max; t+=dt){
    std::cout <<"Inside the main, time step is : " << m.time << "\n";
   
    
    dt = m.Calculate_dt();
    std::cout << "Current time step size: " << dt << "\n";
    //   flux_and_update(m,dt,std::string("XY"));
    m.applyBC();
    m.applyGhost(obj_vel_x,obj_vel_y);
    if(argv[1]== std::string("ghost_moving_square")){
      //m.update_level_set(dt,m.time,obj_vel_x,obj_vel_y,obj_x_c,obj_y_c,level_set_square_update);
       m.advect_level_set(dt,m.time,obj_vel_x,obj_vel_y);
    }  
    m.time++;

    t+=dt;

    dt = m.Calculate_dt();
    //flux_and_update(m,dt,std::string("YX"));
    m.applyBC();
    m.applyGhost(obj_vel_x,obj_vel_y);
    if(argv[1]== std::string("ghost_moving_square")){
      // m.update_level_set(dt,m.time,obj_vel_x,obj_vel_y,obj_x_c,obj_y_c,level_set_square_update);
      m.advect_level_set(dt,m.time,obj_vel_x,obj_vel_y);
   }
     m.time++;
    
    m.save_w_state(Snap);
    m.slice_x_axis(slice_x_axis);
    m.slice_y_axis(slice_y_axis);

  }
  
  m.save_w_state(file_output);
  
}
