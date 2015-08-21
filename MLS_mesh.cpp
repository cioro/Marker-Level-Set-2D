#include<vector>
#include"MLS_cell.hpp"
#include"MLS_mesh.hpp"
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<functional>
#include<algorithm>
#include <unordered_map>
#include <initializer_list>



namespace MLS{

  Mesh::Mesh(){};
  Mesh::Mesh(int AnumM, double AT_max,int arg_ncells, int AnGhost, double Ax_min, double Ax_max,double Ay_min,double Ay_max,double arg_cfl, double (*s_x)(double x, double y, double t, double T), double (*s_y)(double x, double y, double t, double T), Cell (*level_set)(double x, double y, double time, double T_max)) :T_max(AT_max), ncells(arg_ncells),nGhost(AnGhost), x_min(Ax_min),x_max(Ax_max),y_min(Ay_min),y_max(Ay_max),cfl(arg_cfl), iter_counter(0),time(0), speed_x(s_x), speed_y(s_y),numMarkers(AnumM)
  {
     
    dx = (x_max-x_min)/(double)ncells;
    dy = (y_max-y_min)/(double)ncells;
 
    MLS_data.resize(ncells + 2*nGhost,ncells + 2*nGhost);
    xaxis.resize(ncells + 2*nGhost +1);
    yaxis.resize(ncells + 2*nGhost +1);
    x_cell_axis.resize(ncells+2*nGhost);
    y_cell_axis.resize(ncells+2*nGhost);
   
    int x_counter = 0;//Used to correctly calculate value of xaxis (value starts at 0 not nGhost,
    int y_counter = 0;//but first element is nGhost; Difference between array index and represented value
   

    //fill in xaxis and y axis
    for(int i = nGhost; i <(ncells+nGhost+1); i++){
      xaxis(i) = x_min + x_counter*dx;
      x_counter++;
      for(int j = nGhost; j<(ncells+nGhost+1); j++){
	yaxis(j) = y_min + y_counter*dy;
	y_counter++;
      }
      y_counter=0;
     
    }
    
    x_counter=0;
    y_counter=0;
    
    for(int i = nGhost; i <(ncells+nGhost); i++){
      x_cell_axis(i) = dx*0.5 + x_min + x_counter*dx;
      x_counter++;
      for(int j = nGhost; j<(ncells+nGhost); j++){
	y_cell_axis(j) = dy*0.5 + y_min + y_counter*dy;
	y_counter++;
      }
      y_counter=0;
    }
    
    

    for(int i = nGhost; i <(ncells+nGhost); i++){
      for(int j = nGhost; j<(ncells+nGhost); j++){
	//Loop overage 
	//Loop average
	MLS_data(i,j) = level_set(x_cell_axis(i),y_cell_axis(j),time,T_max);
      }
    }

    

    //--------MARKER VECTOR INITALIZATION-------------
    //------------------------------------------------

    
    double r = (3.0/20.0);
    double x_c = 0.5;
    double y_c = 0.75;
    double alpha =(2*M_PI)/double(numMarkers);
    double beta = sin(0.025/r);
    //    std::cout << "The angle is " << beta << "\n";
    double beta_start = (1.5*M_PI-beta);
    double beta_end = (1.5*M_PI+beta);
    //std::cout << "The start angle is : " << beta_start << " the end angle is " << beta_end << "\n";
    //std::cout << "This is a tiny angle: " <<alpha <<"\n";
    double arg_x_coord = x_c;
    double arg_y_coord = y_c;
    double slot_length = sqrt((0.15)*(0.15)-(0.025)*(0.025));
    Particle particle;
    double slot_markers;
    double angle;
    for(int i = 0; i < numMarkers; i++){
      angle = i*alpha;
      //std::cout <<"This is the current angle: " << angle <<"\n";
      if(0 < angle && angle < beta_start){
	arg_x_coord =r*cos(angle);
	arg_y_coord =r*sin(angle);
	//std::cout << "\t x_coord =" << arg_x_coord << "\n";
	//std::cout << "\t y_coord =" << arg_y_coord << "\n";
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;//+y_c;
	MLS_markers.push_back(particle);
      }else if(beta_start <= angle && angle <= beta_end){
	//Leave blank for slot to be created
      }else if(beta_end < angle && angle <= 2*M_PI){
	arg_x_coord =r*cos(angle);
	arg_y_coord =r*sin(angle);
	//std::cout << "\t x_coord =" << arg_x_coord << "\n";
	//std::cout << "\t y_coord =" << arg_y_coord << "\n";
	//std::cout << "AFTER SLOT \n";
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;//+y_c;
	MLS_markers.push_back(particle);
      }
    }
    
    //std::cout << "Current number of markers : " << MLS_markers.size() << std::endl;
      double abs_slot_length = 0.25;
      int x_counter_slot = 1.5*(0.05)/(alpha*r);
      int y_counter_slot = 1.5*(0.25)/(alpha*r);
      double x_step = 0.05/double(x_counter_slot);
      double y_step = abs_slot_length/double(y_counter_slot);
            
      //double y_step = 0.01; 
      //left-side of slot(25)
      for(int i = 0; i <= y_counter_slot; ++i){
	arg_x_coord = -0.025;
	arg_y_coord = -slot_length + i*y_step;
     	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;//+y_c;
	MLS_markers.push_back(particle);

      }
      //top of slot(10)
      //double x_step = 0.005;
      for(int i = 0; i <= x_counter_slot; ++i){
	arg_y_coord = 0.25-slot_length;
	arg_x_coord = -0.025 + x_step*i;
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;//+y_c;
	MLS_markers.push_back(particle);
      }

      //righ-side of slot(25)
      for(int i = 0; i <= y_counter_slot; ++i){
	arg_x_coord = 0.025;
	arg_y_coord = -slot_length + i*y_step;
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;//+y_c;
	MLS_markers.push_back(particle);

      }
    
    std::cout << "The numer of markers is : " << MLS_markers.size() << "\n";
    //------------------------------------------------
    

    /*
    //--------MARKER VECTOR INITALIZATION-------------
    //------------------------------------------------

    
    double r = 1;
    double x_c = 0.49;
    double y_c = 0.75;
    //double c = 0.699;
    double a = 0.7;
    double b = 0.0875;//sqrt(a*a-c*c);  
    double t_alpha = atan2(a*tan(M_PI_2),b); 
    double t_step =(t_alpha)/double(0.25*numMarkers);
    double alpha = (2*M_PI)/double(numMarkers);
    //std::cout << "t_alpha " << t_alpha << " t_step " << t_step << "\n";
    //std::cout << "This is a tiny angle: " <<alpha <<"\n";
    double arg_x_coord = x_c;
    double arg_y_coord = y_c;
    Particle particle;
    double angle;
    double t;
    for(int i = 0; i < numMarkers; i++){
      angle = i*alpha;
      //std::cout <<"This is the current angle: " << angle <<"\n";
      if(0 <= angle && angle < M_PI_2){
	//t = atan2(a*tan(angle),b);
	t = i*t_step;
	arg_x_coord = a*r*cos(t);
	arg_y_coord = b*r*sin(t);
	//std::cout << "\t a  =" << a << "\n";
	//std::cout << "\t b  =" << b << "\n";
	//std::cout << "\t tan(angle)  =" << tan(angle) << "\n";
	//std::cout << "\t a*tan(angle)  =" << a*tan(angle) << "\n";
	//std::cout << "\t t  =" << t << "\n";
	//std::cout << "\t x_coord =" << arg_x_coord << "\n";
	//std::cout << "\t y_coord =" << arg_y_coord << "\n";
	//std::cout << "\n";
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;
	MLS_markers.push_back(particle);
      }else if(M_PI_2 <= angle && angle < M_PI){
	//t = atan2(a*tan(angle),b);
	t = i*t_step;
	arg_x_coord =std::copysign( a*r*cos(t),-1);
	arg_y_coord =std::copysign( b*r*sin(t),1);
	//std::cout << "\t a  =" << a << "\n";
	//std::cout << "\t b  =" << b << "\n";
	//std::cout << "\t tan(angle)  =" << tan(angle) << "\n";
	//std::cout << "\t a*tan(angle)  =" << a*tan(angle) << "\n";
	//std::cout << "\t t  =" << t << "\n";
	//std::cout << "\t x_coord =" << arg_x_coord << "\n";
	//std::cout << "\t y_coord =" << arg_y_coord << "\n";
	//std::cout << "\n";
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;
	MLS_markers.push_back(particle);
      }else if(M_PI <= angle && angle <1.5*M_PI){
	//t = atan2(a*tan(angle),b);
	t = i*t_step;
	arg_x_coord =std::copysign( a*r*cos(t),-1);
	arg_y_coord =std::copysign( b*r*sin(t),-1);
	//std::cout << "\t a  =" << a << "\n";
	//std::cout << "\t b  =" << b << "\n";
	//std::cout << "\t tan(angle)  =" << tan(angle) << "\n";
	//std::cout << "\t a*tan(angle)  =" << a*tan(angle) << "\n";
	//std::cout << "\t t  =" << t << "\n";
	//std::cout << "\t x_coord =" << arg_x_coord << "\n";
	//std::cout << "\t y_coord =" << arg_y_coord << "\n";
	//std::cout << "\n";
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;
	MLS_markers.push_back(particle);
      }else if(1.5*M_PI <= angle && angle < 2*M_PI){
	//t = atan2(a*tan(angle),b);
	t = i*t_step;
	arg_x_coord =std::copysign( a*r*cos(t),+1);
	arg_y_coord =std::copysign( b*r*sin(t),-1);
	//std::cout << "\t a  =" << a << "\n";
	//std::cout << "\t b  =" << b << "\n";
	//std::cout << "\t tan(angle)  =" << tan(angle) << "\n";
	//std::cout << "\t a*tan(angle)  =" << a*tan(angle) << "\n";
	//std::cout << "\t t  =" << t << "\n";
	//std::cout << "\t x_coord =" << arg_x_coord << "\n";
	//std::cout << "\t y_coord =" << arg_y_coord << "\n";
	//std::cout << "\n";
	particle.x_coord = arg_x_coord + x_c;
	particle.y_coord = arg_y_coord + y_c;
	MLS_markers.push_back(particle);
      }
    }
    
    
    std::cout << "The numer of markers is : " << MLS_markers.size() << "\n";
    //------------------------------------------------
    
    */

  };

  Mesh::~Mesh(){};

  void Mesh::Calculate_dt(){

    double speed_x =0.0;
    double speedtemp_x =0.0;
    double speed_y =0.0;
    double speedtemp_y =0.0;
    double min_coef = 0.0;//This is min of (S_x/dx,S_y/dy)-The minimum of the wave speed over the space step.
 
    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
     
	speedtemp_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	//std::cout << speedtemp_x << std::endl;
	speedtemp_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
      
	if(fabs(speedtemp_x) > fabs(speed_x)){
	  speed_x = speedtemp_x;
	}
	if(fabs(speedtemp_y) > fabs(speed_y)){
	  speed_y = speedtemp_y;
	}          
      }
    }
  
    min_coef = std::min(dx/fabs(speed_x),dy/fabs(speed_y));
    
    //If time < 5 then dt = 0.2
  
    if(iter_counter < 10){
      double  cfl_init = 0.2;
      dt = cfl_init*min_coef;
    
      //std::cout << "Inside calculate dt function, cfl = " << cfl_init << " dt " << dt << "\n"; 
    
    }else{
      //std::cout << "Inside calculate dt function, cfl = " << cfl << " dt " << dt << "\n";
  
      dt = cfl*min_coef;
    }

  
  };
  
  void Mesh::advect_level_set(){

    // std::cout << "ADVECTING" << "\n";

    double phi;
    double phi_xdir,phi_ydir;
    double phi_xdir_plus,phi_xdir_minus;
    double phi_ydir_plus,phi_ydir_minus;
    blitz::Array<double,2>phi_new;
    double obj_speed_x, obj_speed_y;
    
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    //The advectin is calculated using Strang-splitting
    //First Calculate XY
    // X-
    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){

	obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	phi = MLS_data(col,row).phi;
	phi_xdir_plus = MLS_data(col+1,row).phi;
	phi_xdir_minus = MLS_data(col-1,row).phi;

	if(obj_speed_x < 0){
	  phi_xdir = (phi_xdir_plus-phi)/dx; 
	}else if(obj_speed_x == 0){
	  phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  phi_xdir = (phi-phi_xdir_minus)/dx; 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
	
	phi = phi - dt*(obj_speed_x*phi_xdir);
	phi_new(col,row) = phi;
	
      }
    }
    applyBC(phi_new);
    //Y
    blitz::Array<double,2> phi_new_2;
    phi_new_2.resize(ncells+2*nGhost,ncells+2*nGhost);
    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){

	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	phi = phi_new(col,row);
	phi_ydir_plus = phi_new(col,row+1);
	phi_ydir_minus = phi_new(col,row-1);
    
	if(obj_speed_y < 0){
	  phi_ydir = (phi_ydir_plus-phi)/dy;
	}else if(obj_speed_y == 0){
	  phi_ydir = 0;
	}else if(obj_speed_y > 0){
	  phi_ydir = (phi-phi_ydir_minus)/dy;
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	phi = phi - dt*obj_speed_y*phi_ydir;
	phi_new_2(col,row) = phi;
	
      }
    }
    applyBC(phi_new_2);
    //recalculate time -step. In this case it is not too relevant since time step is constant as V only depends on x not t. 
    this->Calculate_dt();
    
    //Y
    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){

	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	phi = phi_new_2(col,row);
	phi_ydir_plus = phi_new_2(col,row+1);
	phi_ydir_minus = phi_new_2(col,row-1);
    
	if(obj_speed_y < 0){
	  phi_ydir = (phi_ydir_plus-phi)/dy;
	}else if(obj_speed_y == 0){
	  phi_ydir = 0;
	}else if(obj_speed_y > 0){
	  phi_ydir = (phi-phi_ydir_minus)/dy;
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	phi = phi - dt*obj_speed_y*phi_ydir;
	phi_new(col,row) = phi;
	
	
      }
    }
    applyBC(phi_new);
    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){

	obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	phi = phi_new(col,row);
	phi_xdir_plus = phi_new(col+1,row);
	phi_xdir_minus = phi_new(col-1,row);

	if(obj_speed_x < 0){
	  phi_xdir = (phi_xdir_plus-phi)/dx; 
	}else if(obj_speed_x == 0){
	  phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  phi_xdir = (phi-phi_xdir_minus)/dx; 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
	
	phi = phi - dt*(obj_speed_x*phi_xdir);
	phi_new_2(col,row) = phi;

      }
    }
    applyBC(phi_new_2);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	MLS_data(col,row).phi = phi_new_2(col,row);
      }
    }


};

  void Mesh::save_to_file(std::string filename)const{
    
     
    std::string dir = "data/";
    std::stringstream ss;
    ss << dir << filename << iter_counter <<"_"<< ncells;
    std::string tmppath = ss.str();
  
    //std::cout << "CREATING FILE \n";
    FILE * outfile = fopen(tmppath.c_str(),"w");
  
    for(int i = nGhost; i < ncells+nGhost; i++)
      {
	for(int j = nGhost; j < ncells+nGhost; j++){

	  fprintf(outfile, "%.8f \t %.8f \t %.8f \t %.8f \t %.8f  \n", xaxis(i), yaxis(j),MLS_data(i,j).phi,MLS_data(i,j).phi_u,MLS_data(i,j).phi_v);
	}

	fprintf(outfile,"\n");
      
      }
    fclose(outfile);
    
    std::string marker = "marker";
    std::stringstream ssm;
    ssm << dir << filename << marker << iter_counter;
    std::string tmppathssm = ssm.str();
  
    FILE * markerfile = fopen(tmppathssm.c_str(),"w");
    for(int i = 0; i < MLS_markers.size(); ++i){
      fprintf(markerfile, "%.8f \t %.8f \n", MLS_markers[i].x_coord,MLS_markers[i].y_coord);
    }

  };

  void Mesh::applyBC(){

    for(int x_dir = 0; x_dir < nGhost; x_dir++){
      //Loop over each row
      for(int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){
	MLS_data(x_dir,y_dir).phi= MLS_data(nGhost,y_dir).phi;
      }
    }

    //Right Boundary of Square Mesh
    //Loop over Ghost right columns
    for(int x_dir = 0; x_dir < nGhost; x_dir++ ){
      //Loop over each row
      for (int y_dir = nGhost; y_dir < ncells+nGhost; y_dir++){
	MLS_data(ncells+nGhost+x_dir,y_dir).phi= MLS_data(ncells+nGhost-1,y_dir).phi;
      }
    }

    //Top Boundary of Square Mesh
    //Loop over Ghost top rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++){
      //Loop over each column
      for(int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	MLS_data(x_dir,y_dir+ncells+nGhost).phi = MLS_data(x_dir,ncells+nGhost-1).phi;
      }
    }

    //Bottom Boundary of Square Mesh
    //Loop over Ghost bottom rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++ ){
      //Loop over each column
      for (int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	MLS_data(x_dir,y_dir).phi = MLS_data(x_dir,nGhost).phi;
      }
    }

  }

void Mesh::applyBC_periodic(){

    for(int x_dir = 0; x_dir < nGhost; x_dir++){
      //Loop over each row
      for(int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){
	MLS_data(x_dir,y_dir).phi= MLS_data(ncells+x_dir,y_dir).phi;
      }
    }

    //Right Boundary of Square Mesh
    //Loop over Ghost right columns
    for(int x_dir = 0; x_dir < nGhost; x_dir++ ){
      //Loop over each row
      for (int y_dir = nGhost; y_dir < ncells+nGhost; y_dir++){
	MLS_data(nGhost+ncells+x_dir,y_dir).phi= MLS_data(nGhost+x_dir,y_dir).phi;
      }
    }

    //Top Boundary of Square Mesh
    //Loop over Ghost top rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++){
      //Loop over each column
      for(int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	MLS_data(x_dir,nGhost+ncells+y_dir).phi = MLS_data(x_dir,nGhost+y_dir).phi;
      }
    }

    //Bottom Boundary of Square Mesh
    //Loop over Ghost bottom rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++ ){
      //Loop over each column
      for (int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	MLS_data(x_dir,y_dir).phi = MLS_data(x_dir,ncells+y_dir).phi;
      }
    }

  }


  double Mesh::D_minus(int i, int j, std::string dir,const blitz::Array<double,2> & input ){

    double phi=0;
    double phi_minus=0;
    double phi_result=0;
    if(dir == "x_dir"){
      phi = input(i,j);
      phi_minus = input((i-1),j);
      phi_result = (phi - phi_minus )/dx;
    }else if (dir == "y_dir"){
      phi = input(i,j);
      phi_minus = input(i,(j-1));
      phi_result = (phi - phi_minus )/dy;
    }
    return phi_result;

  }

  double Mesh::D_plus(int i, int j, std::string dir, const blitz::Array<double,2> & input){

    double phi=0;
    double phi_plus=0;
    double phi_result=0;
    if(dir == "x_dir"){
      phi = input(i,j);
      phi_plus = input((i+1),j);
      phi_result = (phi_plus - phi)/dx;
    }else if (dir == "y_dir"){
      phi = input(i,j);
      phi_plus = input(i,(j+1));
      phi_result = (phi_plus - phi)/dy;
    }
    return phi_result;

  }

  double Mesh::WENO(int i, int j, std::string stencil, std::string dir, const blitz::Array<double,2> & input){
    
    double S1,S2,S3;
    double a1,a2,a3;
    double w1,w2,w3;
    double v1,v2,v3,v4,v5;
    double phi1,phi2,phi3;
    double phi_result;
    double eps; //epsilon
   
    if(stencil == "minus" && dir == "x_dir"){ 
      //std::cout << "x_dir minus" << "\n"; 
      //std::cout << "The phi values are " << input((i-2),j) << "\t" << input((i-1),j) << "\t" << input((i),j) << "\t" \
	// << input((i+1),j) << "\t" << input((i+2),j) <<"\n";
      v1 = D_minus((i-2),j,dir,input);
      v2 = D_minus((i-1),j,dir,input);
      v3 = D_minus((i),j,dir,input);
      v4 = D_minus((i+1),j,dir,input);
      v5 = D_minus((i+2),j,dir,input);
    }else if(stencil == "plus" && dir == "x_dir"){
      v1 = D_plus((i+2),j,dir,input);
      v2 = D_plus((i+1),j,dir,input);
      v3 = D_plus((i),j,dir,input);
      v4 = D_plus((i-1),j,dir,input);
      v5 = D_plus((i-2),j,dir,input);
    }else if(stencil == "minus" && dir == "y_dir"){ 
      //std::cout << "y_dir minus" << "\n";
      //std::cout << "The phi values are " << input(i,(j-2)) << "\t" << input(i,(j-1)) << "\t" << input(i,j) << "\t" \
	//<< input(i,(j+1)) << "\t" << input(i,(j+2)) <<"\n";
      v1 = D_minus(i,(j-2),dir,input);
      v2 = D_minus(i,(j-1),dir,input);
      v3 = D_minus(i,j,dir,input);
      v4 = D_minus(i,(j+1),dir,input);
      v5 = D_minus(i,(j+2),dir,input);
    }else if(stencil == "plus" && dir == "y_dir"){
      v1 = D_plus(i,(j+2),dir,input);
      v2 = D_plus(i,(j+1),dir,input);
      v3 = D_plus(i,j,dir,input);
      v4 = D_plus(i,(j-1),dir,input);
      v5 = D_plus(i,(j-2),dir,input);
    }else{
      std::cout << "No upwinding differencing stencil selected" << "\n";
    }
   
    // std::cout << "The v values are " << v1 << "\t" << v2 << "\t"<< v3 << "\t"<< v4 << "\t"<< v5 << "\n";

    phi1 =  (1.0/3.0)*v1 - (7.0/6.0)*v2 + (11.0/6.0)*v3;
    phi2 = (-1.0/6.0)*v2 + (5.0/6.0)*v3 + (1.0/3.0)*v4;
    phi3 =  (1.0/3.0)*v3 + (5.0/6.0)*v4 - (1.0/6.0)*v5;
 
    S1 = (13.0/12.0)*(v1-2*v2+v3)*(v1-2*v2+v3) + \
      (1.0/4.0)*(v1-4*v2+3*v3)*(v1-4*v2+3*v3);

    S2 = (13.0/12.0)*(v2-2*v3+v4)*(v2-2*v3+v4) + \
      (1.0/4.0)*(v2-v4)*(v2-v4);

    S3 = (13.0/12.0)*(v3-2*v4+v5)*(v3-2*v4+v5) + \
      (1.0/4.0)*(3*v3-4*v4+v5)*(3*v3-4*v4+v5);
      
    std::vector<double> v = {(v1*v1),(v2*v2),(v3*v3),(v4*v4),(v5*v5)};
    eps = 1e-6*(*std::max_element(v.begin(),v.end()))+1e-99;

    a1 = 0.1/((S1+eps)*(S1+eps));
    a2 = 0.6/((S2+eps)*(S2+eps));
    a3 = 0.3/((S3+eps)*(S3+eps));

    double sum_a = a1+a2+a3;
      
    w1 = a1/sum_a;
    w2 = a2/sum_a;
    w3 = a3/sum_a;

    phi_result = w1*phi1 + w2*phi2 + w3*phi3;

    //std::cout << "The derivative is : " << phi_result << "\n";
    return phi_result;
 
  }



  blitz::Array<double,2> Mesh::spatial_WENO(blitz::Array<double,2> input ){

    // std::cout << "ADVECTING" << "\n";
    applyBC(input);
    double phi;
    double phi_xdir,phi_ydir;
    blitz::Array<double,2>phi_new;
    double obj_speed_x, obj_speed_y;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	
	phi = 0;

	if(obj_speed_x < 0){
	  stencil = "plus";
	  dir = "x_dir";
	  phi_xdir = WENO(col,row,stencil,dir,input); 
	}else if(obj_speed_x == 0){
	   phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  stencil = "minus";
	  dir = "x_dir";
	  phi_xdir = WENO(col,row,stencil,dir,input); 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
      
	if(obj_speed_y < 0){
	  stencil = "plus";
	  dir = "y_dir";
	  phi_ydir = WENO(col,row,stencil,dir,input);
	}else if(obj_speed_y == 0){
	  phi_ydir = 0.0;
	}else if(obj_speed_y > 0){
	  stencil = "minus";
	  dir = "y_dir";
	  phi_ydir = WENO(col,row,stencil,dir,input);
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	
	phi = obj_speed_x*phi_xdir+obj_speed_y*phi_ydir;
	phi_new(col,row) = phi;
      }
    }

    return phi_new;

  }

 blitz::Array<double,2> Mesh::spatial_WENO_X(blitz::Array<double,2> input ){

    // std::cout << "ADVECTING" << "\n";
    applyBC(input);
    double phi;
    double phi_xdir;
    blitz::Array<double,2>phi_new;
    double obj_speed_x;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	
	phi = 0;

	if(obj_speed_x < 0){
	  stencil = "plus";
	  dir = "x_dir";
	  phi_xdir = WENO(col,row,stencil,dir,input); 
	}else if(obj_speed_x == 0){
	   phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  stencil = "minus";
	  dir = "x_dir";
	  phi_xdir = WENO(col,row,stencil,dir,input); 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
	
	phi = obj_speed_x*phi_xdir;
	phi_new(col,row) = phi;
      }
    }
    applyBC(phi_new);
    return phi_new;

  }

 blitz::Array<double,2> Mesh::spatial_WENO_X_periodic(blitz::Array<double,2> input ){

    // std::cout << "ADVECTING" << "\n";
    applyBC_periodic(input);
    double phi;
    double phi_xdir;
    blitz::Array<double,2>phi_new;
    double obj_speed_x;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	
	phi = 0;

	if(obj_speed_x < 0){
	  stencil = "plus";
	  dir = "x_dir";
	  phi_xdir = WENO(col,row,stencil,dir,input); 
	}else if(obj_speed_x == 0){
	   phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  stencil = "minus";
	  dir = "x_dir";
	  phi_xdir = WENO(col,row,stencil,dir,input); 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
	
	phi = obj_speed_x*phi_xdir;
	phi_new(col,row) = phi;
      }
    }
    applyBC_periodic(phi_new);
    return phi_new;

  }


  blitz::Array<double,2> Mesh::spatial_WENO_Y(blitz::Array<double,2> input ){

    // std::cout << "ADVECTING" << "\n";
    applyBC(input);
    double phi;
    double phi_ydir;
    blitz::Array<double,2>phi_new;
    double obj_speed_y;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	
	phi = 0;
      
	if(obj_speed_y < 0){
	  stencil = "plus";
	  dir = "y_dir";
	  phi_ydir = WENO(col,row,stencil,dir,input);
	}else if(obj_speed_y == 0){
	  phi_ydir = 0.0;
	}else if(obj_speed_y > 0){
	  stencil = "minus";
	  dir = "y_dir";
	  phi_ydir = WENO(col,row,stencil,dir,input);
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	
	phi = obj_speed_y*phi_ydir;
	phi_new(col,row) = phi;
      }
    }

    applyBC(phi_new);
    return phi_new;

  }
      
   blitz::Array<double,2> Mesh::spatial_WENO_Y_periodic(blitz::Array<double,2> input ){

    // std::cout << "ADVECTING" << "\n";
    applyBC_periodic(input);
    double phi;
    double phi_ydir;
    blitz::Array<double,2>phi_new;
    double obj_speed_y;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	
	phi = 0;
      
	if(obj_speed_y < 0){
	  stencil = "plus";
	  dir = "y_dir";
	  phi_ydir = WENO(col,row,stencil,dir,input);
	}else if(obj_speed_y == 0){
	  phi_ydir = 0.0;
	}else if(obj_speed_y > 0){
	  stencil = "minus";
	  dir = "y_dir";
	  phi_ydir = WENO(col,row,stencil,dir,input);
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	
	phi = obj_speed_y*phi_ydir;
	phi_new(col,row) = phi;
      }
    }

    applyBC_periodic(phi_new);
    return phi_new;

  }
      
     


  blitz::Array<double,2> Mesh::spatial_first(blitz::Array<double,2> input){

    //std::cout << "ADVECTING" << "\n";
    applyBC(input);
    double phi;
    double phi_xdir,phi_ydir;
    double phi_xdir_plus,phi_xdir_minus;
    double phi_ydir_plus,phi_ydir_minus;
    blitz::Array<double,2>phi_new;
    double obj_speed_x, obj_speed_y;
    
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	phi = input(col,row);
	phi_xdir_plus = input(col+1,row);
	phi_xdir_minus = input(col-1,row);
	phi_ydir_plus = input(col,row+1);
	phi_ydir_minus = input(col,row-1);

	if(obj_speed_x < 0){
	  phi_xdir = (phi_xdir_plus-phi)/dx; 
	}else if(obj_speed_x == 0){
	  phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  phi_xdir = (phi-phi_xdir_minus)/dx; 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
      
	if(obj_speed_y < 0){
	  phi_ydir = (phi_ydir_plus-phi)/dy;
	}else if(obj_speed_y == 0){
	  phi_ydir = 0;
	}else if(obj_speed_y > 0){
	  phi_ydir = (phi-phi_ydir_minus)/dy;
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	// phi_ydir = 0;
	phi = obj_speed_x*phi_xdir+obj_speed_y*phi_ydir;
	phi_new(col,row) = phi;
      }
    }
    applyBC(phi_new);
    return phi_new; 

  };

 blitz::Array<double,2> Mesh::spatial_first_X(blitz::Array<double,2> input){
   
   applyBC(input);
   double phi;
   double phi_xdir;
   double phi_xdir_plus,phi_xdir_minus;
   blitz::Array<double,2>phi_new;
   double obj_speed_x;
   
   phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);
   
   for(int row = nGhost; row < nGhost+ncells; row++){
     for(int col = nGhost; col < nGhost+ncells; col++){
       obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);

       phi = input(col,row);
       phi_xdir_plus = input(col+1,row);
       phi_xdir_minus = input(col-1,row);
       
       if(obj_speed_x < 0){
	 phi_xdir = (phi_xdir_plus-phi)/dx; 
       }else if(obj_speed_x == 0){
	  phi_xdir = 0.0;
       }else if(obj_speed_x > 0){
	 phi_xdir = (phi-phi_xdir_minus)/dx; 
       }else{
	 std::cout << "Advection of level_set. Something went wrong" <<"\n";
       }
       
       phi = obj_speed_x*phi_xdir;
       phi_new(col,row) = phi;
     }
    }
   applyBC(phi_new);
   return phi_new; 
   
  };

  blitz::Array<double,2> Mesh::spatial_first_Y(blitz::Array<double,2> input){

    //std::cout << "ADVECTING" << "\n";
    applyBC(input);
    double phi;
    double phi_ydir;
    double phi_ydir_plus,phi_ydir_minus;
    blitz::Array<double,2>phi_new;
    double obj_speed_y;
    
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){

	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	phi = input(col,row);
	phi_ydir_plus = input(col,row+1);
	phi_ydir_minus = input(col,row-1);

      
	if(obj_speed_y < 0){
	  phi_ydir = (phi_ydir_plus-phi)/dy;
	}else if(obj_speed_y == 0){
	  phi_ydir = 0;
	}else if(obj_speed_y > 0){
	  phi_ydir = (phi-phi_ydir_minus)/dy;
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	// phi_ydir = 0;
	phi = obj_speed_y*phi_ydir;
	phi_new(col,row) = phi;
      }
    }
    applyBC(phi_new);
    return phi_new; 

  };

  void Mesh::advect_RK(){

    blitz::Array<double, 2> phi_n; 
    blitz::Array<double, 2> phi_n_one; 
    blitz::Array<double, 2> phi_n_two; 
    blitz::Array<double, 2> phi_n_half; 
    blitz::Array<double, 2> phi_n_three_half;

    blitz::Array<double, 2> V_phi_n; 
    blitz::Array<double, 2> V_phi_n_one; 
    blitz::Array<double, 2> V_phi_n_half; 

    phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_two.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_three_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    
    V_phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	phi_n(row,col)= MLS_data(row,col).phi;
      }
    }

      //-----------------------------------------------------------------

      V_phi_n = spatial_first(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
    
      V_phi_n_one = spatial_first(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }

      //----------------------------------------------------------------------------
      V_phi_n_half = spatial_first(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }
  
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  MLS_data(row,col).phi= phi_n_one(row,col);
	}
      }
      
  
  }


  void Mesh::applyBC(blitz::Array<double,2> & input){

     for(int x_dir = 0; x_dir < nGhost; x_dir++){
      //Loop over each row
      for(int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){
	input(x_dir,y_dir)= input(nGhost,y_dir);
      }
    }
    //Right Boundary of Square Mesh
    //Loop over Ghost right columns
    for(int x_dir = 0; x_dir < nGhost; x_dir++ ){
      //Loop over each row
      for (int y_dir = nGhost; y_dir < ncells+nGhost; y_dir++){
	input(ncells+nGhost+x_dir,y_dir)= input(ncells+nGhost-1,y_dir);
      }
    }

    //Top Boundary of Square Mesh
    //Loop over Ghost top rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++){
      //Loop over each column
      for(int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	input(x_dir,y_dir+ncells+nGhost) = input(x_dir,ncells+nGhost-1);
      }
    }

    //Bottom Boundary of Square Mesh
    //Loop over Ghost bottom rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++ ){
      //Loop over each column
      for (int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	input(x_dir,y_dir) = input(x_dir,nGhost);
      }
    }

    
  }

 void Mesh::applyBC_periodic(blitz::Array<double,2> & input){

     for(int x_dir = 0; x_dir < nGhost; x_dir++){
      //Loop over each row
      for(int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){
	input(x_dir,y_dir)= input(ncells+x_dir,y_dir);
      }
    }
    //Right Boundary of Square Mesh
    //Loop over Ghost right columns
    for(int x_dir = 0; x_dir < nGhost; x_dir++ ){
      //Loop over each row
      for (int y_dir = nGhost; y_dir < ncells+nGhost; y_dir++){
	input(ncells+nGhost+x_dir,y_dir)= input(nGhost+x_dir,y_dir);
      }
    }

    //Top Boundary of Square Mesh
    //Loop over Ghost top rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++){
      //Loop over each column
      for(int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	input(x_dir,y_dir+ncells+nGhost) = input(x_dir,nGhost+y_dir);
      }
    }

    //Bottom Boundary of Square Mesh
    //Loop over Ghost bottom rows
    for(int y_dir = 0; y_dir < nGhost; y_dir++ ){
      //Loop over each column
      for (int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	input(x_dir,y_dir) = input(x_dir,ncells+y_dir);
      }
    }

    
  }


  void Mesh::advect_RK_WENO_split(){

    //------------------------------------------------------------------


    blitz::Array<double, 2> phi_n; 
    blitz::Array<double, 2> phi_n_one; 
    blitz::Array<double, 2> phi_n_two; 
    blitz::Array<double, 2> phi_n_half; 
    blitz::Array<double, 2> phi_n_three_half;

    blitz::Array<double, 2> V_phi_n; 
    blitz::Array<double, 2> V_phi_n_one; 
    blitz::Array<double, 2> V_phi_n_half; 

    phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_two.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_three_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    
    V_phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	phi_n(row,col)= MLS_data(row,col).phi;
      }
    }

    if(iter_counter%2 == 0){
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------

      blitz::Array<double,2> phi_n_one_star;
      phi_n_one_star.resize(2*nGhost+ncells,2*nGhost+ncells);

      V_phi_n = spatial_WENO_X(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_WENO_Y(phi_n_one_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n_one_star(row,col) - dt*V_phi_n(row,col);
	}
      }

      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
      blitz::Array<double,2> phi_n_two_star;
      phi_n_two_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_one = spatial_WENO_X(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_WENO_Y(phi_n_two_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_two_star(row,col) - dt*V_phi_n_one(row,col);
	}
      }
      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }
      
      //Calculate phi^{n+3/2}
      //----------------------------------------------------------------------------
      blitz::Array<double,2> phi_n_three_half_star;
      phi_n_three_half_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_half = spatial_WENO_X(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_WENO_Y(phi_n_three_half_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_three_half_star(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }

    }else if(iter_counter%2 == 1){
       
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------

      blitz::Array<double,2> phi_n_one_star;
      phi_n_one_star.resize(2*nGhost+ncells,2*nGhost+ncells);

      V_phi_n = spatial_WENO_Y(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_WENO_X(phi_n_one_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n_one_star(row,col) - dt*V_phi_n(row,col);
	}
      }

      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
      blitz::Array<double,2> phi_n_two_star;
      phi_n_two_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_one = spatial_WENO_Y(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_WENO_X(phi_n_two_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_two_star(row,col) - dt*V_phi_n_one(row,col);
	}
      }
      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }

      //----------------------------------------------------------------------------
      blitz::Array<double,2> phi_n_three_half_star;
      phi_n_three_half_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_half = spatial_WENO_Y(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_WENO_X(phi_n_three_half_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_three_half_star(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }
    }

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	MLS_data(row,col).phi= phi_n_one(row,col);
      }
    }


    //-----------------------------------------------------------------
   
  }

void Mesh::advect_RK_WENO(){

    //------------------------------------------------------------------


    blitz::Array<double, 2> phi_n; 
    blitz::Array<double, 2> phi_n_one; 
    blitz::Array<double, 2> phi_n_two; 
    blitz::Array<double, 2> phi_n_half; 
    blitz::Array<double, 2> phi_n_three_half;

    blitz::Array<double, 2> V_phi_n; 
    blitz::Array<double, 2> V_phi_n_one; 
    blitz::Array<double, 2> V_phi_n_half; 

    phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_two.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_three_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    
    V_phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	phi_n(row,col)= MLS_data(row,col).phi;
      }
    }
       
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------
           
      V_phi_n = spatial_WENO(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }
      
      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
                
      V_phi_n_one = spatial_WENO(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }
     
      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }

      //----------------------------------------------------------------------------
         
      V_phi_n_half = spatial_WENO(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }
    

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	MLS_data(row,col).phi= phi_n_one(row,col);
      }
    }


    //-----------------------------------------------------------------
   
  }


void Mesh::advect_RK_WENO_periodic(){

    //------------------------------------------------------------------


    blitz::Array<double, 2> phi_n; 
    blitz::Array<double, 2> phi_n_one; 
    blitz::Array<double, 2> phi_n_two; 
    blitz::Array<double, 2> phi_n_half; 
    blitz::Array<double, 2> phi_n_three_half;

    blitz::Array<double, 2> V_phi_n; 
    blitz::Array<double, 2> V_phi_n_one; 
    blitz::Array<double, 2> V_phi_n_half; 

    phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_two.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_three_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    
    V_phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	phi_n(row,col)= MLS_data(row,col).phi;
      }
    }

    if(iter_counter%2 == 0){
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------

      blitz::Array<double,2> phi_n_one_star;
      phi_n_one_star.resize(2*nGhost+ncells,2*nGhost+ncells);

      V_phi_n = spatial_WENO_X_periodic(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_WENO_Y_periodic(phi_n_one_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n_one_star(row,col) - dt*V_phi_n(row,col);
	}
      }

      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
      blitz::Array<double,2> phi_n_two_star;
      phi_n_two_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_one = spatial_WENO_X_periodic(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_WENO_Y_periodic(phi_n_two_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_two_star(row,col) - dt*V_phi_n_one(row,col);
	}
      }
      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }
      
      //Calculate phi^{n+3/2}
      //----------------------------------------------------------------------------
      blitz::Array<double,2> phi_n_three_half_star;
      phi_n_three_half_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_half = spatial_WENO_X_periodic(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_WENO_Y_periodic(phi_n_three_half_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_three_half_star(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }

    }else if(iter_counter%2 == 1){
       
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------

      blitz::Array<double,2> phi_n_one_star;
      phi_n_one_star.resize(2*nGhost+ncells,2*nGhost+ncells);

      V_phi_n = spatial_WENO_Y_periodic(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_WENO_X_periodic(phi_n_one_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n_one_star(row,col) - dt*V_phi_n(row,col);
	}
      }

      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
      blitz::Array<double,2> phi_n_two_star;
      phi_n_two_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_one = spatial_WENO_Y_periodic(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_WENO_X_periodic(phi_n_two_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_two_star(row,col) - dt*V_phi_n_one(row,col);
	}
      }
      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }

      //----------------------------------------------------------------------------
      blitz::Array<double,2> phi_n_three_half_star;
      phi_n_three_half_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_half = spatial_WENO_Y_periodic(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_WENO_X_periodic(phi_n_three_half_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_three_half_star(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }
    }

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	MLS_data(row,col).phi= phi_n_one(row,col);
      }
    }


    //-----------------------------------------------------------------
  
}

  void Mesh::RK(Particle & p){

    double k1x = this-> speed_x(p.x_coord,p.y_coord,time,T_max);
    double k1y = this-> speed_y(p.x_coord,p.y_coord,time,T_max);
    
    double k2x = this-> speed_x(p.x_coord + 0.5*dt*k1x, p.y_coord + 0.5*dt*k1y,time + 0.5*dt,T_max);
    double k2y = this-> speed_y(p.x_coord + 0.5*dt*k1x, p.y_coord + 0.5*dt*k1y,time + 0.5*dt,T_max);
    
    double k3x = this-> speed_x(p.x_coord + 0.5*dt*k2x, p.y_coord + 0.5*dt*k2y,time + 0.5*dt,T_max);
    double k3y = this-> speed_y(p.x_coord + 0.5*dt*k2x, p.y_coord + 0.5*dt*k2y,time + 0.5*dt,T_max);
    
    double k4x = this-> speed_x(p.x_coord + dt*k3x, p.y_coord + dt*k3y,time + dt,T_max);
    double k4y = this-> speed_y(p.x_coord + dt*k3x, p.y_coord + dt*k3y,time + dt,T_max);

    double x = p.x_coord + (dt/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
    double y = p.y_coord + (dt/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
   
    p.x_coord = x;
    p.y_coord = y;
  
  }

  void Mesh::advect_markers(){

    for(int i = 0; i < MLS_markers.size(); i++){
      RK(MLS_markers[i]);
    }

  }

  void Mesh::vtk_output(std::string filename)const{
    
    std::string dir = "data/";
    std::string vtk = ".vtk";
    std::stringstream ss;
    ss << dir << filename << iter_counter << vtk ;
    std::string tmppath = ss.str();

    std::ofstream outFile;
    outFile.open(tmppath,std::ofstream::out);
    
    outFile.precision(8);
    outFile << "# vtk DataFile Version 2.0" << std::endl;
    outFile << "vtk output" << std::endl;
    outFile << "ASCII" << std::endl;
    outFile << "DATASET STRUCTURED_GRID" << std::endl;
    outFile << "DIMENSIONS " << (ncells+1) << " " << (ncells+1) << " " << 1.0 << std::endl;
    outFile << "POINTS" << " " << (ncells+1)*(ncells+1) << " " << "double" << std::endl;
    for(int i = nGhost; i < ncells+nGhost+1; i++){
	for(int j = nGhost; j < ncells+nGhost+1; j++){
	  outFile << xaxis(j) << "\t" <<  yaxis(i) << "\t" << 0.0  << std::endl;
	}
      }
    
    outFile << "POINT_DATA " << (ncells+1)*(ncells+1) << std::endl;
    outFile << "SCALARS Phi_Nodes float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(int j = nGhost; j < ncells+nGhost+1; j++){
      for(int i = nGhost; i < ncells+nGhost+1; i++){
	//Shouldn't there be an array access issue-not corner values are emty.
	outFile << 0.25*(MLS_data(i-1,j-1).phi+MLS_data(i-1,j).phi+MLS_data(i,j).phi+MLS_data(i,j-1).phi)  << std::endl;
	}
    }
    
    outFile << "CELL_DATA " << ncells*ncells << std::endl;
    outFile << "SCALARS Phi_Cells float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(int i = nGhost; i < ncells+nGhost; i++){
	for(int j = nGhost; j < ncells+nGhost; j++){
	  outFile << MLS_data(j,i).phi  << std::endl;
	}
      }
          
    outFile.close();
     
  }


  void Mesh::vtk_output_marker(std::string filename)const{

    std::string dir = "data/";
    std::string vtk = ".vtk";
    std::stringstream ss;
    ss << dir << filename << iter_counter << "_markers" << vtk ;
    std::string tmppath = ss.str();

    std::ofstream outFile;
    outFile.open(tmppath,std::ofstream::out);
    
    outFile.precision(8);
    outFile << "# vtk DataFile Version 2.0" << std::endl;
    outFile << "vtk output" << std::endl;
    outFile << "ASCII" << std::endl;
     
    outFile<< "DATASET POLYDATA" << std::endl;
    outFile << "POINTS" << " " << MLS_markers.size()<< " " << "float" << std::endl;
    for(int i = 0; i < MLS_markers.size(); i++){
      outFile<< MLS_markers[i].x_coord << "\t" << MLS_markers[i].y_coord << "\t" << 0.0 << std::endl;
      }
   
    outFile<<"POINT_DATA "<< MLS_markers.size()<< std::endl;
    outFile <<"SCALARS markers float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(int i = 0; i < MLS_markers.size(); i++){
      outFile<<  0.0 << std::endl;
    }
    
    outFile.close();

  }
  
  
  void Mesh::vtk_output_marker(std::vector<Particle> particle_vec, std::string filename)const{

    std::string dir = "data/";
    std::string vtk = ".vtk";
    std::stringstream ss;
    ss << dir << filename << iter_counter << "_markers" << vtk ;
    std::string tmppath = ss.str();

    std::ofstream outFile;
    outFile.open(tmppath,std::ofstream::out);
    
    outFile.precision(8);
    outFile << "# vtk DataFile Version 2.0" << std::endl;
    outFile << "vtk output" << std::endl;
    outFile << "ASCII" << std::endl;
     
    outFile<< "DATASET POLYDATA" << std::endl;
    outFile << "POINTS" << " " << particle_vec.size()<< " " << "float" << std::endl;
    for(int i = 0; i < particle_vec.size(); i++){
      outFile<< particle_vec[i].x_coord << "\t" << particle_vec[i].y_coord << "\t" << 0.0 << std::endl;
      }
   
    outFile<<"POINT_DATA "<< particle_vec.size()<< std::endl;
    outFile <<"SCALARS markers float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(int i = 0; i < particle_vec.size(); i++){
      outFile<<  0.0 << std::endl;
    }
    
    outFile.close();

  }

  void Mesh::correction1(){
    //std::cout <<"Inside correction function" << std::endl;
    std::vector<MLS::Particle> interface_nodes;
    std::unordered_map<std::pair<int,int>,double>  new_dist;

    Particle node_particle;
    int x_cell_index = 0;
    int y_cell_index = 0;
    double current_phi;
   
    double x_dist = 0;
    double y_dist = 0;
    double distance = 0;
    double x_quad,y_quad;

    for(auto &marker : MLS_markers){
       double new_phi;
      //-------------------UPDATE Distance PROCEDURE--------------------------------
      
      //for each node phi(i,j) near mk-assume it is within the same cell.
      //In the future- could use quadrants to update the closest neigbouring cells.
      x_cell_index = int(floor((marker.x_coord-x_min)/dx)) + nGhost;
      y_cell_index = int(floor((marker.y_coord-y_min)/dy)) + nGhost;
      
      x_dist = x_cell_axis(x_cell_index);
      y_dist = y_cell_axis(y_cell_index);
      
      std::cout << marker.x_coord << "\t" << marker.y_coord << std::endl;
      std::cout << x_dist << "\t" << y_dist << std::endl;
      std::cout << x_cell_index << "\t" << y_cell_index << std::endl;
      std::cout << dx*0.5 << "\t" << dy*0.5 << std::endl;
      
     
      for(int x_it = -1; x_it < 2; ++x_it){
	for(int y_it = -1; y_it < 2; ++y_it ){
	  std::pair<int,int> key(x_cell_index+x_it,y_cell_index+y_it); 
	  
	  current_phi = MLS_data(x_cell_index+x_it,y_cell_index+y_it).phi;

	  distance = sqrt((x_dist+x_it*dx- marker.x_coord)*(x_dist+x_it*dx - marker.x_coord) \
			  + (y_dist+y_it*dy - marker.y_coord)*(y_dist+y_it*dy - marker.y_coord));
	 	 
	  if(current_phi > 0){
	    double phi_pos = std::min(distance,fabs(MLS_data(x_cell_index+x_it,y_cell_index+y_it).phi));
	    new_phi = std::copysign(phi_pos,+1);
	    //MLS_data(x_cell_index+x_it,y_cell_index+y_it).phi = new_phi;
	   	   
	   if(current_phi != new_phi){
	     new_dist.insert(std::make_pair(key,new_phi));
	    }

	  }else if(current_phi < 0){
	    double phi_neg = std::min(distance,fabs(MLS_data(x_cell_index+x_it,y_cell_index+y_it).phi));
	    new_phi = std::copysign(phi_neg,-1);
	    //MLS_data(x_cell_index+x_it,y_cell_index+y_it).phi = new_phi;
	    if(current_phi != new_phi){
	      new_dist.insert(std::make_pair(key,new_phi));
	    }
	  }
	  
	   std::cout << x_dist+x_it*dx << "\t" << y_dist+y_it*dy << std::endl;
	   std::cout << distance << "\n";
	   std::cout << current_phi << std::endl;
	   std::cout <<"\n";
     

	}
      }
      std::cout <<"\n";
    }

    std::cout <<"The number of cells corrected has been " << new_dist.size() << std::endl;

    std::unordered_map<std::pair<int,int>,double>::iterator it;
     for(it = new_dist.begin(); it != new_dist.end(); ++it){
       int i = it->first.first;
       int j = it->first.second;
       MLS_data(i,j).phi = it->second; 
       // std::cout << "AT cell (" << i << ","<< j << ") with value " << it->second << std::endl;
       }
     new_dist.clear();
    
  }//End of correction fcn

  void Mesh::signCorrection(){
    
    //Loop through the interface map. Check surrounding cells to see if they are part of either +deltax or -delta x map. If yes calculate distance. 
    Z_area = 0;
    marchingSquares();
        
    std::unordered_map<std::pair<int,int>,std::vector<std::pair<Particle,Particle>>>::iterator it;
    //std::cout << "Correcting sign" << "\n";
    //std::cout << " There are : " << interface.bucket_count() << " interfacial cells" << "\n";
   
    for(it = interface.begin(); it != interface.end(); ++it){
      
	  int x_iter = it->first.first;
	  int y_iter = it->first.second; 

	  double node_x = x_cell_axis(x_iter);
	  double node_y = y_cell_axis(y_iter);
	  Particle node(node_x,node_y);
	  Particle p1,p2,p3,p4;

	  Particle pi1,pi2;
	  pi1.x_coord = it->second[0].first.x_coord;
	  pi1.y_coord = it->second[0].first.y_coord;
	  pi2.x_coord = it->second[0].second.x_coord;
	  pi2.y_coord = it->second[0].second.y_coord;	   ;
	   

	  double d_plus=1;
	  double d_plus_temp =1;//Distance from interface to plus delta x level set
	  double d_plus_amb = 1;
	  double d_minus = 1;
	  double d_minus_temp = 1;//Distance from interface to minus delta x level set
	  double d_minus_amb = 1;
      
	  double DLS_plus = 1;
	  double DLS_plus_temp = 1;
	  double DLS_minus = 1;
	  double DLS_minus_temp =1;
      
	  std::pair<int,int> min_key_cell; 
	  std::cout << "Correcting sign in cell " <<  x_iter << " " << y_iter << "\n";
	  std::cout << "Correcting sign at location " <<  x_cell_axis(x_iter) << " " << y_cell_axis(y_iter) << "\n";
      
	  //-------------------------------------------------------------------------------
	  for(int j = -2; j < 3; j++ ){
	    for(int i = -2; i < 3; i++){
	  
	      int x_cell = x_iter+i;
	      int y_cell = y_iter+j;
	  
	      std::pair<int,int> coords(x_cell,y_cell);
	      auto plus_it = plus_dx_LS.find(coords);
	 
	      if(plus_it == plus_dx_LS.end()){
		std::cout << " key not valid LS_plus" << "\n";
	      }else{
		std::cout <<"PLUS LS Valid key " << coords.first << "\t" << coords.second << "\n";
		if(plus_it->second.size() == 1){
		  p1.x_coord = plus_it->second[0].first.x_coord;
		  p1.y_coord = plus_it->second[0].first.y_coord;
		  p2.x_coord = plus_it->second[0].second.x_coord;
		  p2.y_coord = plus_it->second[0].second.y_coord;
		}else if(plus_it->second.size() == 2){
		  p1.x_coord = plus_it->second[0].first.x_coord;
		  p1.y_coord = plus_it->second[0].first.y_coord;
		  p2.x_coord = plus_it->second[0].second.x_coord;
		  p2.y_coord = plus_it->second[0].second.y_coord;
		  p3.x_coord = plus_it->second[1].first.x_coord;
		  p3.y_coord = plus_it->second[1].first.y_coord;
		  p4.x_coord = plus_it->second[1].second.x_coord;
		  p4.y_coord = plus_it->second[1].second.y_coord;
		}
	   	  
		d_plus_temp = minimum_distance(p1,p2,node);
		if(plus_it->second.size()==2){
		  d_plus_amb = minimum_distance(p3,p4,node);
		  d_plus_temp = std::min(d_plus_temp,d_plus_amb);
		}
		if(fabs(d_plus_temp) < fabs(d_plus)){
		  d_plus = d_plus_temp;
		  min_key_cell.first = x_cell;
		  min_key_cell.second = y_cell;
		}
	      }

	    }//End of LS_plus loop 
	  }//End of LS_plus loop 

	  
	  std::cout << "The LS_plus key is " << min_key_cell.first << " " << min_key_cell.second << "\n";
	  auto minimum_it = plus_dx_LS.find(min_key_cell);
	  if(minimum_it != plus_dx_LS.end()){
	  p1.x_coord = minimum_it->second[0].first.x_coord;
	  p1.y_coord = minimum_it->second[0].first.y_coord;
	  p2.x_coord = minimum_it->second[0].second.x_coord;
	  p2.y_coord = minimum_it->second[0].second.y_coord;
	  
	  Particle pproj = minimum_distance_point(p1,p2,node);
	  std::pair<Particle,Particle> S1(pi1,pi2);
	  std::pair<Particle,Particle> S2(node,pproj);
	  Particle i_point= inter_point(S1,S2);
	  DLS_plus_temp = dist(i_point,pproj);
	  
	  if(fabs(DLS_plus_temp) < fabs(DLS_plus)){
	    DLS_plus = DLS_plus_temp;
	    }	    
	  }
	  //----------------------------------------------------------------

	  //----------------------------------------------------------------
	  for(int j = -2; j < 3; j++ ){
	    for(int i = -2; i < 3; i++){
	      int x_cell = x_iter+i;
	      int y_cell = y_iter+j;
	  
	      std::pair<int,int> coords(x_cell,y_cell);
	      auto minus_it = minus_dx_LS.find(coords);

	      if(minus_it == minus_dx_LS.end()){
		std::cout << " key not valid LS_minus" << "\n";
	      }else{
		std::cout <<"MINUS LS Valid key " << coords.first << "\t" << coords.second << "\n";
		if(minus_it->second.size()==1){
		  p1.x_coord = minus_it->second[0].first.x_coord;
		  p1.y_coord = minus_it->second[0].first.y_coord;
		  p2.x_coord = minus_it->second[0].second.x_coord;
		  p2.y_coord = minus_it->second[0].second.y_coord;
		}else if (minus_it->second.size()==2){
		  p1.x_coord = minus_it->second[0].first.x_coord;
		  p1.y_coord = minus_it->second[0].first.y_coord;
		  p2.x_coord = minus_it->second[0].second.x_coord;
		  p2.y_coord = minus_it->second[0].second.y_coord;
		  p3.x_coord = minus_it->second[1].first.x_coord;
		  p3.y_coord = minus_it->second[1].first.y_coord;
		  p4.x_coord = minus_it->second[1].second.x_coord;
		  p4.y_coord = minus_it->second[1].second.y_coord;
		}
	  	    	  
		d_minus_temp = minimum_distance(p1,p2,node);
		if(minus_it->second.size()==2){
		  d_minus_amb = minimum_distance(p3,p4,node);
		  d_minus_temp = std::min(d_minus_temp,d_minus_amb);
		}
		if(fabs(d_minus_temp) < fabs(d_minus)){
		  d_minus = d_minus_temp;
		  min_key_cell.first = x_cell;
		  min_key_cell.second = y_cell;
		}
	    
	      }

	    }
	  }
	  std::cout << "The LS_minus key is " << min_key_cell.first << " " << min_key_cell.second << "\n";
	  auto min_it_minus = minus_dx_LS.find(min_key_cell);
	  if(min_it_minus != minus_dx_LS.end()){
	  
	    p1.x_coord = min_it_minus->second[0].first.x_coord;
	    p1.y_coord = min_it_minus->second[0].first.y_coord;
	    p2.x_coord = min_it_minus->second[0].second.x_coord;
	    p2.y_coord = min_it_minus->second[0].second.y_coord;
	    
	    Particle pproj_minus = minimum_distance_point(p1,p2,node);
	    
	    std::pair<Particle,Particle> S3(pi1,pi2);
	    std::pair<Particle,Particle> S4(node,pproj_minus);
	    Particle i_point_m= inter_point(S3,S4);
	    DLS_minus_temp = dist(i_point_m,pproj_minus);
	
	    if(fabs(DLS_minus_temp) < fabs(DLS_minus)){
	      DLS_minus = DLS_minus_temp;
	  }	    
	  //---------------------------------------------------------------------
	  }
	  /*
	  if(x_iter == 35 && y_iter == 130){
	      
	    std::cout << "\n";
	    std::cout << "Current iter is " << iter_counter << "\n";
	    std::cout << "At cell " << x_iter << " " << y_iter << "\n";
	    std::cout << "dx: " << dx << "\n";
	    std::cout << "The d_minus is : " << d_minus << " and d_plus: " << d_plus << " \n";
	    std::cout << "The distance between minus_dx_LS and 0LS is : " << DLS_minus << " \n";
	    std::cout << " The coeef is " << d_minus - DLS_minus << "\n";
	    std::cout << "The node is : " << node.x_coord << "," << node.y_coord << "\n";
	    std::cout << "The interface segment is decribed by point : " << "\n";
	    std::cout << "(" << pi1.x_coord << "," << pi1.y_coord << " )" << "  " << "(" << pi2.x_coord << "," << pi2.y_coord << " )" << "\n";
	    std::cout << " The LS_minus segment is defined by \n";
	    std::cout <<  "(" << p1.x_coord << "," << p1.y_coord << " )" << "  " << "(" << p2.x_coord << "," << p2.y_coord << " )" << "\n"; 
	  }
	  */


	  if(d_minus == 1 || d_plus == 1 ) {
	  }else {
	    double  coeff1,coeff2;
	    double phi1,phi2,phi3,phi4;
	    if(d_minus < d_plus){
	      coeff1 = d_minus - DLS_minus;
	      phi1 = MLS_data(x_iter,y_iter).phi;
	      //Overall result should be negattive independent of current sign value of data
	      MLS_data(x_iter,y_iter).phi=std::copysign( MLS_data(x_iter,y_iter).phi,coeff1);
	      phi2 = MLS_data(x_iter,y_iter).phi;
	
	      if (phi1 != phi2){
#ifdef DEBUG
		Sign_change_mark.push_back(node);
#endif
	      }
	
	    }else{
	      phi3 = MLS_data(x_iter,y_iter).phi;
	      coeff2 = DLS_plus-d_plus;
	      //Overall result should be positive independent of current sign value of data
	      MLS_data(x_iter,y_iter).phi=std::copysign( MLS_data(x_iter,y_iter).phi,coeff2);
	      phi4 = MLS_data(x_iter,y_iter).phi;
	      if (phi3 != phi4){
	 
#ifdef DEBUG
		Sign_change_mark.push_back(node);
#endif
		//	   if(x_iter == 138 && y_iter == 107){
		/*
		  std::cout << "\n";
		  std::cout << "Current iter is " << iter_counter << "\n";
		  std::cout << "There has been a sign change from :" << phi3 <<" to " << phi4 <<"\n";
		  std::cout << "At cell " << x_iter << " " << y_iter << "\n";
		  std::cout << "The coefficient is : " << coeff2 << "\n";
		  std::cout << "dx: " << dx << "\n";
		  std::cout << "The d1 is : " << d_minus << " and d2: " << d_plus << " \n";
		  std::cout << "The distance between plus_dx_LS and 0LS is : " << DLS_minus << " \n";
		  std::cout << "The node is : " << node.x_coord << "," << node.y_coord << "\n";
		  std::cout << "The interface segment is decribed by point : " << "\n";
		  std::cout << "(" << pi1.x_coord << "," << pi1.y_coord << " )" << "  " << "(" << pi2.x_coord << "," << pi2.y_coord << " )" << "\n";
		  std::cout << " The LS_plus segment is defined by \n";
		  std::cout <<  "(" << p1.x_coord << "," << p1.y_coord << " )" << "  " << "(" << p2.x_coord << "," << p2.y_coord << " )" << "\n"; 
		*/
	      }
	    }
	  }

    }//End of loop through interfacial cells; 

    interface.clear();
    plus_dx_LS.clear();
    minus_dx_LS.clear();

    std::string Snap_zero = "Snap_zero";
    std::string Snap_minus = "Snap_minus";
    std::string Snap_plus = "Snap_plus";
    vtk_output_marker(LS_zero,Snap_zero);
    vtk_output_marker(LS_minus,Snap_minus);
    vtk_output_marker(LS_plus,Snap_plus);
    LS_zero.clear();
    LS_minus.clear();
    LS_plus.clear();
#ifdef DEBUG
    std::string Snap_sign = "Snap_sign_check";
    vtk_output_marker(Sign_change_mark, Snap_sign);
    Sign_change_mark.clear();
#endif

 }
	  
   void Mesh::marchingSquares(){
    
    //std::cout << "Inside marchingSquares function" << std::endl;
    double zero_threshold = 0.0;
    double plus_delta_x_threshold = 1*dx;
    double minus_delta_x_threshold = -1*dx;
    int case_zero;
    int case_plus_x;
    int case_minus_x;
    std::vector<std::pair<Particle,Particle>> zero_pair,plus_pair,minus_pair;
    double vert1,vert2,vert3,vert4; // There are the corners of each cell name anticlockwise top right = 1.
    // std::cout << "Inside marching squares function" << std::endl;
    //Go through each of the squares.
    for(int j = nGhost; j < ncells+nGhost; j++){
      for(int i = nGhost; i < ncells+nGhost; i++){
	
	
	vert1 = 0.25*(MLS_data(i,j).phi+MLS_data(i+1,j).phi+MLS_data(i,j+1).phi+MLS_data(i+1,j+1).phi);
	vert2 = 0.25*(MLS_data(i,j).phi+MLS_data(i-1,j).phi+MLS_data(i,j+1).phi+MLS_data(i-1,j+1).phi);
	vert3 = 0.25*(MLS_data(i,j).phi+MLS_data(i-1,j).phi+MLS_data(i,j-1).phi+MLS_data(i-1,j-1).phi);
	vert4 = 0.25*(MLS_data(i,j).phi+MLS_data(i+1,j).phi+MLS_data(i,j-1).phi+MLS_data(i+1,j-1).phi);
	
	//Determine the case for each square
	case_zero = Case(zero_threshold, vert1,vert2,vert3,vert4);
	case_minus_x = Case(minus_delta_x_threshold, vert1,vert2,vert3,vert4);
	case_plus_x = Case(plus_delta_x_threshold, vert1,vert2,vert3,vert4);
	
	//Test for the case and set the key.
	zero_pair = squareReconstruction(case_zero,i,j,vert1,vert2,vert3,vert4,zero_threshold);
	minus_pair = squareReconstruction(case_minus_x,i,j,vert1,vert2,vert3,vert4,minus_delta_x_threshold);
	plus_pair = squareReconstruction(case_plus_x,i,j,vert1,vert2,vert3,vert4,plus_delta_x_threshold);
	
	std::pair<int,int> key_index(i,j);
	if(case_zero != 0 && case_zero !=15){
	  interface.insert({key_index,zero_pair});
	}
	if(case_minus_x != 0 && case_minus_x != 15){
	  minus_dx_LS.insert({key_index,minus_pair});
	}
	if(case_plus_x != 0 && case_plus_x != 15){
	  plus_dx_LS.insert({key_index,plus_pair});
	}

	}
    }

  }//End of fcn

  int Mesh::Case(double threshold, double vert1,double vert2,double vert3,double vert4){
   
    int case_output = 0;
    
    if(vert2 > threshold){ case_output += 1 ;}
    if(vert1 > threshold){ case_output += 2 ;}
    if(vert3 > threshold){ case_output += 8 ;}  
    if(vert4 > threshold){ case_output += 4 ;}
  
   
    return case_output;
  }

  std::vector<std::pair<Particle,Particle>> Mesh::squareReconstruction(int Case, int x_iter ,int y_iter,double v1,double v2,double v3,double v4,double threshold){

    //std::cout <<"Inside squareReconstruction" << std::endl;
    std::vector<std::pair<Particle,Particle>> output;
    //the first member of the pair is the gradient of the line, the second one a point on the line.
    std::pair<Particle,Particle> values;
    double a1,a2,a3,a4;
    Particle p1,p2,p3,p4;

    if(Case == 0){
      //std::cout << "This is case 0 " << "\n";
      output.resize(1);
      p1.x_coord = -1;
      p1.y_coord = -1;
      values.first = p1;
      values.second = p1;
      output.push_back(values);
      if(threshold == 0){
      Z_area += 0.0;
      }
    }else if(Case == 15){
      //std::cout << " Case fifteen" << "\n";
      output.resize(1);
      p1.x_coord = -1;
      p1.y_coord = -1;
      values.first = p1;
      values.second = p1;
      output.push_back(values);
      if(threshold == 0){      Z_area += dx*dy;}
    }else if(Case == 1){
      //std::cout << " Case one" << "\n";
      a1 = (v2-threshold)/(v2-v1);
      a2 = (v2-threshold)/(v2-v3);
       
      p1.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a1*dx;
      p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);
      
      p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
      p2.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a2*dy;
          
      values.first = p1;
      values.second = p2;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      if(threshold == 0){      Z_area += dx*dy*0.5*a1*a2;}

      if(threshold == 0){
	LS_zero.push_back(p1);
	LS_zero.push_back(p2);
      }else if(threshold == dx){
	LS_plus.push_back(p1);
	LS_plus.push_back(p2);
      }else if(threshold == -dx){
	LS_minus.push_back(p1);
	LS_minus.push_back(p2);
      }
      
    }else if(Case == 2){
      //std::cout << " Case two" << "\n";
      a1 = (v1-threshold)/(v1-v2);
      a4 = (v1-threshold)/(v1-v4);
       
      p1.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a1*dx;
      p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);

      p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
      p4.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a4*dy;
                  
      values.first = p1;
      values.second = p4;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      
      if(threshold == 0){ Z_area += dx*dy*0.5*a1*a4;}
      
      if(threshold == 0){
	LS_zero.push_back(p1);
	LS_zero.push_back(p4);
      }else if(threshold == dx){
	LS_plus.push_back(p1);
	LS_plus.push_back(p4);
      }else if(threshold == -dx){
	LS_minus.push_back(p1);
	LS_minus.push_back(p4);
      }
      

    }else if(Case == 3){
      //std::cout << " Case three" << "\n";
      a2 = (v2-threshold)/(v2-v3);
      a4 = (v1-threshold)/(v1-v4);
       
      p1.x_coord = (x_cell_axis(x_iter)+0.5*dx);
      p1.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a4*dy;

      p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
      p2.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a2*dy;
                  
      values.first = p1;
      values.second = p2;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      
      if(threshold == 0){ Z_area += dx*dy*(a2 + 0.5*(a4-a2));}
      
      if(threshold == 0){
	LS_zero.push_back(p1);
	LS_zero.push_back(p2);
      }else if(threshold == dx){
	LS_plus.push_back(p1);
	LS_plus.push_back(p2);
      }else if(threshold == -dx){
	LS_minus.push_back(p1);
	LS_minus.push_back(p2);
      }
      
      
    }else if(Case == 4){
      //std::cout << " Case four" << "\n";
      a4 = (v4-threshold)/(v4-v1);
      a3 = (v4-threshold)/(v4-v3);
       
      p3.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a3*dx;
      p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);

      p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
      p4.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a4*dy;
                  
      values.first = p3;
      values.second = p4;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      
      if(threshold == 0){ Z_area += dx*dy*0.5*a3*a4;}
      
      if(threshold == 0){
	LS_zero.push_back(p3);
	LS_zero.push_back(p4);
      }else if(threshold == dx){
	LS_plus.push_back(p3);
	LS_plus.push_back(p4);
      }else if(threshold == -dx){
	LS_minus.push_back(p3);
	LS_minus.push_back(p4);
      }
      
    }else if(Case == 5){
      
      //std::cout << "In AMBIGUOUS case five." << "\n";
       if(threshold == 0){
	 std::cout << "Ambiguous case five for zero level set \n";
	 std::cout << "Iteration is : " << iter_counter << "\n";
	 std::cout << "The cell is : " << x_iter << " " << y_iter << "\n";
	 std::cout << "the coords are : " << x_cell_axis(x_iter) << " " << y_cell_axis(y_iter) << "\n";
       }else if(threshold == dx){
	 std::cout << "Ambiguous case five for plus delta x level set \n";
	 std::cout << "Iteration is : " << iter_counter << "\n";
	 std::cout << "The cell is : " << x_iter << " " << y_iter << "\n";
	 std::cout << "the coords are : " << x_cell_axis(x_iter) << " " << y_cell_axis(y_iter) << "\n";
       }else if(threshold == -dx){
	 std::cout << "Ambiguous case five for minus delta x level set \n";
	 std::cout << "Iteration is : " << iter_counter << "\n";
	 std::cout << "The cell is : " << x_iter << " " << y_iter << "\n";
	 std::cout << "the coords are : " << x_cell_axis(x_iter) << " " << y_cell_axis(y_iter) << "\n";
       }

      if( MLS_data(x_iter,y_iter).phi > 0.0 ){
	a2 = (threshold-v3)/(v2-v3);
	a3 = (threshold-v3)/(v4-v3);
       
	p1.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a3*dx;
	p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
	p1.y_coord = (y_cell_axis(y_iter)-0.5*dy);
	p2.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a2*dy;
                  
	values.first = p1;
	values.second = p2;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);

	a1 = (threshold-v1)/(v2-v1);
	a4 = (threshold-v1)/(v4-v1);
       
	p3.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a1*dx;
	p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
	p3.y_coord = (y_cell_axis(y_iter)+0.5*dy);
	p4.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a4*dy;
                  
	values.first = p3;
	values.second = p4;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);
	
      if(threshold == 0){
	Z_area += dx*dy*(1-0.5*(1-a1)*(1-a4)-0.5*(1-a2)*(1-a3));
      }
      /*
	std::cout << "Case 5-stripe " << "\n";
	std::cout << "The vertices are " << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";
	std::cout << "The value at the cell centre is : " << MLS_data(x_iter,y_iter).phi << "\n";
	std::cout << "The weights are " << a1 << " "<< a2 << " "<< a3 << " " << a4 << "\n";
	std::cout << "The cell centre is " << x_cell_axis(x_iter)<< " "<< y_cell_axis(y_iter) << "\n";
	std::cout << "p1 : " << p1.x_coord << " " << p1.y_coord << "\n";
	std::cout << "p2 : " << p2.x_coord << " " << p2.y_coord << "\n";
	std::cout << "p3 : " << p3.x_coord << " " << p3.y_coord << "\n";
	std::cout << "p4 : " << p4.x_coord << " " << p4.y_coord << "\n";
	*/
      }else if(MLS_data(x_iter,y_iter).phi <  0.0){
	
	//----------------------------------------------------
	a1 = (v2-threshold)/(v2-v1);
	a2 = (v2-threshold)/(v2-v3);
       
	p1.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a1*dx;
	p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
	p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);
	p2.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a2*dy;
                  
	values.first = p1;
	values.second = p2;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);
	//---------------------------------------------------
	a3 = (v4-threshold)/(v4-v1);
	a4 = (v4-threshold)/(v4-v3);
       
	p3.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a3*dx;
	p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
	p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);
	p4.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a4*dy;
                  
	values.first = p3;
	values.second = p4;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);

	if(threshold == 0){Z_area += dx*dy*(0.5*a1*a2 + 0.5*a3*a4);}
	/*
	std::cout << "Case 5- two corners " << "\n";
	std::cout << "The vertices are " << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";
	std::cout << "The value at the cell centre is : " << MLS_data(x_iter,y_iter).phi << "\n";
	std::cout << "The weights are " << a1 << " "<< a2 << " "<< a3 << " " << a4 << "\n";
	std::cout << "The cell centre is " << x_cell_axis(x_iter)<< " "<< y_cell_axis(y_iter) << "\n";
	std::cout << "p1 : " << p1.x_coord << " " << p1.y_coord << "\n";
	std::cout << "p2 : " << p2.x_coord << " " << p2.y_coord << "\n";
	std::cout << "p3 : " << p3.x_coord << " " << p3.y_coord << "\n";
	std::cout << "p4 : " << p4.x_coord << " " << p4.y_coord << "\n";
	*/
      }


    }else if(Case == 6){
      //std::cout << " case six" << "\n";
      a1 = (v1-threshold)/(v1-v2);
      a3 = (v4-threshold)/(v4-v3);
       
      p1.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a1*dx;
      p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);

      p3.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a3*dx;
      p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);
                  
      values.first = p1;
      values.second = p3;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);

      if(threshold == 0){Z_area += dx*dy*(a3 +0.5*(a1-a3));}
      
      if(threshold == 0){
	LS_zero.push_back(p1);
	LS_zero.push_back(p3);
      }else if(threshold == dx){
	LS_plus.push_back(p1);
	LS_plus.push_back(p3);
      }else if(threshold == -dx){
	LS_minus.push_back(p1);
	LS_minus.push_back(p3);
      }
      
      
    }else if(Case == 7){
      //std::cout << " Case seven" << "\n";
      a2 = (v2-threshold)/(v2-v3);
      a3 = (v4-threshold)/(v4-v3);
       
      p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
      p2.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a2*dy;

      p3.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a3*dx;
      p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);
                  
      values.first = p2;
      values.second = p3;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      
      if(threshold == 0){Z_area += dx*dy*(1-0.5*(1-a2)*(1-a3));}
      
      if(threshold == 0){
	LS_zero.push_back(p3);
	LS_zero.push_back(p2);
      }else if(threshold == dx){
	LS_plus.push_back(p3);
	LS_plus.push_back(p2);
      }else if(threshold == -dx){
	LS_minus.push_back(p3);
	LS_minus.push_back(p2);
      }
      
      /*
      std::cout << "Case seven " << "\n";
      std::cout << "The vertices are " << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";
      std::cout << "The weights are " << a2 << " " << a4 << "\n";
      std::cout << "The cell centre is " << x_cell_axis(x_iter)<< " "<< y_cell_axis(y_iter) << "\n";
      std::cout << "p1 : " << p1.x_coord << " " << p1.y_coord << "\n";
      std::cout << "p2 : " << p2.x_coord << " " << p2.y_coord << "\n";
      std::cout << "The gradient is " << grad << " the norm : " << norm << "\n";
      */
    }else if(Case == 8){
      //std::cout << " Case eight" << "\n";
      a2 = (v3-threshold)/(v3-v2);
      a3 = (v3-threshold)/(v3-v4);
       
      p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
      p2.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a2*dy;
    
      p3.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a3*dx;
      p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);
                  
      values.first = p2;
      values.second = p3;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      if(threshold == 0){Z_area += dx*dy*(0.5*a2*a3);}
      
      if(threshold == 0){
	LS_zero.push_back(p3);
	LS_zero.push_back(p2);
      }else if(threshold == dx){
	LS_plus.push_back(p3);
	LS_plus.push_back(p2);
      }else if(threshold == -dx){
	LS_minus.push_back(p3);
	LS_minus.push_back(p2);
      }
      
    }else if(Case == 9){
      //std::cout << " Case nine" << "\n";
      a1 = (v2-threshold)/(v2-v1);
      a3 = (v3-threshold)/(v3-v4);
       
      p1.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a1*dx;
      p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);

      p3.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a3*dx;
      p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);
                  
      values.first = p1;
      values.second = p3;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      

      if(threshold == 0){ Z_area += dx*dy*(a1 -0.5*(a3-a1));}
      if(threshold == 0){
	LS_zero.push_back(p1);
	LS_zero.push_back(p3);
      }else if(threshold == dx){
	LS_plus.push_back(p1);
	LS_plus.push_back(p3);
      }else if(threshold == -dx){
	LS_minus.push_back(p1);
	LS_minus.push_back(p3);
      }
      
    }else if(Case == 10){
      std::cout << "In ambiguous case 10 " << "\n";
       if(threshold == 0){
	 std::cout << "Ambiguous case ten for zero level set \n";
	 std::cout << "Iteration is : " << iter_counter << "\n";
	 std::cout << "The cell is : " << x_iter << " " << y_iter << "\n";
	 std::cout << "the coords are : " << x_cell_axis(x_iter) << " " << y_cell_axis(y_iter) << "\n";
       }else if(threshold == dx){
	 std::cout << "Ambiguous case ten for plus delta x level set \n";
	 std::cout << "Iteration is : " << iter_counter << "\n";
	 std::cout << "The cell is : " << x_iter << " " << y_iter << "\n";
	 std::cout << "the coords are : " << x_cell_axis(x_iter) << " " << y_cell_axis(y_iter) << "\n";
       }else if(threshold == -dx){
	 std::cout << "Ambiguous case ten for minus delta x level set \n";
	 std::cout << "Iteration is : " << iter_counter << "\n";
	 std::cout << "The cell is : " << x_iter << " " << y_iter << "\n";
	 std::cout << "the coords are : " << x_cell_axis(x_iter) << " " << y_cell_axis(y_iter) << "\n";
       }
      if(MLS_data(x_iter,y_iter).phi > 0.0){
	a1 = (threshold-v2)/(v1-v2);
	a2 = (threshold-v2)/(v4-v2);
       
	p1.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a1*dx;
	p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
	p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);
	p2.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a2*dy;

	values.first = p1;
	values.second = p2;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);
      
	a3 = (threshold-v4)/(v3-v4);
	a4 = (threshold-v4)/(v1-v4);
       
	p3.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a3*dx;
	p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
	p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);
	p4.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a4*dy;
                  
	values.first = p3;
	values.second = p4;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);
      if(threshold == 0){
	Z_area += dx*dy*(1-0.5*((1-a1)*(1-a2))-0.5*((1-a3)*(1-a4)));
      }	
	/*
	std::cout << "Case 10 cenral stripe " << "\n";
	std::cout << "The vertices are " << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";
	std::cout << "The value at the cell centre is : " << MLS_data(x_iter,y_iter).phi << "\n";
	std::cout << "The weights are " << a1 << " "<< a2 << " "<< a3 << " " << a4 << "\n";
	std::cout << "The cell centre is " << x_cell_axis(x_iter)<< " "<< y_cell_axis(y_iter) << "\n";
	std::cout << "p1 : " << p1.x_coord << " " << p1.y_coord << "\n";
	std::cout << "p2 : " << p2.x_coord << " " << p2.y_coord << "\n";
	std::cout << "p3 : " << p3.x_coord << " " << p3.y_coord << "\n";
	std::cout << "p4 : " << p4.x_coord << " " << p4.y_coord << "\n";
	*/
	
      }else if(MLS_data(x_iter,y_iter).phi < 0.0){
	
	a1 = (v1-threshold)/(v1-v2);
	a4 = (v1-threshold)/(v1-v4);
       
	p1.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a1*dx;
	p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
	p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);
	p4.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a4*dy;

	values.first = p1;
	values.second = p4;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);
      
	a2 = (v3-threshold)/(v3-v2);
	a3 = (v3-threshold)/(v3-v4);
       
	p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
	p2.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a3*dy;

	p3.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a4*dx;
	p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);
                  
	values.first = p3;
	values.second = p4;//This is a particle object not a double -- need to change signature of function. 
	output.push_back(values);
	if(threshold == 0){Z_area += dx*dy*(0.5*a1*a4+0.5*a2*a3);}
	/*
	std::cout << "Case 10- two corners " << "\n";
	std::cout << "The vertices are " << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";
	std::cout << "The value at the cell centre is : " << MLS_data(x_iter,y_iter).phi << "\n";
	std::cout << "The weights are " << a1 << " "<< a2 << " "<< a3 << " " << a4 << "\n";
	std::cout << "The cell centre is " << x_cell_axis(x_iter)<< " "<< y_cell_axis(y_iter) << "\n";
	std::cout << "p1 : " << p1.x_coord << " " << p1.y_coord << "\n";
	std::cout << "p2 : " << p2.x_coord << " " << p2.y_coord << "\n";
	std::cout << "p3 : " << p3.x_coord << " " << p3.y_coord << "\n";
	std::cout << "p4 : " << p4.x_coord << " " << p4.y_coord << "\n";
	*/
      }
      
    }else if(Case == 11){
      //std::cout << " Case eleven" << "\n";
      a3 = (v3-threshold)/(v3-v4);
      a4 = (v1-threshold)/(v1-v4);
       
      p3.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a3*dx;
      p3.y_coord = (y_cell_axis(y_iter)-0.5*dy);

      p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
      p4.y_coord = (y_cell_axis(y_iter)+0.5*dy)-a4*dy;
                  
      values.first = p3;
      values.second = p4;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      if(threshold == 0){ Z_area += dx*dy*(1-0.5*(1-a4)*(1-a3));}

      if(threshold == 0){
	LS_zero.push_back(p4);
	LS_zero.push_back(p3);
      }else if(threshold == dx){
	LS_plus.push_back(p4);
	LS_plus.push_back(p3);
      }else if(threshold == -dx){
	LS_minus.push_back(p4);
	LS_minus.push_back(p3);
      }
      


      //std::cout << "In case eleven: " << "\n";
    }else if(Case == 12){

      a2 = (v3-threshold)/(v3-v2);
      a4 = (v4-threshold)/(v4-v1);
       
      p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
      p2.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a2*dy;

      p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
      p4.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a4*dy;
                  
      values.first = p2;
      values.second = p4;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      if(threshold == 0){ Z_area += dx*dy*(a2 + 0.5*(a4-a2));}
      if(threshold == 0){
	LS_zero.push_back(p2);
	LS_zero.push_back(p4);
      }else if(threshold == dx){
	LS_plus.push_back(p2);
	LS_plus.push_back(p4);
      }else if(threshold == -dx){
	LS_minus.push_back(p2);
	LS_minus.push_back(p4);
      }
      
    }else if(Case == 13){
      //std::cout << " Case thirteen" << "\n";
      a1 = (v2-threshold)/(v2-v1);
      a4 = (v4-threshold)/(v4-v1);
       
      p1.x_coord = (x_cell_axis(x_iter)-0.5*dx)+a1*dx;
      p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);      
      
      p4.x_coord = (x_cell_axis(x_iter)+0.5*dx);
      p4.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a4*dy;
                  
      values.first = p1;
      values.second = p4;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      if(threshold == 0){ Z_area += dx*dy*(1-0.5*(1-a1)*(1-a4));}
      
      if(threshold == 0){
	LS_zero.push_back(p1);
	LS_zero.push_back(p4);
      }else if(threshold == dx){
	LS_plus.push_back(p1);
	LS_plus.push_back(p4);
      }else if(threshold == -dx){
	LS_minus.push_back(p1);
	LS_minus.push_back(p4);
      }
      
    }else if(Case == 14){
      //std::cout << " Case fourteen" << "\n";

      a1 = (v1-threshold)/(v1-v2);
      a2 = (v3-threshold)/(v3-v2);
       
      p1.x_coord = (x_cell_axis(x_iter)+0.5*dx)-a1*dx;
      p1.y_coord = (y_cell_axis(y_iter)+0.5*dy);

      p2.x_coord = (x_cell_axis(x_iter)-0.5*dx);
      p2.y_coord = (y_cell_axis(y_iter)-0.5*dy)+a2*dy;
                  
      values.first = p1;
      values.second = p2;//This is a particle object not a double -- need to change signature of function. 
      output.push_back(values);
      if(threshold == 0){ Z_area +=dx*dy*(1-0.5*(1-a1)*(1-a2));}
      
      if(threshold == 0){
	LS_zero.push_back(p1);
	LS_zero.push_back(p2);
      }else if(threshold == dx){
	LS_plus.push_back(p1);
	LS_plus.push_back(p2);
      }else if(threshold == -dx){
	LS_minus.push_back(p1);
	LS_minus.push_back(p2);
      }
      
    }
    //std::cout <<"The current area is: " << Z_area << "\n";
    return output;
    
  }

  double Mesh::minimum_distance (Particle p1, Particle p2, Particle node){
    //Return minimum distance between line segment p1p2 and the point node 
    Particle a,b;
    const double l2 = distance_squared(p1,p2);
    //std::cout << "\n Length squared " << l2 << "\n";
    a.x_coord = node.x_coord-p1.x_coord;
    a.y_coord = node.y_coord-p1.y_coord;

    b.x_coord = p2.x_coord-p1.x_coord;
    b.y_coord = p2.y_coord-p1.y_coord;

    //Consider the line extending the segment parameterized as v+t(w-v)
    //We find the projection of point p onto the line
    //It falls where t = [(p-v).(W-V)]/|w-v|^2;
    const double t = dot_prod(a,b)/l2;
    //std::cout << " t factor " << t << "\n";
    if(t < 0.0) return dist(node,p1);//Beyond the "v" end of the segment
    if(t > 1.0) return dist(node,p2);//Beyond the "w" end of the segment
    Particle projection;
    projection.x_coord = p1.x_coord + t *(p2.x_coord-p1.x_coord);
    projection.y_coord = p1.y_coord + t *(p2.y_coord-p1.y_coord);
    return dist(node,projection);
  }
  
  double Mesh::dist(Particle p1, Particle p2){
    return sqrt(distance_squared(p1,p2)); 
  
  }
  
  double Mesh::distance_squared(Particle p1, Particle p2){
    return ((p1.x_coord-p2.x_coord)*(p1.x_coord-p2.x_coord) + (p1.y_coord-p2.y_coord)*(p1.y_coord-p2.y_coord)); 
  }

  double Mesh::dot_prod(Particle p1,Particle p2){
    return (p1.x_coord*p2.x_coord + p1.y_coord*p2.y_coord);
  }

  

  double Mesh::dist2D_Segment_to_Segment( std::pair<Particle,Particle> S1, std::pair<Particle,Particle> S2)
{
  double SMALL_NUM= 0.00000001; // anything that avoids division overflow
  //std::cout << "In distance function \n";
  Particle p0,p1,q0,q1;
  p0 = S1.first;
  p1 = S1.second;
  q0 = S2.first;
  q1 = S2.second;
  Particle   u;
  u.x_coord = p1.x_coord - p0.x_coord;
  u.y_coord = p1.y_coord - p0.y_coord;
 
  //Vector   v = S2.P1 - S2.P0;
  Particle v;
  v.x_coord = q1.x_coord-q0.x_coord;
  v.y_coord = q1.y_coord-q0.y_coord;
  
  //  Vector   w = S1.P0 - S2.P0;
  Particle w;
  w.x_coord = p0.x_coord - q0.x_coord;
  w.y_coord = p0.y_coord - q0.y_coord;

  double    a = dot_prod(u,u);         // always >= 0
  double    b = dot_prod(u,v);
  double    c = dot_prod(v,v);         // always >= 0
  double    d = dot_prod(u,w);
  double    e = dot_prod(v,w);
  double    D = a*c - b*b;        // always >= 0
  double    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
  double    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0
  
    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
	//std::cout << "Parallel lines \n";
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
	    //std::cout << " S=0 tc=t_0 \n";
        }else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
	    //std::cout << " S=1 tc=t_0 \n";
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0){
	  sN = 0.0;
	  //std::cout << " Sc=0 Tc=0 \n";
        }else if (-d > a){
	  sN = sD;
	  //std::cout << " Sc=1 Tc=0 \n";
	}else{
	  sN = -d;
	  sD = a;
	  //std::cout << " Sc=S_0 Tc=0 \n";
        }
    }else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0){
            sN = 0;
	    //std::cout << " Sc=0 Tc=1 \n";
        }else if ((-d + b) > a){
            sN = sD;
	    //std::cout << " Sc=1 Tc=1 \n";
        }else {
            sN = (-d +  b);
            sD = a;
	    //std::cout << " Sc=s_0 Tc=1 \n";
	}
    }
    
    // finally do the division to get sc and tc
    sc = (fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);
    //std::cout << "Final value Sc = " << sc << " and Tc = " << tc << "\n";

    // get the difference of the two closest points
    Particle dP;
    dP.x_coord = w.x_coord + (sc * u.x_coord) - (tc * v.x_coord);  // =  S1(sc) - S2(tc)
    dP.y_coord = w.y_coord + (sc * u.y_coord) - (tc * v.y_coord);  // =  S1(sc) - S2(tc)
    //Particle Origin (0,0);
    
    return sqrt(dot_prod(dP,dP));   // return the closest distance
}

 Particle Mesh::inter_point(std::pair<Particle,Particle> S1, std::pair<Particle,Particle> S2){

   double SMALL_NUM= 0.00000001; // anything that avoids division overflow
   //std::cout << "In distance function \n";
   Particle p0,p1; //Segment describing zero Ls
   Particle q0,q1;//Segment describing node-p_projection
   p0 = S1.first;
   p1 = S1.second;
   q0 = S2.first;
   q1 = S2.second;
   Particle   u;
   u.x_coord = p1.x_coord - p0.x_coord;
   u.y_coord = p1.y_coord - p0.y_coord;
 
   //Vector   v = S2.P1 - S2.P0;
   Particle v;
   v.x_coord = q1.x_coord-q0.x_coord;
   v.y_coord = q1.y_coord-q0.y_coord;
  
   //  Vector   w = S1.P0 - S2.P0;
   Particle w;
   w.x_coord = p0.x_coord - q0.x_coord;
   w.y_coord = p0.y_coord - q0.y_coord;

    
   double D = perp(u,v);
   //std::cout << "The perp product is " << D << "\n";

  
   double distance = 0; 
   // test if  they are parallel (includes either being a point)
   if (fabs(D) < SMALL_NUM) {           // S1 and S2 are parallel
     //std::cout << "The lines are parallel \n";
     Particle origin(100,100);
     return origin;
   }
  
   // the segments are skew and may intersect in a point
   // get the intersect parameter for S1
   double     sI = perp(v,w) / D;
   double     tI = perp(u,w) / D;

   //std::cout << " sI is " << sI << " and tI is " << tI << "\n";

   Particle P_s;
   P_s.x_coord = p0.x_coord + sI*u.x_coord;  
   P_s.y_coord = p0.y_coord + sI*u.y_coord;  
  
   Particle Q_t;
   Q_t.x_coord = q0.x_coord + tI*v.x_coord;  
   Q_t.y_coord = q0.y_coord + tI*v.y_coord; 

   //std::cout << "The location of the intersection point is: " << P_s.x_coord << "," << P_s.y_coord << "\n";
  
  

   if(P_s.x_coord == Q_t.x_coord && P_s.x_coord == Q_t.x_coord){
    
     //std::cout << "Point coindice" <<std::endl;
     //Calculate distance.Distance between point on LS and Point of intersection with LS 0.  
     return P_s;
   }
 }
 //===================================================================

 Particle Mesh::minimum_distance_point (Particle p1, Particle p2, Particle node){
   //Return minimum distance between line segment p1p2 and the point node 
   Particle a,b;
   const double l2 = distance_squared(p1,p2);
   ////std::cout << "\n Length squared " << l2 << "\n";
   a.x_coord = node.x_coord-p1.x_coord;
   a.y_coord = node.y_coord-p1.y_coord;

   b.x_coord = p2.x_coord-p1.x_coord;
   b.y_coord = p2.y_coord-p1.y_coord;

   //Consider the line extending the segment parameterized as v+t(w-v)
   //We find the projection of point p onto the line
   //It falls where t = [(p-v).(W-V)]/|w-v|^2;
   const double t = dot_prod(a,b)/l2;
   ////std::cout << " t factor " << t << "\n";
   if(t < 0.0) return p1;//Beyond the "v" end of the segment
   if(t > 1.0) return p2;//Beyond the "w" end of the segment
   Particle projection;
   projection.x_coord = p1.x_coord + t *(p2.x_coord-p1.x_coord);
   projection.y_coord = p1.y_coord + t *(p2.y_coord-p1.y_coord);
   return projection;
 }
  

 
 double Mesh::perp(Particle u, Particle v){
   return (u.x_coord * v.y_coord - u.y_coord * v.x_coord);
 }
 

  void Mesh::advect_RK_HOUC(){

    //------------------------------------------------------------------


    blitz::Array<double, 2> phi_n; 
    blitz::Array<double, 2> phi_n_one; 
    blitz::Array<double, 2> phi_n_two; 
    blitz::Array<double, 2> phi_n_half; 
    blitz::Array<double, 2> phi_n_three_half;

    blitz::Array<double, 2> V_phi_n; 
    blitz::Array<double, 2> V_phi_n_one; 
    blitz::Array<double, 2> V_phi_n_half; 

    phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_two.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    phi_n_three_half.resize(2*nGhost+ncells,2*nGhost+ncells);
    
    V_phi_n.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_one.resize(2*nGhost+ncells,2*nGhost+ncells);
    V_phi_n_half.resize(2*nGhost+ncells,2*nGhost+ncells);

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	phi_n(row,col)= MLS_data(row,col).phi;
      }
    }

    if(iter_counter%2 == 0){
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------

      blitz::Array<double,2> phi_n_one_star;
      phi_n_one_star.resize(2*nGhost+ncells,2*nGhost+ncells);

      V_phi_n = spatial_HOUC_X(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_HOUC_Y(phi_n_one_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n_one_star(row,col) - dt*V_phi_n(row,col);
	}
      }

      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
      blitz::Array<double,2> phi_n_two_star;
      phi_n_two_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_one = spatial_HOUC_X(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_HOUC_Y(phi_n_two_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_two_star(row,col) - dt*V_phi_n_one(row,col);
	}
      }
      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }
      
      //Calculate phi^{n+3/2}
      //----------------------------------------------------------------------------
      blitz::Array<double,2> phi_n_three_half_star;
      phi_n_three_half_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_half = spatial_HOUC_X(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_HOUC_Y(phi_n_three_half_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_three_half_star(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }

    }else if(iter_counter%2 == 1){
       
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------

      blitz::Array<double,2> phi_n_one_star;
      phi_n_one_star.resize(2*nGhost+ncells,2*nGhost+ncells);

      V_phi_n = spatial_HOUC_Y(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_HOUC_X(phi_n_one_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = phi_n_one_star(row,col) - dt*V_phi_n(row,col);
	}
      }

      //-----------------------------------------------------------------


      //Calculate phi^(n+2)
      //-----------------------------------------------------------------
      blitz::Array<double,2> phi_n_two_star;
      phi_n_two_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_one = spatial_HOUC_Y(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_HOUC_X(phi_n_two_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two(row,col) = phi_n_two_star(row,col) - dt*V_phi_n_one(row,col);
	}
      }
      //------------------------------------------------------------------

      double a_coeff =  0.75;
      double b_coeff = 0.25;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_half(row,col) = a_coeff*phi_n(row,col) + b_coeff*phi_n_two(row,col) ;
	}
      }

      //----------------------------------------------------------------------------
      blitz::Array<double,2> phi_n_three_half_star;
      phi_n_three_half_star.resize(2*nGhost+ncells,2*nGhost+ncells);
    
      V_phi_n_half = spatial_HOUC_Y(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_HOUC_X(phi_n_three_half_star);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half(row,col) = phi_n_three_half_star(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      //----------------------------------------------------------------------------
      double c_coeff = 1.0/3.0;
      double d_coeff = 2.0/3.0;

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one(row,col) = c_coeff*phi_n(row,col) + d_coeff*phi_n_three_half(row,col) ;
	}
      }
    }

    for(int row = nGhost; row < nGhost+ncells; ++row){
      for(int col = nGhost; col < nGhost+ncells; ++col){
	MLS_data(row,col).phi= phi_n_one(row,col);
      }
    }


    //-----------------------------------------------------------------
   
  }

    blitz::Array<double,2> Mesh::spatial_HOUC_Y(blitz::Array<double,2> input ){

    // std::cout << "ADVECTING" << "\n";
    applyBC(input);
    double phi;
    double phi_ydir;
    blitz::Array<double,2>phi_new;
    double obj_speed_y;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_y = this->speed_y(x_cell_axis(col),y_cell_axis(row),time,T_max);
	
	phi = 0;
      
	if(obj_speed_y < 0){
	  stencil = "right";
	  dir = "y_dir";
	  phi_ydir = HOUC(col,row,stencil,dir,input);
	}else if(obj_speed_y == 0){
	  phi_ydir = 0.0;
	}else if(obj_speed_y > 0){
	  stencil = "left";
	  dir = "y_dir";
	  phi_ydir = HOUC(col,row,stencil,dir,input);
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	
	phi = obj_speed_y*phi_ydir;
	phi_new(col,row) = phi;
      }
    }

    applyBC(phi_new);
    return phi_new;

  }

  blitz::Array<double,2> Mesh::spatial_HOUC_X(blitz::Array<double,2> input ){

    // std::cout << "ADVECTING" << "\n";
    applyBC(input);
    double phi;
    double phi_xdir;
    blitz::Array<double,2>phi_new;
    double obj_speed_x;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_x = this->speed_x(x_cell_axis(col),y_cell_axis(row),time,T_max);
	
	phi = 0;

	if(obj_speed_x < 0){
	  stencil = "right";
	  dir = "x_dir";
	  phi_xdir = HOUC(col,row,stencil,dir,input); 
	}else if(obj_speed_x == 0){
	   phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  stencil = "left";
	  dir = "x_dir";
	  phi_xdir = HOUC(col,row,stencil,dir,input); 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
	
	phi = obj_speed_x*phi_xdir;
	phi_new(col,row) = phi;
      }
    }
    applyBC(phi_new);
    return phi_new;

  }

  double Mesh::HOUC(int i, int j, std::string stencil, std::string dir, const blitz::Array<double,2> & input){
    double phi_output;
    if(dir == "y_dir" && stencil == "right"){
        phi_output = (2*input(i,j+3) - 15*input(i,j+2) + 60*input(i,j+1) - 20*input(i,j)- 30*input(i,j-1)+3*input(i,j-2))/(60*dx);   

    }else if(dir == "y_dir" && stencil == "left"){
      phi_output = (-2*input(i,j-3) + 15*input(i,j-2) - 60*input(i,j-1) + 20*input(i,j)+ 30*input(i,j+1)-3*input(i,j+2))/(60*dx); 
   
    }else if(dir == "x_dir" && stencil == "right"){
      
      phi_output = (2*input(i+3,j) - 15*input(i+2,j) + 60*input(i+1,j) - 20*input(i,j)- 30*input(i-1,j)+3*input(i-2,j))/(60*dx);

    }else if(dir == "x_dir" && stencil == "left"){
    
      phi_output = (-2*input(i-3,j) + 15*input(i-2,j) - 60*input(i-1,j) + 20*input(i,j)+ 30*input(i+1,j)-3*input(i+2,j))/(60*dx);
    } 
    return phi_output;
  }

 
}//End of NAMESPACE MLS

