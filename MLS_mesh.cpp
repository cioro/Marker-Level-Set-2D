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

    if(iter_counter%2 == 0){
      //Calculate phi^(n+1)
      //-----------------------------------------------------------------

      blitz::Array<double,2> phi_n_one_star;
      phi_n_one_star.resize(2*nGhost+ncells,2*nGhost+ncells);

      V_phi_n = spatial_first_X(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_first_Y(phi_n_one_star);

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
    
      V_phi_n_one = spatial_first_X(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_first_Y(phi_n_two_star);

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
    
      V_phi_n_half = spatial_first_X(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_first_Y(phi_n_three_half_star);

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

      V_phi_n = spatial_first_Y(phi_n);
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_one_star(row,col) = phi_n(row,col) - dt*V_phi_n(row,col);
	}
      }

      V_phi_n = spatial_first_X(phi_n_one_star);

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
    
      V_phi_n_one = spatial_first_Y(phi_n_one);

      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_two_star(row,col) = phi_n_one(row,col) - dt*V_phi_n_one(row,col);
	}
      }

      V_phi_n_one = spatial_first_X(phi_n_two_star);

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
    
      V_phi_n_half = spatial_first_Y(phi_n_half); 
    
      for(int row = nGhost; row < nGhost+ncells; ++row){
	for(int col = nGhost; col < nGhost+ncells; ++col){
	  phi_n_three_half_star(row,col) = phi_n_half(row,col) - dt*V_phi_n_half(row,col);
	}
      }

      V_phi_n_half = spatial_first_X(phi_n_three_half_star);

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

  void Mesh::correction1(){
    //std::cout <<"Inside correction function" << std::endl;
    std::vector<MLS::Particle> interface_nodes;
    Particle node_particle;
    int x_cell_index = 0;
    int y_cell_index = 0;
    double current_phi;
    double new_phi;
    double x_dist = 0;
    double y_dist = 0;
    double distance = 0;
    double x_quad,y_quad;
    //Update Distance
    //for each marker mk
    for(auto &marker : MLS_markers){
      
      //-------------------UPDATE Distance PROCEDURE--------------------------------
      
      //for each node phi(i,j) near mk-assume it is within the same cell.
      //In the future- could use quadrants to update the closest neigbouring cells.
      x_cell_index = int(floor((marker.x_coord-x_min)/dx)) + nGhost;
      y_cell_index = int(floor((marker.y_coord-y_min)/dy)) + nGhost;
      
      x_dist = x_cell_axis(x_cell_index);
      y_dist = y_cell_axis(y_cell_index);
     
      distance = sqrt((x_dist- marker.x_coord)*(x_dist - marker.x_coord) + (y_dist - marker.y_coord)*(y_dist - marker.y_coord));
      current_phi = MLS_data(x_cell_index,y_cell_index).phi;
      
      //abs(phi(i,j))=min(abs(phi(i,j),dist(node,mk)  

      if( MLS_data(x_cell_index,y_cell_index).phi < 0){
	new_phi =-1* double(std::min(distance,fabs(current_phi)));    
	MLS_data(x_cell_index,y_cell_index).phi = new_phi;
      } else if(MLS_data(x_cell_index,y_cell_index).phi > 0){
	new_phi = double(std::min(distance,fabs(current_phi)));    
	MLS_data(x_cell_index,y_cell_index).phi = new_phi;
      }
      //-------------------------------------------------------------------------
      
      //----------FINDING NEIGHBOURING CELLS-------------------------------------
      
      //std::cout << "The current cell has index : i " <<x_cell_index << " and j " << y_cell_index << std::endl;
      //std::cout << "This corresponds to node at x : " << x_dist << " and y : " << y_dist << std::endl; 
      x_quad = marker.x_coord-x_dist;
      y_quad = marker.y_coord-y_dist;

      if(x_quad >= 0 && y_quad >= 0){
	//std::cout << "marker is in first quadrant : " << x_quad << "\t" << y_quad << std::endl;
	//std::cout << " The current cel has index : i " << x_cell_index+1 << std::endl;
	//std::cout << "Are these two the same : x_cell_index " << x_cell_axis(x_cell_index+1) << "x_dist+dx" <<x_dist+dx << std::endl; 
	//Find distance to three uper corner nodes.
	//right
	distance = sqrt((x_dist+dx- marker.x_coord)*(x_dist+dx - marker.x_coord) + (y_dist - marker.y_coord)*(y_dist - marker.y_coord));
	if(MLS_data((x_cell_index+1),y_cell_index).phi < 0){
	MLS_data((x_cell_index+1),y_cell_index).phi = -1*std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index).phi));
	}else if(MLS_data((x_cell_index+1),y_cell_index).phi > 0){
	MLS_data((x_cell_index+1),y_cell_index).phi = std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index).phi));
	}
	
	//corner
	distance = sqrt((x_dist+dx- marker.x_coord)*(x_dist+dx - marker.x_coord) + (y_dist+dy - marker.y_coord)*(y_dist+dy - marker.y_coord));
	if(MLS_data((x_cell_index+1),y_cell_index+1).phi < 0 ){
	MLS_data((x_cell_index+1),y_cell_index+1).phi = -1*std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index+1).phi));
	}else if(MLS_data((x_cell_index+1),y_cell_index+1).phi > 0){
	MLS_data((x_cell_index+1),y_cell_index+1).phi = std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index+1).phi));
	}

	//upper
	
	distance = sqrt((x_dist- marker.x_coord)*(x_dist - marker.x_coord) + (y_dist+dy - marker.y_coord)*(y_dist+dy - marker.y_coord));
	if( MLS_data(x_cell_index,y_cell_index+1).phi < 0 ){
	  MLS_data(x_cell_index,y_cell_index+1).phi =-1* std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index+1).phi));
	}else if( MLS_data(x_cell_index,y_cell_index+1).phi > 0){
	  MLS_data(x_cell_index,y_cell_index+1).phi = std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index+1).phi));
	}
	
      }else if(x_quad < 0 && y_quad > 0 ){
	//std::cout << "marker is in second quadrant : " << x_quad << "\t" << y_quad << std::endl;
	//left
	distance = sqrt((x_dist-dx- marker.x_coord)*(x_dist-dx - marker.x_coord) + (y_dist - marker.y_coord)*(y_dist - marker.y_coord));
	if( MLS_data((x_cell_index-1),y_cell_index).phi < 0){
	  MLS_data((x_cell_index-1),y_cell_index).phi =-1*std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index).phi));
	}else if( MLS_data((x_cell_index-1),y_cell_index).phi > 0){
	  MLS_data((x_cell_index-1),y_cell_index).phi = std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index).phi));
	}

	//corner
	distance = sqrt((x_dist-dx- marker.x_coord)*(x_dist-dx - marker.x_coord) + (y_dist+dy - marker.y_coord)*(y_dist+dy - marker.y_coord));
	if( MLS_data(x_cell_index-1,y_cell_index+1).phi < 0){
	  MLS_data(x_cell_index-1,y_cell_index+1).phi =-1*std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index+1).phi));
	}else if( MLS_data(x_cell_index-1,y_cell_index+1).phi > 0){
	  MLS_data(x_cell_index-1,y_cell_index+1).phi = std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index+1).phi));
	}
	
	//upper
	distance = sqrt((x_dist- marker.x_coord)*(x_dist - marker.x_coord) + (y_dist+dy - marker.y_coord)*(y_dist+dy - marker.y_coord));
	if( MLS_data(x_cell_index,y_cell_index+1).phi < 0){
	  MLS_data(x_cell_index,y_cell_index+1).phi = -1*std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index+1).phi));
	}else if( MLS_data(x_cell_index,y_cell_index+1).phi > 0){
	  MLS_data(x_cell_index,y_cell_index+1).phi = std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index+1).phi));
	}
	
	
      }else if(x_quad < 0 && y_quad < 0){
	//left
	distance = sqrt((x_dist-dx- marker.x_coord)*(x_dist-dx - marker.x_coord) + (y_dist - marker.y_coord)*(y_dist - marker.y_coord));
	if( MLS_data((x_cell_index-1),y_cell_index).phi < 0){
	  MLS_data((x_cell_index-1),y_cell_index).phi =-1*std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index).phi));
	}else if( MLS_data((x_cell_index-1),y_cell_index).phi > 0){
	  MLS_data((x_cell_index-1),y_cell_index).phi = std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index).phi));
	}
	

	//lowerleft corner
	distance = sqrt((x_dist-dx- marker.x_coord)*(x_dist-dx - marker.x_coord) + (y_dist-dy - marker.y_coord)*(y_dist-dy - marker.y_coord));
	if(MLS_data((x_cell_index-1),y_cell_index-1).phi < 0){
	  MLS_data((x_cell_index-1),y_cell_index-1).phi =-1* std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index-1).phi));
	}else if(MLS_data((x_cell_index-1),y_cell_index-1).phi > 0){
	  MLS_data((x_cell_index-1),y_cell_index-1).phi = std::min(distance,fabs(MLS_data((x_cell_index-1),y_cell_index-1).phi));
	}

	//lower
	distance = sqrt((x_dist- marker.x_coord)*(x_dist - marker.x_coord) + (y_dist-dy - marker.y_coord)*(y_dist-dy - marker.y_coord));
	if( MLS_data(x_cell_index,y_cell_index-1).phi < 0){
	  MLS_data(x_cell_index,y_cell_index-1).phi = -1*std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index-1).phi));	
	}else if( MLS_data(x_cell_index,y_cell_index-1).phi > 0){
	  MLS_data(x_cell_index,y_cell_index-1).phi = std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index-1).phi));
	}
	

	//std::cout << "marker is in third  quadrant : " << x_quad << "\t" << y_quad << std::endl;
      }else if(x_quad > 0 && y_quad < 0){
	//std::cout << "marker is in fourth quadrant : " << x_quad << "\t" << y_quad << std::endl;
	//right
	distance = sqrt((x_dist+dx- marker.x_coord)*(x_dist+dx - marker.x_coord) + (y_dist - marker.y_coord)*(y_dist - marker.y_coord));
	if( MLS_data((x_cell_index+1),y_cell_index).phi < 0 ){
	  MLS_data((x_cell_index+1),y_cell_index).phi = -1*std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index).phi));
	}else if( MLS_data((x_cell_index+1),y_cell_index).phi > 0){
	  MLS_data((x_cell_index+1),y_cell_index).phi = std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index).phi));
	}
	
	//corner
	distance = sqrt((x_dist+dx- marker.x_coord)*(x_dist+dx - marker.x_coord) + (y_dist-dy - marker.y_coord)*(y_dist-dy - marker.y_coord));
	if(  MLS_data((x_cell_index+1),y_cell_index-1).phi < 0){
	MLS_data((x_cell_index+1),y_cell_index-1).phi = -1*std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index-1).phi));
	}else if(  MLS_data((x_cell_index+1),y_cell_index-1).phi > 0){
	  MLS_data((x_cell_index+1),y_cell_index-1).phi = std::min(distance,fabs(MLS_data((x_cell_index+1),y_cell_index-1).phi));
	}
	

	//lower
	distance = sqrt((x_dist- marker.x_coord)*(x_dist - marker.x_coord) + (y_dist-dy - marker.y_coord)*(y_dist-dy - marker.y_coord));
	if( MLS_data(x_cell_index,y_cell_index-1).phi < 0){
	  MLS_data(x_cell_index,y_cell_index-1).phi =-1*std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index-1).phi));
	}else if( MLS_data(x_cell_index,y_cell_index-1).phi > 0){
	  MLS_data(x_cell_index,y_cell_index-1).phi = std::min(distance,fabs(MLS_data((x_cell_index),y_cell_index-1).phi));
	}
	

      }
      
      
      //-------------------------------------------------------------------------
      
      //-------------CREATING VECTOR OF INTERFACE CELL NODES ---------------
      node_particle.x_coord = x_dist;
      node_particle.y_coord = y_dist;
      interface_nodes.push_back(node_particle);
      //--------------------------------------------------------------------
      
    }
    
    //---------------UPDATE SIGN OF INTERFACIAL NODES ----------------------------
    /*
    double d1,d2;
    double x_node_index,y_node_index;
    //Update Sign
    for(auto &node : interface_nodes){
      x_node_index =  int(floor((marker.x_coord-x_min)/dx)) + nGhost;
      y_node_index = ();
      d1;// = sqrt(node.x_coord-
    
      }*/
    //----------------------------------------------------------------------------

  }//End of correction fcn


  void Mesh::marchingSquares(){
    double zero_threshold = 0.0;
    double plus_delta_x_threshold = dx;
    double minus_delta_x_threshold = -1*dx;
    int case_zero;
    int case_plus_x;
    int case_minus_x;
    std::vector<std::pair <double,double>> zero_pair,plus_pair,minus_pair;
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

	interface.insert({key_index,zero_pair});
	minus_dx_LS.insert({key_index,minus_pair});
	plus_dx_LS.insert({key_index,plus_pair});

	}
    }
    //Check for each threshold-determined with case we are dealing with.
    //
    //Check the different cases
    //update the maps as required.
  }

  int Mesh::Case(double threshold, double vert1,double vert2,double vert3,double vert4, double threshold){
   
    int case_output = 0;
    
    if(vert1 > threshold){ case_output += 1 ;}
    if(vert2 > threshold){ case_output += 2 ;}
    if(vert3 > threshold){ case_output += 8 ;}  
    if(vert4 > threshold){ case_output += 4 ;}
  
   
    return case_output;
  }

  std::vector<std::pair<double,double>> Mesh::squareReconstruction(int Case, int i ,int j,double v1,double v2,double v3,double v4,double threshold){

    std::vector<std::pair<double,double>> output;
    //the first member of the pair is the gradient of the line, the second one a point on the line.
    std::pair<double,double> values;
    double a1,a2,a3,a4;

    if(Case == 0){
      output.resize(1);
      values.first = -1;
      values.second = -1;
      output.push_back(values);
    }else if(Case == 15){
      output.resize(1);
      values.first = -1;
      values.second = -1;
      output.push_back(values);
    }else if(Case == 1){
      a1 = (v1-threshold)/(v1-v2);
    }else if(Case == 2){
    
    }else if(Case == 3){
    }else if(Case == 4){
    }else if(Case == 5){
    }else if(Case == 6){    
    }else if(Case == 7){
    }else if(Case == 8){
    }else if(Case == 9){
    }else if(Case == 10){
    }else if(Case == 11){
    }else if(Case == 12){
    }else if(Case == 13){
    }else if(Case == 14){
    }
    return output;
    
  }

}//End of NAMESPACE MLS
