#include<vector>
#include"MLS_cell.hpp"
#include"MLS_mesh.hpp"
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>

namespace MLS{

  Mesh::Mesh(){};
  Mesh::Mesh(double AT_max,int arg_ncells, int AnGhost, double Ax_min, double Ax_max,double Ay_min,double Ay_max,double arg_cfl, double (*s_x)(double x, double y, double t, double T), double (*s_y)(double x, double y, double t, double T), Cell (*level_set)(double x, double y, double time, double T_max)) :T_max(AT_max), ncells(arg_ncells),nGhost(AnGhost), x_min(Ax_min),x_max(Ax_max),y_min(Ay_min),y_max(Ay_max),cfl(arg_cfl), iter_counter(0),time(0), speed_x(s_x), speed_y(s_y)
  {
     
    dx = (x_max-x_min)/(double)ncells;
    dy = (y_max-y_min)/(double)ncells;
 
    MLS_data.resize(ncells + 2*nGhost,ncells + 2*nGhost);
    xaxis.resize(ncells + 2*nGhost);
    yaxis.resize(ncells + 2*nGhost);
   
    int x_counter = 0;//Used to correctly calculate value of xaxis (value starts at 0 not nGhost,
    int y_counter = 0;//but first element is nGhost; Difference between array index and represented value
   

    //fill in xaxis
    for(int i = nGhost; i <(ncells+nGhost); i++){
      xaxis(i) = x_min + x_counter*dx;
      x_counter++;
     
      for(int j = nGhost; j<(ncells+nGhost); j++){
	yaxis(j) = y_min + y_counter*dy;
	y_counter++;

	MLS_data(i,j) = level_set(xaxis(i),yaxis(j),time,T_max);

      }
      y_counter=0;
     
    }
  
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
     
	speedtemp_x = this->speed_x(xaxis(col),yaxis(row),time,T_max);
	//std::cout << speedtemp_x << std::endl;
	speedtemp_y = this->speed_y(xaxis(col),yaxis(row),time,T_max);
      
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
    
      std::cout << "Inside calculate dt function, cfl = " << cfl_init << " dt " << dt << "\n"; 
    
    }else{
      std::cout << "Inside calculate dt function, cfl = " << cfl << " dt " << dt << "\n";
  
      dt = cfl*min_coef;
    }

  
  };
  
  void Mesh::advect_level_set(){

     std::cout << "ADVECTING" << "\n";

    double phi;
    double phi_xdir,phi_ydir;
    double phi_xdir_plus,phi_xdir_minus;
    double phi_ydir_plus,phi_ydir_minus;
    blitz::Array<double,2>phi_new;
    double obj_speed_x, obj_speed_y;
    
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_x = this->speed_x(xaxis(col),yaxis(row),time,T_max);
	obj_speed_y = this->speed_y(xaxis(col),yaxis(row),time,T_max);
	phi = MLS_data(col,row).phi;
	phi_xdir_plus = MLS_data(col+1,row).phi;
	phi_xdir_minus = MLS_data(col-1,row).phi;
	phi_ydir_plus = MLS_data(col,row+1).phi;
	phi_ydir_minus = MLS_data(col,row-1).phi;

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
	phi = phi - dt*(obj_speed_x*phi_xdir+obj_speed_y*phi_ydir);
	phi_new(col,row) = phi;
      }
    }

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	MLS_data(col,row).phi = phi_new(col,row);
      }
    }


};

  void Mesh::save_to_file(std::string filename)const{

    std::string dir = "data/";
    std::stringstream ss;
    ss << dir << filename << iter_counter;
    std::string tmppath = ss.str();
  
    std::cout << "CREATING FILE \n";
    FILE * outfile = fopen(tmppath.c_str(),"w");
  
    for(int i=nGhost; i<ncells+nGhost; i++)
      {
	for(int j=nGhost; j<ncells+nGhost; j++){

	  fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f  \n", xaxis(i), yaxis(j),MLS_data(i,j).phi,MLS_data(i,j).phi_u,MLS_data(i,j).phi_v);
	}

	fprintf(outfile,"\n");
      
      }
    fclose(outfile);

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

  double Mesh::D_minus(int i, int j, std::string dir){

    double phi;
    double phi_minus;
    double phi_result;
    if(dir == "x_dir"){
      phi = MLS_data(i,j).phi;
      phi_minus = MLS_data((i-1),j).phi;
      phi_result = (phi - phi_minus )/dx;
    }else if (dir == "y_dir"){
      phi = MLS_data(i,j).phi;
      phi_minus = MLS_data(i,(j-1)).phi;
      phi_result = (phi - phi_minus )/dy;
    }
    return phi_result;

  }

  double Mesh::D_plus(int i, int j, std::string dir){

    double phi;
    double phi_plus;
    double phi_result;
    if(dir == "x_dir"){
      phi = MLS_data(i,j).phi;
      phi_plus = MLS_data((i+1),j).phi;
      phi_result = (phi_plus - phi)/dx;
    }else if (dir == "y_dir"){
      phi = MLS_data(i,j).phi;
      phi_plus = MLS_data(i,(j+1)).phi;
      phi_result = (phi_plus - phi)/dy;
    }
    return phi_result;

  }

  double Mesh::WENO(int i, int j, std::string stencil, std::string dir){
    
    double S1,S2,S3;
    double a1,a2,a3;
    double w1,w2,w3;
    double v1,v2,v3,v4,v5;
    double phi1,phi2,phi3;
    double phi_result;
    double eps; //epsilon
   
    if(stencil == "minus" && dir == "x_dir"){ 
      v1 = D_minus((i-2),j,dir);
      v2 = D_minus((i-1),j,dir);
      v3 = D_minus((i),j,dir);
      v4 = D_minus((i+1),j,dir);
      v5 = D_minus((i+2),j,dir);
    }else if(stencil == "plus" && dir == "x_dir"){
      v1 = D_plus((i+2),j,dir);
      v2 = D_plus((i+1),j,dir);
      v3 = D_plus((i),j,dir);
      v4 = D_plus((i-1),j,dir);
      v5 = D_plus((i-2),j,dir);
    }else if(stencil == "minus" && dir == "y_dir"){ 
      v1 = D_minus(i,(j-2),dir);
      v2 = D_minus(i,(j-1),dir);
      v3 = D_minus(i,j,dir);
      v4 = D_minus(i,(j+1),dir);
      v5 = D_minus(i,(j+2),dir);
    }else if(stencil == "plus" && dir == "y_dir"){
      v1 = D_plus(i,(j+2),dir);
      v2 = D_plus(i,(j+1),dir);
      v3 = D_plus(i,j,dir);
      v4 = D_plus(i,(j-1),dir);
      v5 = D_plus(i,(j-2),dir);
    }else{
      std::cout << "No upwinding differencing stencil selected" << "\n";
    }
    
    phi1 =  (1.0/3.0)*v1 - (7.0/6.0)*v2 + (11.0/6.0)*v3;
    phi2 = (-1.0/6.0)*v2 + (5.0/6.0)*v3 + (1.0/3.0)*v4;
    phi3 =  (1.0/3.0)*v3 + (5.0/6.0)*v4 - (1.0/6.0)*v5;
 
    S1 = (13.0/12.0)*(v1-2*v2+v3)*(v1-2*v2+v3) + \
      (1.0/4.0)*(v1-4*v2+3*v3)*(v1-4*v2+3*v3);

    S2 = (13.0/12.0)*(v2-2*v3+v4)*(v2-2*v3+v4) + \
      (1.0/4.0)*(v2-v4)*(v2-v4);

    S3 = (13.0/12.0)*(v3-2*v4+v5)*(v3-2*v4+v5) + \
      (1.0/4.0)*(3*v3-4*v4+v5)*(3*v3-4*v4+v5);
      
    std::vector<double> v = {v1,v2,v3,v4,v5};
    eps = 1e-6*(*std::max_element(v.begin(),v.end()))+1e-99;

    a1 = 0.1/((S1+eps)*(S1+eps));
    a2 = 0.6/((S2+eps)*(S2+eps));
    a3 = 0.3/((S3+eps)*(S3+eps));

    double sum_a = a1+a2+a3;
      
    w1 = a1/sum_a;
    w2 = a2/sum_a;
    w3 = a3/sum_a;

    phi_result = w1*phi1 + w2*phi2 + w3*phi3;

    return phi_result;
 
  }



  void Mesh::advect_WENO(){

     std::cout << "ADVECTING" << "\n";

    double phi;
    double phi_xdir,phi_ydir;
    double phi_xdir_plus,phi_xdir_minus;
    double phi_ydir_plus,phi_ydir_minus;
    blitz::Array<double,2>phi_new;
    double obj_speed_x, obj_speed_y;
    std::string stencil;
    std::string dir;
    phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	obj_speed_x = this->speed_x(xaxis(col),yaxis(row),time,T_max);
	obj_speed_y = this->speed_y(xaxis(col),yaxis(row),time,T_max);
	phi = MLS_data(col,row).phi;
	phi_xdir_plus = MLS_data(col+1,row).phi;
	phi_xdir_minus = MLS_data(col-1,row).phi;
	phi_ydir_plus = MLS_data(col,row+1).phi;
	phi_ydir_minus = MLS_data(col,row-1).phi;

	if(obj_speed_x < 0){
	  stencil = "plus";
	  dir = "x_dir";
	  phi_xdir = WENO(row,col,stencil,dir); 
	}else if(obj_speed_x == 0){
	   phi_xdir = 0.0;
	}else if(obj_speed_x > 0){
	  stencil = "minus";
	  dir = "x_dir";
	  phi_xdir = WENO(row,col,stencil,dir); 
	}else{
	  std::cout << "Advection of level_set. Something went wrong" <<"\n";
	}
      
	if(obj_speed_y < 0){
	  stencil = "plus";
	  dir = "y_dir";
	  phi_ydir = WENO(row,col,stencil,dir);
	}else if(obj_speed_y == 0){
	  phi_ydir = 0;
	}else if(obj_speed_y > 0){
	  stencil = "minus";
	  dir = "y_dir";
	  phi_ydir = WENO(row,col,stencil,dir);
	}else{
	  std::cout << "Advection of Level_set" <<"\n";
	}
	// phi_ydir = 0;
	phi = phi - dt*(obj_speed_x*phi_xdir+obj_speed_y*phi_ydir);
	phi_new(col,row) = phi;
      }
    }

    for(int row = nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
	MLS_data(col,row).phi = phi_new(col,row);
      }
    }


  }
}
