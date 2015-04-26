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
  Mesh::Mesh(int arg_ncells, int AnGhost, double Ax_min, double Ax_max,double Ay_min,double Ay_max,double arg_cfl, Cell (*s_x)(double x, double y, double t), Cell (*s_y)(double x, double y, double t), Cell (*level_set)(double x, double y, double Atime)) : ncells(arg_ncells),nGhost(AnGhost), x_min(Ax_min),x_max(Ax_max),y_min(Ay_min),y_max(Ay_max),cfl(arg_cfl), iter_counter(0),time(0), speed_x(s_x), speed_y(s_y)
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

	MLS_data(i,j) = level_set(xaxis(i),yaxis(j),time);

      }
      y_counter=0;
     
    }
  
  };

  Mesh::~Mesh(){};

  void Mesh::Calculate_dt(){

    double speed_x=0.0;
    double speedtemp_x=0.0;
    double speed_y=0.0;
    double speedtemp_y=0.0;
    double min_coef = 0.0;//This is min of (S_x/dx,S_y/dy)-The minimum of the wave speed over the space step.
 
    for(int row=nGhost; row < nGhost+ncells; row++){
      for(int col = nGhost; col < nGhost+ncells; col++){
     
	speedtemp_x = MLS_data(row,col).phi_u;
      
	speedtemp_y = MLS_data(row,col).phi_v;
      
	if(fabs(speedtemp_x) > fabs(speed_x)){
	  speed_x = speedtemp_x;
	}
	if(fabs(speedtemp_y) > fabs(speed_y)){
	  speed_y = speedtemp_y;
	}          
      }
    }
  
    min_coef = std::min(fabs(dx/speed_x),fabs(dy/speed_y));
    double dt;
    //If time < 5 then dt = 0.2
  
    if(iter_counter < 10){
      double  cfl_init = 0.2;
      dt = cfl_init*min_coef;
    
      std::cout << "Inside calculate dt function, cfl = " << cfl_init << "\n"; 
    
    }else{
      std::cout << "Inside calculate dt function, cfl = " << cfl << "\n";
  
      dt = cfl*min_coef;
    }

  
  };
  
  void Mesh::advect_level_set(){};

  void Mesh::save_to_file(std::string filename)const{

    std::string dir = "data/";
    std::stringstream ss;
    ss << dir << filename << time;
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

};
