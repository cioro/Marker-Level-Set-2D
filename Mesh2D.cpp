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

//Initialises the parameters of the grid and fills the vector of primitives with the initial conditions
Mesh::Mesh(){};
Mesh::Mesh(int Ancells, double Ax_min, double Ax_max,double Ay_min,double Ay_max,double Acfl, Euler::U_state (*f)(double x,double y), Euler::W_state (*b1)(Euler::W_state  w), Euler::W_state (*b2)(Euler::W_state w), int AnGhost, std::string arg_BC, MLS (*level_set)(double x, double y)) : ncells(Ancells), x_min(Ax_min),x_max(Ax_max),y_min(Ay_min),y_max(Ay_max),cfl(Acfl), time(0), boundary1(b1), boundary2(b2), nGhost(AnGhost), BC(arg_BC)
 { 
   //Set boundary conditions:
   
   ptr_euler = new Euler();
   dx = (x_max-x_min)/(double)ncells;
   dy = (y_max-y_min)/(double)ncells;

   Bdata.resize(ncells + 2*nGhost,ncells + 2*nGhost);
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
       Bdata(i,j) = f(xaxis(i),yaxis(j));
       MLS_data(i,j) = level_set(xaxis(i),yaxis(j));
     }
     y_counter=0;
     
   }
  
 }

//Destructor
Mesh::~Mesh()
{
  delete ptr_euler;
}

void Mesh::reset(Euler::U_state (*f)(double x,double y)){

  int x_counter = 0;//Used to correctly calculate value of xaxis (value starts at 0 not nGhost,
  int y_counter = 0;//but first element is nGhost; Difference between array index and represented value
    
  for(int i = nGhost; i <(ncells+nGhost); i++){
    xaxis(i) = x_min + x_counter*dx;//fill in x-axis
    x_counter++;
    
    for(int j = nGhost; j<(ncells+nGhost); j++){
      yaxis(j) = y_min + y_counter*dy;//fill in y-axis
      y_counter++;
      Bdata(i,j) = f(xaxis(i),yaxis(j));
     
    }

    y_counter=0;
    
  } 

};

void Mesh::reset(Euler::U_state (*f)(double x,double y), MLS (*level_set)(double x, double y)){

  int x_counter = 0;//Used to correctly calculate value of xaxis (value starts at 0 not nGhost,
  int y_counter = 0;//but first element is nGhost; Difference between array index and represented value
    
  for(int i = nGhost; i <(ncells+nGhost); i++){
    xaxis(i) = x_min + x_counter*dx;//fill in x-axis
    x_counter++;
    
    for(int j = nGhost; j<(ncells+nGhost); j++){
      yaxis(j) = y_min + y_counter*dy;//fill in y-axis
      y_counter++;
      if(MLS_data(i,j).phi <= 0){
      Bdata(i,j) = f(xaxis(i),yaxis(j));
      MLS_data(i,j) = level_set(xaxis(i),yaxis(j));
      }else if(MLS_data(i,j).phi > 0){
	Bdata(i,j).rho = 1.0;
	Bdata(i,j).moment_u = 0.0;
	Bdata(i,j).moment_v = 0.0;
	Bdata(i,j).energy = 1.0;
      }      
    }

    y_counter=0;
    
  } 

};


void Mesh::reset_BC(Euler::W_state (*b1)(Euler::W_state w),Euler::W_state (*b2)(Euler::W_state w),std::string arg_BC){
      BC = arg_BC;
      boundary1=b1;
      boundary2=b2;
      
    }
 

//Prints vector of conserved variables to screen
void Mesh::print()const
{
  for (int i = nGhost; i<(ncells+nGhost); i++){
    for (int j = nGhost; j<(ncells+nGhost); j++){
      std::cout << xaxis(i) << "\t" << yaxis(j) << "\t" << Bdata(i,j).rho <<"\t"<< Bdata(i,j).moment_u <<"\t" << Bdata(i,j).moment_v <<"\t" << Bdata(i,j).energy << "\n";
    }
    std::cout <<"\n";
  }
    
}

//Print to a file the 2D matrix of conserved variables and the exact solution
void Mesh::save_u_state(std::string filename)const
{
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  
  std::cout << "CREATING FILE \n";
  FILE * outfile = fopen(tmppath.c_str(),"w");

  for(int i=nGhost; i<ncells+nGhost; i++)
    {
      for(int j=nGhost; j<ncells+nGhost; j++){
	fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n", xaxis(i), yaxis(j),Bdata(i,j).rho, Bdata(i,j).moment_u,Bdata(i,j).moment_v,Bdata(i,j).energy, MLS_data(i,j).phi);
      }
      fprintf(outfile,"\n");
  }
  fclose(outfile);
  
}

//Prints to a file the 2D matrix of primitive variables  
void Mesh::save_w_state(std::string filename)const
{
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  
  std::cout << "CREATING FILE \n";
  FILE * outfile = fopen(tmppath.c_str(),"w");
  double vel_mag;
  Euler::W_state w_print;
  
  for(int i=nGhost; i<ncells+nGhost; i++)
    {
      for(int j=nGhost; j<ncells+nGhost; j++){
	w_print = ptr_euler->PfromC(Bdata(i,j));
	vel_mag = std::sqrt(w_print.v*w_print.v +w_print.u*w_print.u);
	fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f  \t %.4f  \t %.4f \n", xaxis(i), yaxis(j),w_print.rho,w_print.u,w_print.v,w_print.P,MLS_data(i,j).phi,vel_mag);
      }
      fprintf(outfile,"\n");
      
    }
      fclose(outfile);
  
}

//1D slice of u_states parallel to the x-axis.
void Mesh::slice_x_axis(std::string filename)const
{
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  
  std::cout << "CREATING FILE \n";
  FILE * outfile = fopen(tmppath.c_str(),"w");
  Euler::W_state w_print;  

  int y_0 = floor(ncells/2.0);
  std::cout << "slice at : " << y_0 << "\n";
  
  for(int x_axis = nGhost; x_axis < ncells+nGhost; x_axis++){

    w_print = ptr_euler->PfromC(Bdata(x_axis,y_0));
    
    fprintf(outfile, "%.6f \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f \n", xaxis(x_axis), w_print.rho,w_print.u,w_print.v,w_print.P, ptr_euler->int_energy(w_print));
        
  }
  fclose(outfile);
  
}

//1D slice of u_states parallel to the y-axis.
void Mesh::slice_y_axis(std::string filename)const
{
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  
  std::cout << "CREATING FILE \n";
  FILE * outfile = fopen(tmppath.c_str(),"w");
  Euler::W_state w_print;  
  int x_0 = int(ncells/2);
  
  for(int y_axis = nGhost; y_axis < ncells+nGhost; y_axis++){
	
    w_print = ptr_euler->PfromC(Bdata(x_0,y_axis));
    
    fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n", yaxis(y_axis), w_print.rho, w_print.u,w_print.v, w_print.P, ptr_euler->int_energy(w_print));
        
  }
  fclose(outfile);
  
}


//Implements the boundary conditions. The actual boundary condition function should be in the main file
void Mesh::applyBC(){

  if(BC == std::string("Transmissive")){
      //Left Boundary of Square Mesh
      // Loop over Ghost left columns
      for(int x_dir = 0; x_dir < nGhost; x_dir++){

	//Loop over each row
	for(int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){
	  /*
	    Euler::W_state w_Left_End = ptr_euler -> PfromC(Bdata(x_dir+nGhost,y_dir));
	    Euler::W_state w_BC_Left = boundary1(w_Left_End);
	    Bdata(nGhost-1-x_dir,y_dir) = ptr_euler-> CfromP(w_BC_Left);
	  */
     
	  Bdata(x_dir,y_dir)= Bdata(nGhost,y_dir);
      
	}
      }

      //Right Boundary of Square Mesh
      //Loop over Ghost right columns
      for(int x_dir = 0; x_dir < nGhost; x_dir++ ){
	//Loop over each row
	for (int y_dir = nGhost; y_dir < ncells+nGhost; y_dir++){
	  /*
	    Euler::W_state w_Right_End = ptr_euler -> PfromC(Bdata ((nGhost+ncells-1)-x_dir,y_dir));
	    Euler::W_state w_BC_Right = boundary2(w_Right_End);
	    Bdata(x_dir+(nGhost+ncells),y_dir)= ptr_euler -> CfromP(w_BC_Right);
	  */
     
	  Bdata(ncells+nGhost+x_dir,y_dir)= Bdata(ncells+nGhost-1,y_dir);

	}
      }

      //Top Boundary of Square Mesh
      //Loop over Ghost top rows
      for(int y_dir = 0; y_dir < nGhost; y_dir++){

	//Loop over each column
	for(int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	  /*
	    Euler::W_state w_Left_End = ptr_euler -> PfromC(Bdata(x_dir, (nGhost+ncells-1)-y_dir));
	    Euler::W_state w_BC_Left = boundary1(w_Left_End);
	    Bdata(x_dir,y_dir+(nGhost+ncells)) = ptr_euler-> CfromP(w_BC_Left);
	  */

	  Bdata(x_dir,y_dir+ncells+nGhost) = Bdata(x_dir,ncells+nGhost-1);

	}
      }

      //Bottom Boundary of Square Mesh
      //Loop over Ghost bottom rows
      for(int y_dir = 0; y_dir < nGhost; y_dir++ ){
	//Loop over each column
	for (int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	  /*
	    Euler::W_state w_Right_End = ptr_euler -> PfromC(Bdata(x_dir, nGhost + y_dir));
	    Euler::W_state w_BC_Right = boundary2(w_Right_End);
	    Bdata(x_dir, nGhost-1-y_dir)= ptr_euler -> CfromP(w_BC_Right);
	  */

	  Bdata(x_dir,y_dir) = Bdata(x_dir,nGhost);

	}
      }
  }else if(BC==std::string("Reflective")){

	//Left Boundary of Square Mesh
	// Loop over Ghost left columns
	for(int x_dir = 0; x_dir < nGhost; x_dir++){

	  //Loop over each row
	  for(int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){
	    Euler::W_state w_Left_End = ptr_euler -> PfromC(Bdata(x_dir+nGhost,y_dir));
	    Euler::W_state w_BC_Left = boundary1(w_Left_End);
	    Bdata(nGhost-1-x_dir,y_dir) = ptr_euler-> CfromP(w_BC_Left);
	    
	  }
	}
	
	//Right Boundary of Square Mesh
	//Loop over Ghost right columns
	for(int x_dir = 0; x_dir < nGhost; x_dir++ ){
	  //Loop over each row
	  for (int y_dir = nGhost; y_dir < ncells+nGhost; y_dir++){
	    Euler::W_state w_Right_End = ptr_euler -> PfromC(Bdata ((nGhost+ncells-1)-x_dir,y_dir));
	    Euler::W_state w_BC_Right = boundary2(w_Right_End);
	    Bdata(x_dir+(nGhost+ncells),y_dir)= ptr_euler -> CfromP(w_BC_Right);
	    
	  }
	}

	Euler::W_state w_permute_temp;
	//Top Boundary of Square Mesh
	//Loop over Ghost top rows
	for(int y_dir = 0; y_dir < nGhost; y_dir++){

	  //Loop over each column
	  for(int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	    Euler::W_state w_Left_End = ptr_euler -> PfromC(Bdata(x_dir, (nGhost+ncells-1)-y_dir));
	    //Permute before applying BC.
	    w_permute_temp.u = w_Left_End.u;
	    w_permute_temp.v = w_Left_End.v;
	    w_Left_End.u = w_permute_temp.v;
	    w_Left_End.v = w_permute_temp.u;
	    Euler::W_state w_BC_Left = boundary1(w_Left_End);
	    //Permute after BC has been applied.
	    w_permute_temp.u = w_BC_Left.u;
	    w_permute_temp.v = w_BC_Left.v;
	    w_BC_Left.u = w_permute_temp.v;
	    w_BC_Left.v = w_permute_temp.u;

	    Bdata(x_dir,y_dir+(nGhost+ncells)) = ptr_euler-> CfromP(w_BC_Left);
	  }
	}

	//Bottom Boundary of Square Mesh
	//Loop over Ghost bottom rows
	for(int y_dir = 0; y_dir < nGhost; y_dir++ ){
	  //Loop over each column
	  for (int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	    Euler::W_state w_Right_End = ptr_euler -> PfromC(Bdata(x_dir, nGhost + y_dir));
	    w_permute_temp.u = w_Right_End.u;
	    w_permute_temp.v = w_Right_End.v;
	    w_Right_End.u = w_permute_temp.v;
	    w_Right_End.v = w_permute_temp.u;
	    Euler::W_state w_BC_Right = boundary2(w_Right_End);
	    w_permute_temp.u = w_BC_Right.u;
	    w_permute_temp.v = w_BC_Right.v;
	    w_BC_Right.u = w_permute_temp.v;
	    w_BC_Right.v = w_permute_temp.u;
	    
	    Bdata(x_dir, nGhost-1-y_dir)= ptr_euler -> CfromP(w_BC_Right);
	    
	  }
	}



      }



}

//Calculates the adquate size of the time step dt. See page 183 from Toro(ed.2009)
double Mesh::Calculate_dt(){

  double speed_x=0.0;
  double speedtemp_x=0.0;
  double speed_y=0.0;
  double speedtemp_y=0.0;
  double min_coef = 0.0;//This is min of (S_x/dx,S_y/dy)-The minimum of the wave speed over the space step.
  for(int row=nGhost; row < nGhost+ncells; row++){
    for(int col = nGhost; col < nGhost+ncells; col++){
      if(MLS_data(col,row).phi < 0){
	Euler::W_state w = ptr_euler->PfromC(Bdata(col,row));
	if(w.P < 0 || w.rho < 0){
   
	  std::cout << "The cell whose speed of sound cannot be computed is : " << col << "\t" << row << "\n";
	}
	speedtemp_x = ptr_euler->a(w) + fabs(w.u);
      
	speedtemp_y = ptr_euler->a(w) + fabs(w.v);
    
	if(fabs(speedtemp_x) > fabs(speed_x)){
	  speed_x = speedtemp_x;
	}
	if(fabs(speedtemp_y) > fabs(speed_y)){
	  speed_y = speedtemp_y;
	}          
      }
    }

  }

  min_coef = std::min(fabs(dx/speed_x),fabs(dy/speed_y));
  double dt;
  //If time < 5 then dt = 0.2
  
 if(time < 10){
    double  cfl_init = 0.2;
    dt = cfl_init*min_coef;
    
    std::cout << "Inside calculate dt function, cfl = " << cfl_init << "\n"; 

    return dt;

    }
  std::cout << "Inside calculate dt function, cfl = " << cfl << "\n";

  dt = cfl*min_coef;
  return dt;
}

//Need to modify this 1D  WAF.Unsure whether to create WAF_x_dir and WAF_y_dir or just permute v and u and use one fcn. 
//This function is not strictly applicable to the 1D case. It includes the existance of a tangential velocity v, which gives
//rise to a sheer wave, thus a new wave is included in the WAF calculation (plus a new limiter -- see page 553 Toro ed.2009)

blitz::Array<Euler::U_state,1> WAF_1D(blitz::Array<Euler::U_state,1> input_data, double dt, double ds, double ncells,double nGhost,std::string limiter,std::string sweep){
 
  //Total vector of fluxes
  blitz::Array<Euler::U_state,1> flux(ncells+1);
  blitz::Array<Euler::U_state,1> flux_temp(ncells+1);
  blitz::Array<Euler::U_state,1> input_temp(ncells+2*nGhost);
  //Logic to switch u and v depending on the sweep
  //Need to rever this inversion when the result is calculated at the end. 
  //Otherwise, I'll be calculting only the x-sweep twice.
  if (sweep == std::string("x-sweep")){

  }else if(sweep == std::string("y-sweep")){
    for(int i= 0; i < ncells+2*nGhost; i++){
      input_temp(i).moment_u = input_data(i).moment_u;
      input_temp(i).moment_v = input_data(i).moment_v;

      input_data(i).moment_u = input_temp(i).moment_v;
      input_data(i).moment_v = input_temp(i).moment_u;
    }
  }else{
    std::cout <<"Sweep not specified. Fail to compute 1D WAF :" << sweep << "\n";
  }

  //Given that the mesh is no longer an input argument to access Euler class function, an empty euler object is needed
  // to access the different data members and function members
  Euler e;
  double gamma = e.gamma;
  //Variables used

  double P_star,P_L,P_R; //Presures
  double rho_star,rho_L,rho_R;//densities
  double u_L,u_R;//particle/gas speed in cell

  double a_L, a_R;// sound speed in cell
  double S_L, S_R,S_star;// wave speed in cell left wave, right wave, contact wave

  Euler::U_state U_star_L, U_star_R; // Star states of conserved var
  
  double P_pvrs;
  double rho_bar, a_bar;//average density and average sound speed of left and right cells

  Euler::W_state w_temp_left;
  Euler::W_state w_temp_right;

  double q_R,q_L;

  Euler::U_state U_state_L;
  Euler::U_state U_state_R;
  Euler::W_state W_L;
  Euler::W_state W_R;
  Euler::U_state U_state_L_star;
  Euler::U_state U_state_R_star;

  double star_coef_left; // The coeficient in eq. 10.73 from Toro(ed.2009); (S_k-u_k)/(S_k-u_star_k);
  double star_coef_right; // The coeficient in eq. 10.73 from Toro(ed.2009);
  
  //Courant number for each wave speed
  double c_L,c_star,c_R,c_shear;

  //Ratio for each wave speed
  double r_L,r_star,r_R,r_shear;
  double dq_l,dq_l_right_interface, dq_l_left_interface;
  double dq_star, dq_star_right_interface, dq_star_left_interface;
  double dq_r, dq_r_right_interface, dq_r_left_interface;
  double dq_shear, dq_shear_right_interface, dq_shear_left_interface;

  //Limiter functions
  double minmod_l,minmod_star,minmod_r,minmod_shear;
  double superbee_l,superbee_star,superbee_r,superbee_shear;

  blitz::Array<Euler::U_state,1> left_interface;
  blitz::Array<Euler::U_state,1> right_interface;
  left_interface.resize(4);
  right_interface.resize(4);

  //Dummy variable used to calculate shear wave max speed
  double shear_speed;
  
  //Middle state used in WAF calculation to account for the presence of a shear wave.
  Euler::W_state w_middle_state_left;
  Euler::W_state w_middle_state_right;
  Euler::U_state F_ML;
  Euler::U_state F_MR;
  int f_idx = 0;
  //std::cout <<"Inside WAF_1D flux function " << "\n";
  //Loop over whole domain
  for(int i = nGhost-1; i < ncells + nGhost; i++){

    // std::cout << "Inside the WAF_1D loop. Iteration: " << i << "\n";
   
    U_state_L = input_data(i);
    U_state_R = input_data(i+1);
   
    W_L = e.PfromC(U_state_L);
    W_R = e.PfromC(U_state_R);
    
    //---------Pressure estimate-------------------------

    P_L = W_L.P;
    P_R = W_R.P;

    u_L = W_L.u;
    u_R = W_R.u;

    rho_L = W_L.rho;
    rho_R = W_R.rho;
   
    // std::cout << "Inside WAF 1D" << "\n";
    a_L = e.a(W_L);
    a_R = e.a(W_R);

    // std::cout <<"Inside the HLLC function" << "\n";
        
    rho_bar = 0.5*(rho_L + rho_R);
    a_bar = 0.5*(a_L + a_R);

    P_pvrs = 0.5*(P_L + P_R)-0.5*(u_R-u_L)*rho_bar*a_bar;

    P_star = std::max(0.0,P_pvrs);
      
    //----------Wave speed estimate----------------------
    
    //Calculate q_R

    if(P_star <= P_R){
      q_R = 1.0;
    }
    else{
      q_R = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_R)-1));
    }
        
    //Calculate q_L
    if(P_star <= P_L){
      q_L = 1.0;
    }
    else{
      q_L = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_L)-1));
    }
 
    //Calculate S_R and S_L
    S_L = u_L -a_L*q_L;
    S_R = u_R + a_R*q_R;
    
    double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
    double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);

    S_star = numerator/denominator; 

    //---------------------------------------------------

    //---------------HLLC fluxes -------------------------------------------------//

    //Calculate F_L     
      Euler::U_state F_L;
      F_L = e.flux(U_state_L);

     //Calculate F_L_star
      Euler::U_state F_L_star;
      star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
      
      //---Calculate U_state_L_star
      U_state_L_star.rho = star_coef_left;
      U_state_L_star.moment_u = star_coef_left*S_star;
      U_state_L_star.moment_v = star_coef_left*W_L.v;
      U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
//      /* Consider overloading the + operator to write this in one line 
      F_L_star.rho = F_L.rho + S_L*(U_state_L_star.rho -U_state_L.rho);
      F_L_star.moment_u = F_L.moment_u + S_L*(U_state_L_star.moment_u - U_state_L.moment_u);
      F_L_star.moment_v = F_L.moment_v + S_L*(U_state_L_star.moment_v - U_state_L.moment_v);
      F_L_star.energy = F_L.energy + S_L*(U_state_L_star.energy -U_state_L.energy);
      

      //Calculate F_R
      Euler::U_state F_R;
      F_R = e.flux(U_state_R);

      //Calculate F_R_star
      Euler::U_state F_R_star;
      star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
      
      //Calculate U_state_R_star
      U_state_R_star.rho = star_coef_right;
      U_state_R_star.moment_u = star_coef_right*S_star;
      U_state_R_star.moment_v = star_coef_right*W_R.v;
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));


//      /* Consider overloading the + operator to write this in one line 
      F_R_star.rho = F_R.rho + S_R*(U_state_R_star.rho -U_state_R.rho);
      F_R_star.moment_u = F_R.moment_u + S_R*(U_state_R_star.moment_u -U_state_R.moment_u);
      F_R_star.moment_v = F_R.moment_v + S_R*(U_state_R_star.moment_v -U_state_R.moment_v);
      F_R_star.energy = F_R.energy + S_R*(U_state_R_star.energy -U_state_R.energy);

      //---------------END of HLLC Fluxes, F_L, F_L_star, F_R, F_R_star------------------//
      /*
      if( S_L >= 0.0){
	   flux(f_idx).rho = F_L.rho;
	   flux(f_idx).moment_u = F_L.moment_u;
	   flux(f_idx).moment_v = F_L.moment_v;
	   flux(f_idx).energy = F_L.energy;
	   
      }
      if(S_L <= 0.0 && S_star >= 0.0){
	
	flux(f_idx).rho = F_L_star.rho;
	flux(f_idx).moment_u = F_L_star.moment_u;
	flux(f_idx).moment_v = F_L_star.moment_v;
	   flux(f_idx).energy = F_L_star.energy;
	   
      }
      if(S_star <= 0.0 && S_R >= 0.0){
	flux(f_idx).rho = F_R_star.rho;
	flux(f_idx).moment_u = F_R_star.moment_u;
	flux(f_idx).moment_v = F_R_star.moment_v;
	flux(f_idx).energy = F_R_star.energy;
	
      }
      if( S_R <= 0.0){
	flux(f_idx).rho = F_R.rho;
	flux(f_idx).moment_u = F_R.moment_u;
	flux(f_idx).moment_v = F_R.moment_v;
	flux(f_idx).energy = F_R.energy;
      }

      f_idx++;	

      */


      //Compute the courant number
      c_L = dt*S_L/ds;
      c_star = dt*S_star/ds;
      c_R = dt*S_R/ds;
      c_shear = dt*S_star/ds ;

      //Compute the ratios
      //-----compute dq, dq_up, dq_down for each wave at the interface. Up for 
      
      //dq for the interface we are trying to solve for i+1/2
      dq_l = U_state_L.rho-U_state_L_star.rho;
      dq_star = U_state_L_star.rho - U_state_R_star.rho;
      dq_shear = (e.PfromC(U_state_L)).v-(e.PfromC(U_state_R)).v;
      dq_r = U_state_R_star.rho - U_state_R.rho;

      //std::cout << "About to calculate the HLLC values at the right and left interface to our current interface" << "\n";

      left_interface = HLLC_U_state(input_data(i-1), input_data(i));
      right_interface = HLLC_U_state(input_data(i+1),input_data(i+2));
      
      //std::cout << " calculate left and right interfaces " << "\n";

      dq_l_left_interface = left_interface(0).rho-left_interface(1).rho;
      dq_star_left_interface = left_interface(1).rho - left_interface(2).rho;
      dq_r_left_interface = left_interface(2).rho - left_interface(3).rho;
      dq_shear_left_interface = (e.PfromC(left_interface(0))).v-(e.PfromC(left_interface(3))).v;


      dq_l_right_interface = right_interface(0).rho-right_interface(1).rho;
      dq_star_right_interface = right_interface(1).rho - right_interface(2).rho;
      dq_r_right_interface = right_interface(2).rho - right_interface(3).rho;
      dq_shear_right_interface =  (e.PfromC(right_interface(0))).v-(e.PfromC(right_interface(3))).v;
      
     
      //---------Computing ratios from d_q jumps--------------------//      
      if(dq_l > 1e-10){

	if (c_L > 0){
	  r_L = dq_l_left_interface/dq_l;
	}else if(c_L == 0){
	  r_L = 0;
	}else if(c_L < 0){
	  r_L = dq_l_right_interface/dq_l;
	}else{
	  std::cout << "Issue with c_L : " << c_L << "\n";
	}
      }
      else{
	r_L = 0.0;
      }

      if(dq_star > 1e-10){
	
	if (c_star > 0){
	  r_star = dq_star_left_interface/dq_star;
	}else if(c_star == 0){
	  r_star = 0;
	}else if(c_star < 0){
	  r_star = dq_star_right_interface/dq_star;
	}else{
	  std::cout << "Issue with c_star : " << c_star << "\n";
	}
      }
      else{
	r_star = 0.0;
      }

      if(dq_r > 1e-10){
	
	if (c_R > 0){
	  r_R = dq_r_left_interface/dq_r;
	}else if(c_R==0){
	  r_R = 0.0;
	}else if(c_R < 0){
	  r_R = dq_r_right_interface/dq_r;
	}else{
	  std::cout << "Issue with c_r : " << c_R << "\n";
	}
      }
      else{
	r_R = 0.0;
      } 

//----------Shear jump and ratio ---------------------------
	 
      if(fabs(dq_shear) > 0){
	
	if (c_shear > 0){
	  r_shear = dq_shear_left_interface/dq_shear;
	}else if(c_shear == 0){
	  r_shear = 0.0;
	}else if(c_shear < 0){
	  r_shear = dq_shear_right_interface/dq_shear;
	}else{
	  std::cout << "Issue with c_shear : " << c_shear << "\n";
	}
	
      }
      else{
	r_shear = 0.0;
      } 
      
//----------------------------------------------------------

      //Compute limiter functions
      minmod_l = minmod(r_L,fabs(c_L));
      minmod_star = minmod(r_star, fabs(c_star));
      minmod_r = minmod(r_R,fabs(c_R));
      minmod_shear = minmod(r_shear,fabs(c_shear));

      superbee_l = superbee(r_L, fabs(c_L));
      superbee_star = superbee(r_star, fabs(c_star));
      superbee_r = superbee(r_R, fabs(c_R));
      superbee_shear = superbee(r_shear,fabs(c_shear));

       if(minmod_shear < minmod_star){
	 //std::cout << "Inside WAF_1D. WAF Middle Right state calculation" << "\n";
	 
	 w_middle_state_left.rho = U_state_R_star.rho;
	 w_middle_state_left.u = (e.PfromC(U_state_R_star)).u;
	 w_middle_state_left.v = W_L.v;
	 w_middle_state_left.P = (e.PfromC(U_state_R_star)).P;


	 F_ML = e.flux((e.CfromP(w_middle_state_left)));

	 flux(f_idx).rho =  0.5*(F_L.rho+F_R.rho)-0.5*(sign(c_L)*minmod_l*(F_L_star.rho-F_L.rho) +\
						   sign(c_shear)*minmod_shear*(F_ML.rho-F_L_star.rho) +\
						   sign(c_star)*minmod_star*(F_R_star.rho-F_ML.rho)+\
						   sign(c_R)*minmod_r*(F_R.rho-F_R_star.rho));

	 flux(f_idx).moment_u =  0.5*(F_L.moment_u+F_R.moment_u)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_u-F_L.moment_u) + \
								      sign(c_shear)*minmod_shear*(F_ML.moment_u-F_L_star.moment_u) + \
								      sign(c_star)*minmod_star*(F_R_star.moment_u-F_ML.moment_u)+ \
								      sign(c_R)*minmod_r*(F_R.moment_u-F_R_star.moment_u));

	 flux(f_idx).moment_v =  0.5*(F_L.moment_v+F_R.moment_v)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_v-F_L.moment_v) + \
								      sign(c_shear)*minmod_shear*(F_ML.moment_v-F_L_star.moment_v) + \
								      sign(c_star)*minmod_star*(F_R_star.moment_v-F_ML.moment_v)+ \
								      sign(c_R)*minmod_r*(F_R.moment_v-F_R_star.moment_v));
	 

	 flux(f_idx).energy =  0.5*(F_L.energy+F_R.energy)-0.5*(sign(c_L)*minmod_l*(F_L_star.energy-F_L.energy) + \
								sign(c_shear)*minmod_shear*(F_ML.energy-F_L_star.energy) + \
								sign(c_star)*minmod_star*(F_R_star.energy-F_ML.energy)+	\
								sign(c_R)*minmod_r*(F_R.energy-F_R_star.energy));
       }else if(minmod_shear > minmod_star){
	 // std::cout << "Inside WAF_1D. WAF Middle Left state calculation" << "\n";
	 //std::cout << "AFTER .The shear wave is behind the contact wave." << "\n";
	 //std::cout << "The resulting middle state is rho*l and v_R"<< "\n";

	 w_middle_state_right.rho = U_state_L_star.rho;
	 w_middle_state_right.u = (e.PfromC(U_state_L_star)).u;
	 w_middle_state_right.v = W_R.v;
	 w_middle_state_right.P = (e.PfromC(U_state_L_star)).P;



	 F_MR = e.flux((e.CfromP(w_middle_state_right)));

	 flux(f_idx).rho =  0.5*(F_L.rho+F_R.rho)-0.5*(sign(c_L)*minmod_l*(F_L_star.rho-F_L.rho) + \
						       sign(c_shear)*minmod_shear*(F_MR.rho-F_L_star.rho) + \
						       sign(c_star)*minmod_star*(F_R_star.rho-F_MR.rho)+ \
						       sign(c_R)*minmod_r*(F_R.rho-F_R_star.rho));

	 flux(f_idx).moment_u =  0.5*(F_L.moment_u+F_R.moment_u)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_u-F_L.moment_u) + \
								      sign(c_shear)*minmod_shear*(F_MR.moment_u-F_L_star.moment_u) + \
								      sign(c_star)*minmod_star*(F_R_star.moment_u-F_MR.moment_u)+ \
								      sign(c_R)*minmod_r*(F_R.moment_u-F_R_star.moment_u));

	 flux(f_idx).moment_v =  0.5*(F_L.moment_v+F_R.moment_v)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_v-F_L.moment_v) + \
								  sign(c_shear)*minmod_shear*(F_MR.moment_v-F_L_star.moment_v) + \
								  sign(c_star)*minmod_star*(F_R_star.moment_v-F_MR.moment_v)+ \
								  sign(c_R)*minmod_r*(F_R.moment_v-F_R_star.moment_v));

	 flux(f_idx).energy =  0.5*(F_L.energy+F_R.energy)-0.5*(sign(c_L)*minmod_l*(F_L_star.energy-F_L.energy) + \
								sign(c_shear)*minmod_shear*(F_MR.energy-F_L_star.energy) + \
								sign(c_star)*minmod_star*(F_R_star.energy-F_MR.energy)+	\
								sign(c_R)*minmod_r*(F_R.energy-F_R_star.energy));
	 
	 

       }else if(minmod_shear == minmod_star){
	 
	 flux(f_idx).rho = 0.5*(F_L.rho+F_R.rho)-0.5*(sign(c_L)*minmod_l*(F_L_star.rho-F_L.rho) + \
						      sign(c_star)*minmod_star*(F_R_star.rho-F_L_star.rho) + \
						      sign(c_R)*minmod_r*(F_R.rho-F_R_star.rho));

	 flux(f_idx).moment_u=0.5*(F_L.moment_u+F_R.moment_u)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_u-F_L.moment_u) + \
								   sign(c_star)*minmod_star*(F_R_star.moment_u-F_L_star.moment_u) + \
								   sign(c_R)*minmod_r*(F_R.moment_u-F_R_star.moment_u));

	 flux(f_idx).moment_v = 0.5*(F_L.moment_v+F_R.moment_v)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_v-F_L.moment_v) + \
								     sign(c_star)*minmod_star*(F_R_star.moment_v-F_L_star.moment_v) + \
								     sign(c_R)*minmod_r*(F_R.moment_v-F_R_star.moment_v));
						     

	 flux(f_idx).energy = 0.5*(F_L.energy + F_R.energy) - 0.5* (sign(c_L)*minmod_l*(F_L_star.energy-F_L.energy) + \
								    sign(c_star)*minmod_star*(F_R_star.energy-F_L_star.energy) + \
								    sign(c_R)*minmod_r*(F_R.energy-F_R_star.energy));
	 
	 
       }else{
	 std::cout <<"Something went wrong in the limiter calcualtion" << "\n";       
	 }
	 f_idx++;
  }
  //Logic in case of y-sweep;
  
  if (sweep == std::string("x-sweep")){

  }else if(sweep == std::string("y-sweep")){
    for(int i = 0; i < ncells + 1; i++){
      flux_temp(i).moment_v = flux(i).moment_v;
      flux_temp(i).moment_u = flux(i).moment_u;
      flux(i).moment_v = flux_temp(i).moment_u;
      flux(i).moment_u = flux_temp(i).moment_v;
    }
  }else{
    std::cout <<"Sweep not specified. Fail to compute 1D WAF" << "\n";
  }
  
  return flux;
      
}

int sign(double c){
  
  if(c >= 0){
    return 1;
  }
  else{
    return -1;
  }

}
double minmod(double r, double c){

  
  if (r <= 0){
    return 1.0;
  }
  if( r >= 0 && r <= 1){
    return (1-(1-c)*r);
  }
  if (r > 1 ){
    return c;
  }
  else{
    std::cout  << "Error in minmod limiter calculation " << r << "\n";
    return 0;
    }
}

double superbee(double r, double c){
  
  return 1;

  if(r <= 0){
    return 1;
  }
  else{
    double top = (1-c)*2*r;
    return (1 - top/(1+r));
  }

}

//Modified to account for the exatra speed v in 2D. 
blitz::Array<Euler::U_state,1> HLLC_U_state(Euler::U_state U_state_L, Euler::U_state U_state_R){

  blitz::Array<Euler::U_state,1> U_interface; //vector of size 4 with U_L,U_L_star,U_R_star,U_star;
  
  U_interface.resize(4);
  U_interface(0)=U_state_L;
  U_interface(3)=U_state_R;

  Euler e;
  const double gamma = e.gamma;//This is a massive fudge!!!!! 
  
  //Variables used
  double P_star,P_L,P_R; //Presures
  double rho_star,rho_L,rho_R;//densities
  double u_L,u_R;//particle/gas speed in cell

  double a_L, a_R;// sound speed in cell
  double S_L, S_R,S_star;// wave speed in cell left wave, right wave, contact wave
    
  double P_pvrs;
  double rho_bar, a_bar;//average density and average sound speed of left and right cells

  Euler::W_state w_temp_left;
  Euler::W_state w_temp_right;

  double q_R,q_L;//Part of HLLC Calculation
 
  Euler::W_state W_L;
  Euler::W_state W_R;
  Euler::U_state U_state_L_star;
  Euler::U_state U_state_R_star;

  double star_coef_left; // The coeficient in eq. 10.73 from Toro(ed.2009); (S_k-u_k)/(S_k-u_star_k);
  double star_coef_right; // The coeficient in eq. 10.73 from Toro(ed.2009);
    
  W_L = e.PfromC(U_state_L);
  W_R = e.PfromC(U_state_R);

  //---------Pressure estimate-------------------------

  P_L = W_L.P;
  P_R = W_R.P;
  
  u_L = W_L.u;
  u_R = W_R.u;
  
  rho_L = W_L.rho;
  rho_R = W_R.rho;
  
  a_L = e.a(W_L);
  a_R = e.a(W_R);
      
  rho_bar = 0.5*(rho_L + rho_R);
  a_bar = 0.5*(a_L + a_R);
  
  P_pvrs = 0.5*(P_L + P_R)-0.5*(u_R-u_L)*rho_bar*a_bar;
  
  P_star = std::max(0.0,P_pvrs);
  
  //---------------------------------------------------
  
  //----------Wave speed estimate----------------------
  
  //Calculate q_R
  
  if(P_star <= P_R){
    q_R = 1.0;
  }
  else{
    q_R = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_R)-1));
  }
  
  //Calculate q_L
  if(P_star <= P_L){
    q_L = 1.0;
  }
  else{
    q_L = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_L)-1));
  }
  
  //Calculate S_R and S_L
  S_L = u_L -a_L*q_L;
  S_R = u_R + a_R*q_R;
    
  double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
  double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);
  
  S_star = numerator/denominator; 
  
  //---------------------------------------------------
  
  //---------------HLLC U_state -------------------------------------------------//
  
  //--Calculate U_state_L_star
  star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
  
  //---Calculate U_state_L_star
  U_state_L_star.rho = star_coef_left;
  U_state_L_star.moment_u = star_coef_left*S_star;
  U_state_L_star.moment_v = star_coef_left*W_L.v;
  U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
  U_interface(1)=U_state_L_star;
  
  // Euler::U_state F_R_star;
  star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
  
  //Calculate U_state_R_star
  U_state_R_star.rho = star_coef_right;
  U_state_R_star.moment_u = star_coef_right*S_star;
  U_state_R_star.moment_v = star_coef_right*W_R.v;
  U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));
  U_interface(2)=U_state_R_star;
  
  return U_interface;

}

void flux_and_update(Mesh &m,double dt,std::string sweep_order){

  double dt_dx = dt/m.dx;
  double dt_dy = dt/m.dy;
  std::string x_sweep_output = "x_sweep_output";
  std::string y_sweep_output = "y_sweep_output";
  blitz::Array<Euler::U_state,1> input;
  input.resize(m.ncells+2*m.nGhost);
  blitz::Array<Euler::U_state,1> flux_x(m.ncells+1);
  blitz::Array<Euler::U_state,1> flux_y(m.ncells+1);
  std::string limiter = "minmod";
  std::string sweep;
 
  int f_idx = 0; 

  //X-Sweep,Y-sweep
  if(sweep_order == std::string("XY")){
    
    //X-SWEEP
    sweep = std::string("x-sweep");
  
    for(int i = m.nGhost; i < m.ncells + m.nGhost; i++){
      //Selects all the colums of the i th row.
      input = m.Bdata(blitz::Range::all(),i);
      //Calculate x-fluxes
      flux_x = WAF_1D(input,dt,m.dx,m.ncells,m.nGhost,limiter,sweep); 
      f_idx = 0;
    //Update solution U_dt+1/2
      for(int col= m.nGhost; col < m.ncells + m.nGhost; col++){
	m.Bdata(col,i).rho = m.Bdata(col,i).rho - dt_dx*(flux_x(f_idx+1).rho-flux_x(f_idx).rho);
	m.Bdata(col,i).moment_u = m.Bdata(col,i).moment_u - dt_dx*(flux_x(f_idx+1).moment_u - flux_x(f_idx).moment_u);
	m.Bdata(col,i).moment_v = m.Bdata(col,i).moment_v - dt_dx*(flux_x(f_idx+1).moment_v - flux_x(f_idx).moment_v);
	m.Bdata(col,i).energy = m.Bdata(col,i).energy - dt_dx*(flux_x(f_idx+1).energy-flux_x(f_idx).energy);
	f_idx++;
      }

    }  
    
    m.applyBC();
   
    //Y_SWEEP
    sweep = std::string("y-sweep");
    for(int i = m.nGhost; i < m.ncells+m.nGhost; i++){
      input = m.Bdata(i, blitz::Range::all());
      flux_y = WAF_1D(input,dt,m.dy,m.ncells,m.nGhost,limiter,sweep); 
      f_idx = 0;
      for(int row = m.nGhost; row < m.ncells + m.nGhost; row++){
	m.Bdata(i,row).rho = m.Bdata(i,row).rho - dt_dy*(flux_y(f_idx+1).rho-flux_y(f_idx).rho);
	m.Bdata(i,row).moment_u = m.Bdata(i,row).moment_u - dt_dy*(flux_y(f_idx+1).moment_u - flux_y(f_idx).moment_u);
	m.Bdata(i,row).moment_v = m.Bdata(i,row).moment_v - dt_dy*(flux_y(f_idx+1).moment_v - flux_y(f_idx).moment_v);
	m.Bdata(i,row).energy = m.Bdata(i,row).energy - dt_dy*(flux_y(f_idx+1).energy-flux_y (f_idx).energy);
	f_idx++;
      }
    }  

  }else if(sweep_order == std::string("YX")){

       //Y_SWEEP
       sweep = std::string("y-sweep");
       for(int i = m.nGhost; i < m.ncells+m.nGhost; i++){
	 input = m.Bdata(i, blitz::Range::all());
	 flux_y = WAF_1D(input,dt,m.dy,m.ncells,m.nGhost,limiter,sweep); 
	 f_idx=0;
	 for(int row = m.nGhost; row < m.ncells + m.nGhost; row++){

	  m.Bdata(i,row).rho = m.Bdata(i,row).rho - dt_dy*(flux_y(f_idx+1).rho-flux_y(f_idx).rho);
	  m.Bdata(i,row).moment_u = m.Bdata(i,row).moment_u - dt_dy*(flux_y(f_idx+1).moment_u - flux_y(f_idx).moment_u);
	  m.Bdata(i,row).moment_v = m.Bdata(i,row).moment_v - dt_dy*(flux_y(f_idx+1).moment_v - flux_y(f_idx).moment_v);
	  m.Bdata(i,row).energy = m.Bdata(i,row).energy - dt_dy*(flux_y(f_idx+1).energy-flux_y (f_idx).energy);
	   f_idx++;
	 }
    
       }  

       m.applyBC();

       //X_SWEEP
       sweep = std::string("x-sweep");
       for(int i = m.nGhost; i < m.ncells+m.nGhost; i++){
	 input = m.Bdata(blitz::Range::all(),i);
	 flux_x = WAF_1D(input,dt,m.dx,m.ncells,m.nGhost,limiter,sweep); 
	 f_idx=0;
	 //Update solution U_dt+1/2
	 for(int col= m.nGhost; col < m.ncells + m.nGhost; col++){
	  m.Bdata(col,i).rho = m.Bdata(col,i).rho - dt_dx*(flux_x(f_idx+1).rho-flux_x(f_idx).rho);
	  m.Bdata(col,i).moment_u = m.Bdata(col,i).moment_u - dt_dx*(flux_x(f_idx+1).moment_u - flux_x(f_idx).moment_u);
	  m.Bdata(col,i).moment_v = m.Bdata(col,i).moment_v - dt_dx*(flux_x(f_idx+1).moment_v - flux_x(f_idx).moment_v);
	  m.Bdata(col,i).energy = m.Bdata(col,i).energy - dt_dx*(flux_x(f_idx+1).energy-flux_x(f_idx).energy);
	   f_idx++;
	 }
       }
       
      }else{
       std::cout << "Inside flux_and_update function. Wrong sweep input" << "\n";
       
     }
}
//Solid phi > 0, Fluid phi < 0;
void Mesh::applyGhost(double obj_vel_x, double obj_vel_y){
  std::cout << "Inside applyGhost() function" << "\n";
  double d_tau= dx/2;
  double n_x,n_y;//x and y components of normal to level_set
  Euler::W_state w,w_xdir_plus,w_xdir_minus,w_ydir_plus,w_ydir_minus;
  double rho_xdir,rho_ydir,u_xdir,u_ydir,v_xdir,v_ydir,P_xdir,P_ydir;
  double u_t,u_n,v_t,v_n;
  double vel_normal;
  double obj_vel_x_n, obj_vel_y_n, obj_vel_x_t, obj_vel_y_t, vel_normal_obj;
  //Loop over tau - iteration to fill in solid
  for (int tau = 0; tau < 20; tau++){       
    //Loop over rows
    for(int col = nGhost ; col < nGhost+ncells ; col++ ){
      //Loop over columns
      for(int row = nGhost; row < nGhost+ncells; row++){
	//Chek level set
	//std::cout << "The phi value at the beginging is : "<< Bdata(col,row).phi << "\n";
	if (MLS_data(col,row).phi >= 0){
	  //calculate norms, n_x,n_y
	  n_x = ((MLS_data(col+1,row).phi-MLS_data(col-1,row).phi)/(2*dx));
	  n_y = (MLS_data(col,row+1).phi-MLS_data(col,row-1).phi)/(2*dy);
	  w = ptr_euler->PfromC(Bdata(col,row));
	  w_xdir_plus=ptr_euler->PfromC(Bdata(col+1,row));
	  w_xdir_minus=ptr_euler->PfromC(Bdata(col-1,row));
	  w_ydir_plus=ptr_euler->PfromC(Bdata(col,row+1));
	  w_ydir_minus=ptr_euler->PfromC(Bdata(col,row-1));

	  //	  std::cout << "The value of the normal n_x and n_y is : " << n_x << "\t" << n_y << "\n";
	  //Extrapolate : Calculate I_x, I_y, where I= rho, u,v,P
	  
	  //Extrapolate density
	  if(n_x <=0){

	    rho_xdir = (w_xdir_plus.rho-w.rho)/dx; 

	  }else if(n_x >0){
	    rho_xdir = (w.rho-w_xdir_minus.rho)/dx; 
	 
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_X" <<"\n";
	  }
	  
	  if(n_y <=0){
	    rho_ydir = (w_ydir_plus.rho-w.rho)/dy;
	  }else if(n_y >0){
	    rho_ydir = (w.rho-w_ydir_minus.rho)/dy;
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_Y" <<"\n";
	  }
	  
	  w.rho = w.rho - d_tau*(n_x*rho_xdir+n_y*rho_ydir);
	  
	  //Extrapolate u-Bdata is in U_state format! 
	  //Need to find normal and tangential speed and to invert the normal componenet
	  vel_normal_obj= n_x*obj_vel_x + n_y*obj_vel_y;

	  obj_vel_x_n = vel_normal_obj*n_x;
	  obj_vel_x_t = obj_vel_x - obj_vel_x_n;
	  obj_vel_x   = obj_vel_x_n + obj_vel_x_t;
	  
	  obj_vel_y_n = vel_normal_obj*n_y;
	  obj_vel_y_t = obj_vel_y - obj_vel_y_n;
	  obj_vel_y   = obj_vel_y_n + obj_vel_y_t;
	  
	  vel_normal = n_x*w.u + n_y*w.v;

	  u_n = vel_normal*n_x;// - 2*obj_vel_x_n;//compute normal component
	  u_t = w.u-u_n;//compute tangential component
	  w.u = u_n + u_t;// + obj_vel_x;

	  v_n = vel_normal*n_y;// - 2*obj_vel_y_n;
	  v_t = w.v-v_n;
	  w.v = v_n + v_t;// + obj_vel_y;
	  
	  if(MLS_data(col+1,row).phi < 0){
	    vel_normal=n_x*w_xdir_plus.u-n_y*w_ydir_plus.v;
	    
	    u_n = vel_normal*n_x;//compute normal component
	    u_t = w_xdir_plus.u-u_n;//compute tangential component
	    u_n = -1*u_n + 2*obj_vel_x_n;//Reflect normal component
	    w_xdir_plus.u = u_n+u_t;

	    v_n = vel_normal*n_y;//compute normal component
	    v_t = w_xdir_plus.v-v_n;//compute tangential component
	    v_n = -1*v_n + 2*obj_vel_y_n;//Reflect normal component
	    w_xdir_plus.v = v_n+v_t;

	  }

	  if(MLS_data(col-1,row).phi < 0){
	    vel_normal = w_xdir_minus.u*n_x + w_xdir_minus.v*n_y;
	    u_n = vel_normal*n_x;//compute normal component
	    u_t = w_xdir_minus.u-u_n;//compute tangential component
	    u_n = -1*u_n + 2*obj_vel_x_n;//Reflect normal component
	    w_xdir_minus.u = u_n+u_t;

	    v_n = vel_normal*n_y;
	    v_t = w_xdir_minus.v - v_n;
	    v_n = -1*v_n + 2*obj_vel_y_n ;
	    w_xdir_minus.v = v_n+v_t;
	  }

	  if(MLS_data(col,row+1).phi < 0 ){
	    vel_normal = w_ydir_plus.u*n_x + w_ydir_plus.v*n_y;
	    
	    u_n = vel_normal*n_x;//compute normal component
	    u_t = w_ydir_plus.u-u_n;//compute tangential component
	    u_n = -1*u_n + 2*obj_vel_x_n;//Reflect normal component
	    w_ydir_plus.u = u_n+u_t;
	 
	    v_n = vel_normal*n_y;
	    v_t = w_ydir_plus.v-v_n;
	    v_n = -1*v_n + 2*obj_vel_y_n;
	    w_ydir_plus.v = v_n+v_t;
	 
	  }
	  if(MLS_data(col,row-1).phi < 0){
	    vel_normal = w_ydir_minus.u*n_x + w_ydir_minus.v*n_y;
	   
	    u_n = vel_normal*n_x;//compute normal component
	    u_t = w_ydir_minus.u-u_n;//compute tangential component
	    u_n = -1*u_n+ 2*obj_vel_x_n;//Reflect normal component
	    w_ydir_minus.u = u_n+u_t;
	    
	    v_n = vel_normal*n_y;
	    v_t = w_ydir_minus.v-v_n;
	    v_n = -1*v_n+ 2*obj_vel_y_n;
	    w_ydir_minus.v = v_n+v_t;
	  }
	  
	  if(n_x <=0){

	    u_xdir = (w_xdir_plus.u-w.u)/dx; 

	  }else if(n_x >0){
	    u_xdir = (w.u-w_xdir_minus.u)/dx; 
	 
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_X" <<"\n";
	  }
	  
	  if(n_y <=0){
	    u_ydir = (w_ydir_plus.u-w.u)/dy;
	  }else if(n_y >0){
	    u_ydir = (w.u-w_ydir_minus.u)/dy;
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_Y" <<"\n";
	  }
	  
	  w.u = w.u - d_tau*(n_x*u_xdir+n_y*u_ydir);
	 
	  
	  //Extrapolate v-Bdata is in U_state format!
	  if(n_x <=0){

	    v_xdir = (w_xdir_plus.v-w.v)/dx; 

	  }else if(n_x >0){
	    v_xdir = (w.v-w_xdir_minus.v)/dx; 
	 
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_X" <<"\n";
	  }
	  
	  if(n_y <=0){
	    v_ydir = (w_ydir_plus.v-w.v)/dy;
	  }else if(n_y >0){
	    v_ydir = (w.v-w_ydir_minus.v)/dy;
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_Y" <<"\n";
	  }
	  
	  w.v = w.v - d_tau*(n_x*v_xdir+n_y*v_ydir);
	 
	  
	  //Extrapolate P
	  if(n_x <=0){

	    P_xdir = (w_xdir_plus.P-w.P)/dx; 

	  }else if(n_x >0){
	    P_xdir = (w.P-w_xdir_minus.P)/dx; 
	 
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_X" <<"\n";
	  }
	  
	  if(n_y <=0){
	    P_ydir = (w_ydir_plus.P-w.P)/dy;
	  }else if(n_y >0){
	    P_ydir = (w.P-w_ydir_minus.P)/dy;
	  }else{
	    std::cout << "Extrapolation of rho. Something went wrong in calculating N_Y" <<"\n";
	  }
	  
	  w.P = w.P - d_tau*(n_x*P_xdir+n_y*P_ydir);
	  

	  Bdata(col,row) = ptr_euler->CfromP(w);
        
	}else{
	  // std::cout << "Inside applyGhost. This is not a ghost cell " << "\n";
	}
	//std::cout << "The phi value at the end is : "<< Bdata(col,row).phi << "\n";
	
      }
    }
  }
}
/*
void Mesh::update_level_set(double dt, double time, double obj_speed_x, double obj_speed_y,double obj_x_c, double obj_y_c, double (*level_set)(double x, double y, double x_c, double y_c)){
  
  static double x_new,y_new;
  
  if (time == 0){
    y_new = obj_y_c + dt*obj_speed_y;
    x_new = obj_x_c + dt*obj_speed_x;
  }else{
    y_new = y_new +dt*obj_speed_y;
    x_new = x_new +dt*obj_speed_x;
  }
  
  for(int row = nGhost; row < nGhost + ncells; row++ ){
    for(int col= nGhost; col < nGhost + ncells; col++){
      Bdata(row,col).phi = level_set(xaxis(row),yaxis(col), x_new,y_new);
    }
  }
}*/


void Mesh::advect_level_set(double dt, double time, double obj_speed_x,double obj_speed_y){

  std::cout << "ADVECTING" << "\n";

  double phi;
  double phi_xdir,phi_ydir;
  double phi_xdir_plus,phi_xdir_minus;
  double phi_ydir_plus,phi_ydir_minus;
  blitz::Array<double,2>phi_new;
  phi_new.resize(ncells+2*nGhost,ncells+2*nGhost);

  for(int row = nGhost; row < nGhost+ncells; row++){
    for(int col = nGhost; col < nGhost+ncells; col++){

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
      
      if(obj_speed_y <0){
	phi_ydir = (phi_ydir_plus-phi)/dy;
      }else if(obj_speed_y == 0){
	phi_ydir = 0;
      }else if(obj_speed_y >0){
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

}
