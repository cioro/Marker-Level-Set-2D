#include"Euler2D.hpp"
#include<vector>
#include<cmath>
#include<iostream>

Euler::W_state::W_state() : rho(0), u(0), v(0), P(0) {}

Euler::W_state::W_state(double rho_arg, double u_arg, double v_arg, double P_arg) : rho(rho_arg), u(u_arg),v(v_arg), P(P_arg) {}

void Euler::W_state::print()const {
  std::cout << " rho\t " << rho << "\t u \t " << u << "\t v \t" << v << " \t P \t" << P << "\n";

}
double Euler::W_state::get_var(VARIABLE v){
  switch(v){
  case RHO:
    return rho;
    break;
  case U:
    return u;
    break;
  case V:
    return v;
    break; 
  case PRESSURE:
    return P;
    break;
  default:
    return -1;
  }

}

Euler::U_state::U_state() : rho(0), moment_u(0), moment_v(0), energy(0) {}

Euler::U_state::U_state(double rho_arg, double moment_u_arg,double moment_v_arg, double energy_arg) : rho(rho_arg), moment_u(moment_u_arg),moment_v(moment_v_arg), energy(energy_arg) {}

void Euler::U_state::print()const {
  std::cout << " rho \t " << rho << "\t moment_u \t " << moment_u \
	    <<"\t moment_v \t "<< moment_v << " \t energy \t" << energy << "\n";

}

double Euler::U_state::get_var(VARIABLE v){
  switch(v){
  case RHO:
    return rho;
    break;
  case MOMENT_U:
    return moment_u;
    break;
  case MOMENT_V:
    return moment_v;
    break; 
  case ENERGY:
    return energy;
    break;
  default:
    return -1;
  }
}
double Euler::a(){
return 0;
}

double Euler::a(const Euler::W_state& w){
  if(w.P < 0.0 || w.rho < 0.0){
    std::cout << "Error speed  of sound can't be computed. Negative pressure or density" << "\n";
    w.print();
    exit(1);
      }
  double a_result = sqrt(gamma*(w.P)/(w.rho));
  return a_result;
}
  
double Euler::int_energy(const Euler::W_state& w){

  double e = w.P/((gamma-1)*w.rho);
  return e;
  
}

//The Formula come from Toro(ed.2009) p.89
Euler::U_state Euler::flux(const Euler::U_state& U){
  double Pressure = ((gamma-1)*(U.energy-0.5*(U.moment_u*U.moment_u + U.moment_v*U.moment_v)/U.rho));
  double f1 = U.moment_u;
  double f2 = U.moment_u*U.moment_u/U.rho + Pressure;
  double f3 = U.moment_u*U.moment_v/U.rho;
  double f4 = (U.moment_u/U.rho)*(U.energy + Pressure);

  return U_state(f1,f2,f3,f4);
  
}


Euler::U_state Euler::flux_y(const Euler::U_state& U){
 
  double Pressure = ((gamma-1)*(U.energy-0.5*(U.moment_u*U.moment_u + U.moment_v*U.moment_v)/U.rho));
  double f1 = U.moment_v;
  double f2 = U.moment_v*U.moment_u/U.rho;
  double f3 = U.moment_v*U.moment_v/U.rho + Pressure;
  double f4 = (U.moment_v/U.rho)*(U.energy + Pressure);

  return U_state(f1,f2,f3,f4);
  
}


Euler::U_state Euler::CfromP(const Euler::W_state& w){
  double u1 = w.rho;
  double u2 = w.rho*w.u;
  double u3 = w.rho*w.v;
  double u4 = w.rho*(0.5*(w.u*w.u + w.v*w.v) +int_energy(w));

  return U_state(u1,u2,u3,u4);

}


Euler::W_state Euler::PfromC(const Euler::U_state& u){
  double w1 = u.rho;
  double w2 = u.moment_u/u.rho;
  double w3 = u.moment_v/u.rho;
  double w4 = (gamma - 1)*(u.energy - 0.5*((u.moment_u*u.moment_u + u.moment_v*u.moment_v)/u.rho));
 
  return W_state(w1,w2,w3,w4);

}

Euler::Euler(){
  U_state();
  W_state();
}

