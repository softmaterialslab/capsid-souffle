//
// Created by lauren on 1/25/18.
//

#ifndef SOUFFLE_BEAD_H
#define SOUFFLE_BEAD_H

#include <vector>
#include <cmath>
#include "vector3d.h"
#include "thermostat.h"
#include "rand_gen.h"
#include "oligomer.h"

class EDGE;                               //forward declarations
class SUBUNIT;
class FACE;
class BOX;
//class OLIGOMER;




class BEAD {
public:
//member variables
   unsigned int id ;                     //particle ID
   double m;                             //mass of the particle (amu)
   double sigma;                         //diameter of particle (lj reduced unit)
   VECTOR3D bx;                          //vector of box lengths (lj reduced unit)
   VECTOR3D hbx;                         //vector of half box lengths
   int unit ;                            //id of the subunit it belongs to
   int type ;                            //id of type of bead
   VECTOR3D pos;                         //position of the particle (xyz) (lj reduced unit)
   VECTOR3D oldpos;
   VECTOR3D vel;                         //velocity of the particle (xyz) (delt^-1)
   VECTOR3D sforce;                      //forces on the particle (xyz) (amu/delt^2)
   VECTOR3D oldtforce;
   VECTOR3D tforce;
   VECTOR3D ljforce;
   VECTOR3D bforce;
   VECTOR3D eforce;
   VECTOR3D fdrag;
   VECTOR3D fran;
   double se;                            //stretching energy (KbT)
   double ke;                            //kinetic energy (KbT)
   double be;                            //bending energy (KbT)
   double ne;                            //lennard jones energy (KbT)
   double ce;                            //electrostatic energy (KbT)
   double q;                             //charge
   std::vector<EDGE*> itsE;              //its Edges
   std::vector<SUBUNIT*> itsS;           //its Subunit(s)
   std::vector<FACE*> itsF;              //its Faces
   std::vector<int> itsN;                //Its neighbors (pairlist)
   VECTOR3D noise;                       //noise term -- gaussian random number
   std::vector<OLIGOMER> itsSO;         //its sub_oligomer
   OLIGOMER* itsO;
 
//member functions
//constructor fxn
   BEAD(VECTOR3D position_i=VECTOR3D(0,0,0),double mi=0, int id_i=0, int unit_i=0, int type_i=0)
   {                                                       //constructor fxn
      m = mi;
      pos = position_i;
      id = id_i;
      unit = unit_i;
      type = type_i;
   }
   
//position updated full step
   void update_position(double delt){ 
      pos = ( pos + (vel ^ delt) );
      // periodic boundary is accounted for
      if (pos.x > hbx.x) pos.x = pos.x - bx.x;
      else if (pos.x < -hbx.x) pos.x = pos.x + bx.x;
      if (pos.y > hbx.y) pos.y = pos.y - bx.y;
      else if (pos.y < -hbx.y) pos.y = pos.y + bx.y;
      if (pos.z > hbx.z) pos.z = pos.z - bx.z;
      else if (pos.z < -hbx.z) pos.z = pos.z + bx.z;
    }
    
    //position updated full step
    void update_position_brownian(double delta_t, double c2, double fric_zeta){ 
      // pos =  pos + (vel + (tforce ^ (delta_t * (0.5 / m))) + (noise / m)) ^ (delta_t / c2);
      // pos =  pos + ((vel + noise) ^ (delta_t / c2));
       oldpos = pos;
       pos =  pos + (vel ^ (delta_t / c2)) + (noise ^ (0.5 * delta_t / (m * c2))) + (tforce ^ (delta_t * delta_t * 0.5 / (c2 * m)));
       pos = (pos ^ (1 - ((fric_zeta * delta_t) / (c2 * m)) )) + (vel ^ (delta_t / c2));
       // periodic boundary is accounted for
       if (pos.x > hbx.x) pos.x = pos.x - bx.x;
       else if (pos.x < -hbx.x) pos.x = pos.x + bx.x;
       if (pos.y > hbx.y) pos.y = pos.y - bx.y;
       else if (pos.y < -hbx.y) pos.y = pos.y + bx.y;
       if (pos.z > hbx.z) pos.z = pos.z - bx.z;
       else if (pos.z < -hbx.z) pos.z = pos.z + bx.z;
    }

//velocity updated half timestep (No thermostat)
   void update_velocity(double delt){   
      vel = ( vel + ( (tforce) ^ (0.5*delt/m) ) );
   }

// update velocity with integrator that unifies velocity verlet and Nose-Hoover
   void therm_update_velocity(double dt, THERMOSTAT main_bath, double expfac){
      vel = ( ( vel ^ (expfac)  ) + ( tforce ^ (0.5 * dt * (std::sqrt(expfac)) / m) ) );
      //vel = ( ( vel ^ (expfac * expfac)  ) + ( tforce ^ (0.5 * dt * (expfac) / m) ) );
      return;
   }

   void brownian_update_velocity(double delt, double fric_zeta){
      RAND_GEN ugsl;
      vel.x += (vel.x*(-0.5*fric_zeta*delt)) + (tforce.x*(0.5*delt/m)) +
                        sqrt(2*6*delt*fric_zeta/m)*(gsl_rng_uniform(ugsl.r)-0.5);
      vel.y += (vel.y*(-0.5*fric_zeta*delt)) + (tforce.y*(0.5*delt/m)) +
                        sqrt(2*6*delt*fric_zeta/m)*(gsl_rng_uniform(ugsl.r)-0.5);
      vel.z += (vel.z*(-0.5*fric_zeta*delt)) + (tforce.z*(0.5*delt/m)) +
                        sqrt(2*6*delt*fric_zeta/m)*(gsl_rng_uniform(ugsl.r)-0.5);
   }

   void update_kinetic_energy(){
      ke = 0.5 * m * vel.GetMagnitude() * vel.GetMagnitude();
   }

   void update_stretching_energy(double ks);

   void update_stretching_force(double ks);
   
   void compute_fdrag(double damp){
       fdrag = vel*(-1.0 * m / damp);
   }

   void compute_fran(double delta_t, double damp, double rand_num){
       fran.x = sqrt((m*24.0)/(damp * delta_t)) * rand_num;
       fran.y = sqrt((m*24.0)/(damp * delta_t)) * rand_num;
       fran.z = sqrt((m*24.0)/(damp * delta_t)) * rand_num;
   }

   void update_tforce(){
       tforce = tforce + fran + fdrag;
   }


};


#endif //SOUFFLE_BEAD_H
