//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_BEAD_H
#define LEMONSOUFFLE_BEAD_H

#include <vector>
#include <cmath>
#include "vector3d.h"
#include "thermostat.h"
#include "rand_gen.h"


class EDGE;
class SUBUNIT;
class FACE;
class BOX;
class OLIGOMER;
class PAIR;


class BEAD
{
public:

//member variables
    int id ;                         //particle ID
    double m;                           //mass of the particle (amu)
    double sigma;                       //diameter of particle (unitless)
    VECTOR3D bx;
    int unit ;
    int type ;
    VECTOR3D cell_id;                          //cell array index
    std::vector<int> itsN;                    //Its neighbors
    VECTOR3D pos;                       //position of the particle (xyz) (unitless)
    VECTOR3D vel;                       //velocity of the particle (xyz) (delt^-1)
    VECTOR3D sforce;                    //force on the particle (xyz) (amu/delt^2)
    VECTOR3D tforce;
    VECTOR3D ljforce;
    VECTOR3D bforce;
    VECTOR3D eforce;
    double se;                        //potential energy (KbT)
    double ke;                        //kinetic energy (KbT)
    double be;                        //bending energy (KbT)
    double ne;                        //lennard jones energy (KbT)
    double ce;
    double q;                           //charge
    std::vector<EDGE*> itsE;            
    std::vector<SUBUNIT*> itsS;           
    std::vector<FACE*> itsF;            
    OLIGOMER* itsO;
	std::vector<PAIR*> itsP;			// LJ pairs
 


//member functions

    BEAD(VECTOR3D position_i=VECTOR3D(0,0,0),double mi=0, int id_i=0, int unit_i=0, int type_i=0, OLIGOMER* itsO_i=NULL)
    {                                                       //constructor fxn
        m = mi;
        pos = position_i;
		id = id_i;
		unit = unit_i;
		type = type_i;
		itsO = itsO_i;
    }

    void update_position(double delt)                       //position updated full timestep
    {
        pos = ( pos + (vel ^ delt) );
        // periodic boundary
        if (pos.x > bx.x/2.0)
            pos.x = pos.x - bx.x;
        if (pos.x < -bx.x/2.0)
            pos.x = pos.x + bx.x;
        if (pos.y > bx.y/2.0)
            pos.y = pos.y - bx.y;
        if (pos.y < -bx.y/2.0)
            pos.y = pos.y + bx.y;
        if (pos.z > bx.z/2.0)
            pos.z = pos.z - bx.z;
        if (pos.z < -bx.z/2.0)
            pos.z = pos.z + bx.z;
    }

    void update_velocity(double delt)                       //velocity updated half timestep (No thermostat)
    {
        //tforce = sforce + bforce + ljforce + eforce;
        vel = ( vel + ( (tforce) ^ (0.5*delt/m) ) );
    }

    // update velocity with integrator that unifies velocity verlet and Nose-Hoover
    void therm_update_velocity(double dt, THERMOSTAT main_bath, double expfac)
    {
        //tforce = sforce + bforce + ljforce + eforce;
        vel = ( ( vel ^ (expfac)  ) + ( tforce ^ (0.5 * dt * (std::sqrt(expfac)) / m) ) );
        return;
    }

    void brownian_update_velocity(double delt, double fric_zeta){
        RAND_GEN ugsl;
        //tforce = sforce + bforce + ljforce + eforce;
        vel.x += (vel.x*(-0.5*fric_zeta*delt)) + (tforce.x*(0.5*delt/m)) +
                         sqrt(2*6*delt*fric_zeta/m)*(gsl_rng_uniform(ugsl.r)-0.5);
        vel.y += (vel.y*(-0.5*fric_zeta*delt)) + (tforce.y*(0.5*delt/m)) +
                         sqrt(2*6*delt*fric_zeta/m)*(gsl_rng_uniform(ugsl.r)-0.5);
        vel.z += (vel.z*(-0.5*fric_zeta*delt)) + (tforce.z*(0.5*delt/m)) +
                         sqrt(2*6*delt*fric_zeta/m)*(gsl_rng_uniform(ugsl.r)-0.5);
    }

    void update_kinetic_energy()                             //kinetic energy of the particle
    {
        ke = 0.5 * m * vel.GetMagnitude() * vel.GetMagnitude();
    }

    void update_stretching_energy(double ks, double bondlength);

    void update_stretching_force(double ks, double bondlength);


};


#endif //LEMONSOUFFLE_BEAD_H
