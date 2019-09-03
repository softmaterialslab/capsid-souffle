#include<iostream>
#include <fstream>
#include <iomanip>
#include"energies.h"
#include "bead.h"
#include "LJpair.h"
#include "subunit.h"
#include "edge.h"
#include "face.h"
#include "functions.h"

using namespace std;

long double particle_kinetic_energy(vector <BEAD> &subunit_bead) {          //part of thermostat
	for (unsigned int i = 0; i < subunit_bead.size(); i++)
		subunit_bead[i].update_kinetic_energy();
	long double kinetic_energy = 0.0;
	for (unsigned int i = 0; i < subunit_bead.size(); i++)
		kinetic_energy += subunit_bead[i].ke;
	return kinetic_energy;
}

void update_LJ_energies_simplified(vector<BEAD>& subunit_bead, double ecut, vector<vector<int> > lj_a, double elj_att){
	 
	 for (unsigned int i = 0; i < subunit_bead.size(); i++){
		 for (unsigned int j = 0; j < subunit_bead.size(); j++){
			 if (subunit_bead[i].itsS[0]->id != subunit_bead[j].itsS[0]->id){
                      
				 VECTOR3D r_vec = dist( &subunit_bead[i] , &subunit_bead[j] );
				 long double r = r_vec.GetMagnitude();
				 double r6 ;
				 double sigma6;
                                 double sigma12;
				 
				 double sig1 = subunit_bead[i].sigma;
				 double sig2 = subunit_bead[j].sigma;
				 double shc = (sig1 + sig2)/2;
				 double del = (sig1+sig2)/2 - shc;
				 bool lj_attractive = false;
				 for (unsigned int k = 0; k < lj_a[0].size(); k++){
					 if ( subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]  ){
						 lj_attractive = true;
					 }
				 }
				 if ( r < (del+1.12246205*shc) && lj_attractive == false){							//Repulsive
					 sigma6 = shc * shc * shc * shc * shc * shc;
					 double elj = 1;//subunit_bead[j].epsilon;
					 r6 = (r-del) * (r-del) * (r-del) * (r-del) * (r-del) * (r-del);
					 subunit_bead[i].ne += 0.5*((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
				 } else if ( r < (ecut) && lj_attractive == true ){			//Attractive
					 double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
					 double ecut12 = ecut6 * ecut6;
					 sigma6 = sig1 * sig1 * sig1 * sig1 * sig1 * sig1;
                                         sigma12 = sigma6 * sigma6;
					 r6 = (r-del) * (r-del) * (r-del) * (r-del) * (r-del) * (r-del);
					// double elj = 1.8;//subunit_bead[j].epsilon;
					 subunit_bead[i].ne += 0.5*(((4 * elj_att * (sigma6 / r6) * ((sigma6 / r6) - 1)) - (4 * elj_att * ((sigma12 / ecut12) - (sigma6 / ecut6)))));
				 } else {
					 subunit_bead[i].ne += 0;
				 }
			 }
		 } // for j
	 } // for i
 } // energy fxn
 
 void update_ES_energies_simplified(vector<BEAD>& subunit_bead, double lb, double ni, double qs, double ecut_el, double kappa)
 {
	 
	 for (unsigned int i = 0; i < subunit_bead.size(); i++)
	 {
		 for (unsigned int j = 0; j < subunit_bead.size(); j++)
		 {
			 if (subunit_bead[i].id != subunit_bead[j].id && subunit_bead[i].q != 0 && subunit_bead[j].q != 0)
			 {
			    VECTOR3D r_vec = dist( & subunit_bead[i] , & subunit_bead[j] );
			    long double r = r_vec.GetMagnitude();
			    if ( r < ecut_el && ecut_el != 0){
				 subunit_bead[i].ce += 0.5*( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) ) / (r) -
				                              ((subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*ecut_el) ) / (ecut_el)) );
			    } else if (ecut_el == 0) {
			      subunit_bead[i].ce += 0.5*( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) ) / (r) );
			    }
				 
			 }
		 }
	 }
	 
 }
 
 
 
 void update_LJ_ES_energies_simplified(vector<BEAD>& subunit_bead, double ecut, vector<vector<int> > lj_a, double elj_att, double lb, double ni, double qs, double ecut_el, double kappa){
    
   double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
   double ecut12 = ecut6 * ecut6;
   VECTOR3D box = subunit_bead[0].bx;
   VECTOR3D hbox = subunit_bead[0].hbx;
    
   for (unsigned int i = 0; i < subunit_bead.size(); i++){
      for (unsigned int j = i+1; j < subunit_bead.size(); j++){
          
         bool electrostatic = (i != j && subunit_bead[i].q != 0 && subunit_bead[j].q != 0);
         bool lj =(subunit_bead[i].itsS[0]->id != subunit_bead[j].itsS[0]->id);
          
         long double r = 0.0;
         VECTOR3D r_vec;
          
         if (electrostatic || lj){
            r_vec = subunit_bead[i].pos - subunit_bead[j].pos;
            if (r_vec.x > hbox.x) r_vec.x -= box.x;
            else if (r_vec.x < -hbox.x) r_vec.x += box.x;
            if (r_vec.y > hbox.y) r_vec.y -= box.y;
            else if (r_vec.y < -hbox.y) r_vec.y += box.y;
            if (r_vec.z > hbox.z) r_vec.z -= box.z;
            else if (r_vec.z < -hbox.z) r_vec.z += box.z;
            r = r_vec.GetMagnitude();
         }
          
         if (electrostatic && r < ecut_el && ecut_el != 0) {                 //electrostatic energy
            subunit_bead[i].ce += ( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) ) / (r) -
                                   ( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*ecut_el) ) / (ecut_el)) );
         } else if (electrostatic && ecut_el == 0) {
            subunit_bead[i].ce += ( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) ) / (r) );
         }
         
         double sig1 = subunit_bead[i].sigma;
         double sig2 = subunit_bead[j].sigma;
         double shc = (sig1 + sig2)/2;
         
         if (r >= ecut && r >= (1.12246205 * shc)) continue;  
         
         if (lj) {
            bool lj_attractive = false;
            for (unsigned int k = 0; k < lj_a[0].size(); k++){
               if ( subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]  ){
                  lj_attractive = true;
               }
            }
            if (r < (1.12246205*shc) && lj_attractive == false){                                      //Repulsive
               double sigma6 = shc * shc * shc * shc * shc * shc;
               double elj = 1;//subunit_bead[j].epsilon;
               double r6 = r * r * r * r * r * r;
               subunit_bead[i].ne += ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
            } else if ( r < ecut && lj_attractive == true ){               //Attractive
               double sigma6 = sig1 * sig1 * sig1 * sig1 * sig1 * sig1;
               double sigma12 = sigma6 * sigma6;
               double r6 = r * r * r * r * r * r;
               subunit_bead[i].ne += (((4 * elj_att * (sigma6 / r6) * ((sigma6 / r6) - 1)) - (4 * elj_att * ((sigma12 / ecut12) - (sigma6 / ecut6)))));
            } else {
               subunit_bead[i].ne += 0;
            }
         } //if (lj)
         
      } // for j
   } // for i
   
} // energy fxn
