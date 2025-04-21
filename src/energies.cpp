#include<iostream>
#include <fstream>
#include <iomanip>
#include"energies.h"
#include "bead.h"
#include "subunit.h"
#include "edge.h"
#include "face.h"
#include "functions.h"

using namespace std;

long double particle_kinetic_energy(vector <BEAD> &subunit_bead){ //part of thermostat
   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      subunit_bead[i].update_kinetic_energy();
   
   long double kinetic_energy = 0.0;
   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      kinetic_energy += subunit_bead[i].ke;
   
   return kinetic_energy;
}
 
void update_LJ_ES_energies_simplified(vector<BEAD>& subunit_bead, double ecut, vector<vector<int> > lj_a, double elj_att, double lb, double ni, double qs, double ecut_el, double kappa){
   
   double ecut3 = ecut * ecut * ecut;
   double ecut6 = ecut3 * ecut3;
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
            
         if (electrostatic && r < ecut_el && ecut_el != 0){//electrostatic energy
            subunit_bead[i].ce += ( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) ) / (r) -
                                  ( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*ecut_el) ) / (ecut_el)) );
         } else if (electrostatic && ecut_el == 0){
            subunit_bead[i].ce += ( (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) ) / (r) );
         } else{
             subunit_bead[i].ce += 0;
         }
         
         double sig1 = subunit_bead[i].sigma;
         double sig2 = subunit_bead[j].sigma;
         double shc = (sig1 + sig2)/2;
         double sigma3 = shc * shc * shc;
         double sigma6 = sigma3 * sigma3;
         double r3 = r * r * r;
         double r6 = r3 * r3;
         if (r >= ecut && r >= (1.12246205 * shc)) continue;  //for (6,12) potential
         //if (r >= ecut && r >= (1 * shc)) continue;  //for (6,9) potential
         
         if (lj) {
            bool lj_attractive = false;
            for (unsigned int k = 0; k < lj_a[0].size(); k++){
               if (subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]){
                  lj_attractive = true;
                  //elj_att = lj_a[3][k];
                  break;
               }
            }
            if (r < (1.12246205*shc) && lj_attractive == false){
            //if (r < (1*shc) && lj_attractive == false){           //for (6,9) potential, Repulsive
               double elj = 1; //repulsive lj fixed @ 1. Also fixed in forces.cpp
               //subunit_bead[i].ne += ((elj * (sigma6 / r6) * (2*(sigma3 / r3) - 3)) + elj);   //for (6,9) potential
               subunit_bead[i].ne += ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
            } else if ( r < ecut && lj_attractive == true ){            //Attractive
               //subunit_bead[i].ne += (((elj_att * (sigma6 / r6) * (2*(sigma3 / r3) - 3)) - (elj_att * ((sigma6 / ecut6) * (2*(sigma3 / ecut3) - 3))))); //for (6,9) potential
               subunit_bead[i].ne += (((4 * elj_att * (sigma6 / r6) * ((sigma6 / r6) - 1)) - (4 * elj_att * ((sigma6 / ecut6) * ((sigma6 / ecut6) - 1)))));
            } else {
               subunit_bead[i].ne += 0;
            }
         } //if (lj)
      } // for j
   } // for i
} // energy fxn


void update_LJ_ES_energies_bigbeads(vector<BEAD>& big_bead, vector<BEAD>& subunit_bead, double ecut, double lb, double ni, double qs, double ecut_el, double kappa, double shc) {
   VECTOR3D box = subunit_bead[0].bx;
   VECTOR3D hbox = subunit_bead[0].hbx;

   for (unsigned int i = 0; i < big_bead.size(); i++){
      for (unsigned int j = 0; j < subunit_bead.size(); j++){

         bool electrostatic = (subunit_bead[j].q != 0);
         long double r = 0.0;
         double ce = 0.0;
         double ne = 0.0;

         VECTOR3D r_vec;
         r_vec = big_bead[i].pos - subunit_bead[j].pos;
         if (r_vec.x > hbox.x) r_vec.x -= box.x;
         else if (r_vec.x < -hbox.x) r_vec.x += box.x;
         if (r_vec.y > hbox.y) r_vec.y -= box.y;
         else if (r_vec.y < -hbox.y) r_vec.y += box.y;
         if (r_vec.z > hbox.z) r_vec.z -= box.z;
         else if (r_vec.z < -hbox.z) r_vec.z += box.z;
         r = r_vec.GetMagnitude();

         double sig1 = big_bead[i].sigma;
         double sig2 = subunit_bead[j].sigma;
         double sigavg = (sig1 + sig2)/2;
         double delta = sigavg - shc;
         double shc6 = pow(shc, 6);
         double r6 = pow((r-delta), 6);

         double yukawa_scale = exp(kappa * sigavg)/((1+kappa*sig1/2)*(1+kappa*sig2/2));
         if (electrostatic && r < ecut_el) ce = ( (big_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) * yukawa_scale ) / (r)); else ce = 0;
         big_bead[i].ce += ce;
         subunit_bead[j].ce += ce;

         if (r < delta + 1.12246205*shc) ne = ((4 * (shc6 / r6) * ((shc6 / r6) - 1)) + 1); else ne = 0;
         big_bead[i].ne += ne;
         subunit_bead[j].ne += ne;

      }
   }

}