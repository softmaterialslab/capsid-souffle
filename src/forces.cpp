#include<iostream>
#include <fstream>
#include <iomanip>
#include"forces.h"
#include "bead.h"
#include "subunit.h"
#include "edge.h"
#include "face.h"
#include "functions.h"

using namespace std;

void forceCalculation(vector<SUBUNIT> &protein, double lb, double ni, double qs, vector<BEAD> &subunit_bead,
                      double ecut, double ks, double kb, vector<vector<int> > lj_a, double ecut_el, 
                      double kappa, double elj_att, bool updatePairlist, double NListCutoff){
   
   //Common MPI Message objects
   vector<VECTOR3D> forvec(sizFVec * protein[0].itsB.size(), VECTOR3D(0, 0, 0));
   vector<VECTOR3D> forvecGather(subunit_bead.size(), VECTOR3D(0, 0, 0));
   
   //global variables
   unsigned int i, k;
   VECTOR3D box = subunit_bead[0].bx;
   VECTOR3D hbox = subunit_bead[0].hbx;
   double ecut_el2 = ecut_el * ecut_el;
   double ecut2 = ecut * ecut;
   double replj2 = 1.12246205 * 1.12246205;

   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*                               PARALLEL LOOP                                                          */
   //////////////////////////////////////////////////////////////////////////////////////////////////////////

   #pragma omp parallel for schedule(dynamic) default(shared) private(i, k)
   for (i = lowerBound; i <= upperBound; i++){

      VECTOR3D r_vec = VECTOR3D(0, 0, 0);
      long double r2 = 0.0;
      VECTOR3D sForce = VECTOR3D(0, 0, 0);
      VECTOR3D bForce = VECTOR3D(0, 0, 0);



      for (unsigned int ii = 0; ii < protein[i].itsE.size(); ii++){
         protein[i].itsE[ii]->update_length();
      }
      for (unsigned int ii = 0; ii < protein[i].itsF.size(); ii++){
         protein[i].itsF[ii]->update_area_normal();
      }

      for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){
         protein[i].itsB[ii]->update_stretching_force(ks);
         protein[i].itsB[ii]->bforce = VECTOR3D(0, 0, 0);         //resetting bending force here
      }
      for (unsigned int mm = 0; mm < protein[i].itsE.size(); mm++){
         if (protein[i].itsE[mm]->type != 0)                      //if it is a bending edge...
            protein[i].itsE[mm]->update_bending_forces(kb);
      }

     //////////////////////////////////////////////////////////////////////////////////////////////////////////
     /*                               INTERMOLECULAR FORCES                                                  */
     //////////////////////////////////////////////////////////////////////////////////////////////////////////

      if (updatePairlist == true){
         update_pairlist(i, protein, NListCutoff, box, hbox);     //Updating the pairlist
      }

      for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){

         VECTOR3D eForce = VECTOR3D(0, 0, 0);
         VECTOR3D ljForce = VECTOR3D(0, 0, 0);

         for (unsigned int jj = 0; jj < protein[i].itsB[ii]->itsN.size(); jj++) {

            int i_bead = protein[i].itsB[ii]->id;
            int j_bead = protein[i].itsB[ii]->itsN[jj];

            if (j_bead == -1) break; // checking for -1 because this is the "empty" value, an index no bead can have.

            bool electrostatic = (i_bead != j_bead && subunit_bead[i_bead].q != 0 && subunit_bead[j_bead].q != 0);
            bool lj =(protein[i].id != subunit_bead[j_bead].itsS[0]->id);

            long double r = 0.0;

            if (electrostatic || lj) {
               r_vec = subunit_bead[i_bead].pos - subunit_bead[j_bead].pos;
               if (r_vec.x > hbox.x) r_vec.x -= box.x;
               else if (r_vec.x < -hbox.x) r_vec.x += box.x;
               if (r_vec.y > hbox.y) r_vec.y -= box.y;
               else if (r_vec.y < -hbox.y) r_vec.y += box.y;
               if (r_vec.z > hbox.z) r_vec.z -= box.z;
               else if (r_vec.z < -hbox.z) r_vec.z += box.z;
               r2 = r_vec.GetMagnitudeSquared();
            }

            if (electrostatic && r2 < (ecut_el2) && ecut_el != 0){
               r = sqrt(r2);
               eForce += r_vec ^ ((subunit_bead[i_bead].q * subunit_bead[j_bead].q * lb * exp(-kappa * r)
               / (r2)) * (kappa + 1 / r));
            } else if (electrostatic && ecut_el == 0){
               r = sqrt(r2);
               eForce += r_vec ^ ((subunit_bead[i_bead].q * subunit_bead[j_bead].q * lb * exp(-kappa * r)
               / (r2)) * (kappa + 1 / r));
            }

            double sig1 = subunit_bead[i_bead].sigma;
            double sig2 = subunit_bead[j_bead].sigma;
            double shc = (sig1 +sig2)/2;

            if (r2 >= ecut2 && r2 >= (replj2 * shc * shc))
               continue;

            double test;

            if (lj){
               bool lj_attractive = false;
               for (k = 0; k < lj_a[0].size(); k++){
                  if (subunit_bead[i_bead].type == lj_a[1][k] && subunit_bead[j_bead].type == lj_a[2][k]){
                     lj_attractive = true;
                     //elj_att = lj_a[3][k];
                     break;
                  }
               }
               if (lj_attractive == true && r2 < (ecut2)){
                  if (subunit_bead[j_bead].type == 7) test = 50;
                  else test = elj_att;
                  double r6 = r2 * r2 * r2;
                  double sigma6 = sig1 * sig1 * sig1 * sig1 * sig1 * sig1;
                  ljForce += (r_vec ^ (48 * test * ((sigma6 / r6) * ((sigma6 / r6) - 0.5) ) * (1 / (r2))));
               }
               else if (lj_attractive == false && r2 < (replj2 * shc * shc)){
                  double r6 = r2 * r2 * r2;
                  double sigma6 = shc * shc * shc * shc * shc * shc;
                  double elj = 1;//fixed @ 1. Also fixed in energies.cpp
                  ljForce += (r_vec ^ (48 * elj * ((sigma6 / r6) * ((sigma6 / r6) - 0.5) ) * (1 / (r2))));
               }
            } //if lj
         } // for jj
         sForce = protein[i].itsB[ii]->sforce;
         bForce = protein[i].itsB[ii]->bforce;
         VECTOR3D eljforces = eForce + ljForce;
         forvec[(i - lowerBound)*protein[i].itsB.size()+ii] = sForce + bForce + eljforces;
      } // for ii
   } //for i

   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*                               BROADCASTING & ACCUMULATION                                            */
   //////////////////////////////////////////////////////////////////////////////////////////////////////////

   //forvec broadcasting using all gather = gather + broadcast
   if (world.size() > 1) {
      all_gather(world, &forvec[0], forvec.size(), forvecGather);
   } else {
      for (i = lowerBound; i <= upperBound; i++){
         for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
            forvecGather[i*protein[i].itsB.size()+ ii] = forvec[(i - lowerBound)*protein[i].itsB.size()+ii];
         }
      }
   }

   //Total force accumulation

   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      subunit_bead[i].tforce = forvecGather[i];

   forvec.clear();
   forvecGather.clear();

} //void fxn


void forceCalculation_bigbead(vector<BEAD>& big_bead, vector<BEAD>& subunit_bead, double ecut, double lb, double ni, double qs, double ecut_el, double kappa, double shc) {

   VECTOR3D box = subunit_bead[0].bx;
   VECTOR3D hbox = subunit_bead[0].hbx;

   //reset elec and lj force
   for (unsigned int i = 0; i < big_bead.size(); i++) {
      big_bead[i].eforce = VECTOR3D(0, 0, 0);
      big_bead[i].ljforce = VECTOR3D(0, 0, 0);
   }

   for (unsigned int i = 0; i < subunit_bead.size(); i++) {
      subunit_bead[i].eforce = VECTOR3D(0, 0, 0);
      subunit_bead[i].ljforce = VECTOR3D(0, 0, 0);
   }

   for (unsigned int i = 0; i < big_bead.size(); i++){
      for (unsigned int j = 0; j < subunit_bead.size(); j++){

         bool electrostatic = (subunit_bead[j].q != 0);
         long double r = 0.0;
         long double r2 = 0.0;
         VECTOR3D eForce = VECTOR3D(0, 0, 0);
         VECTOR3D ljForce = VECTOR3D(0, 0, 0);

         VECTOR3D r_vec;
         r_vec = big_bead[i].pos - subunit_bead[j].pos;
         if (r_vec.x > hbox.x) r_vec.x -= box.x;
         else if (r_vec.x < -hbox.x) r_vec.x += box.x;
         if (r_vec.y > hbox.y) r_vec.y -= box.y;
         else if (r_vec.y < -hbox.y) r_vec.y += box.y;
         if (r_vec.z > hbox.z) r_vec.z -= box.z;
         else if (r_vec.z < -hbox.z) r_vec.z += box.z;
         r = r_vec.GetMagnitude();
         r2 = r_vec.GetMagnitudeSquared();

         double sig1 = big_bead[i].sigma;
         double sig2 = subunit_bead[j].sigma;
         double sigavg = (sig1 + sig2)/2;
         double delta = sigavg - shc;
         double shc6 = pow(shc, 6);
         double r6 = pow((r-delta), 6);

         double yukawa_scale = exp(kappa * sigavg)/((1+kappa*sig1/2)*(1+kappa*sig2/2));
         r_vec ^ ((big_bead[i].q * subunit_bead[j].q * lb * exp(-kappa * r)/ (r2)) * (kappa + 1 / r));
         if (electrostatic) eForce = r_vec ^ ((big_bead[i].q * subunit_bead[j].q * lb * exp(-kappa * r)/ (r2)) * (kappa + 1 / r) * yukawa_scale); else eForce = VECTOR3D(0, 0, 0);
         big_bead[i].eforce = big_bead[i].eforce + eForce;
         subunit_bead[j].eforce = subunit_bead[j].eforce - eForce;

         if (r2 < (delta + 1.12246205*shc)*(delta + 1.12246205*shc)) ljForce =(r_vec ^ (48 * ((shc6 / r6) * ((shc6 / r6) - 0.5) ) * (1 / r)*(1 / (r-delta)) )); else ljForce = VECTOR3D(0,0,0);
         big_bead[i].ljforce = big_bead[i].ljforce + ljForce;
         subunit_bead[j].ljforce = subunit_bead[j].ljforce - ljForce;

      }
   }

   for (unsigned int i = 0; i < big_bead.size(); i++)
      big_bead[i].tforce = big_bead[i].ljforce + big_bead[i].eforce;

   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      subunit_bead[i].tforce = subunit_bead[i].tforce + subunit_bead[i].ljforce + subunit_bead[i].eforce;


}