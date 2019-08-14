#include<iostream>
#include <fstream>
#include <iomanip>
#include"forces.h"
#include "bead.h"
#include "LJpair.h"
#include "subunit.h"
#include "edge.h"
#include "face.h"
#include "functions.h"


using namespace std;

void forceCalculation(vector<SUBUNIT> &protein, double lb, double ni, double qs, vector<BEAD> &subunit_bead,
                      vector<PAIR> &lj_pairlist, double ecut, double ks, double bondlength, double kb,
                      vector<vector<int> > lj_a, double ecut_el, double kappa, double elj_att, bool updatePairlist, double NListCutoff) {
   
   //ofstream forces("outfiles/forces.out", ios::app);
   
   //Common MPI Message objects
   vector<VECTOR3D> forvec(sizFVec, VECTOR3D(0, 0, 0));
   vector<VECTOR3D> forvecGather(subunit_bead.size() + extraElements, VECTOR3D(0, 0, 0));
   
   //global variables
   unsigned int i, j, k;
   VECTOR3D box = subunit_bead[0].bx;
   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*									INTRA MOLECULAR FORCES												*/
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   //#pragma omp parallel for schedule(dynamic) default(shared) private(i)
   for (unsigned int i = 0; i < protein.size(); i++) {
      for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
         protein[i].itsB[ii]->update_stretching_force(ks, bondlength);
         protein[i].itsB[ii]->bforce = VECTOR3D(0, 0, 0);        //resetting bending force here
      }
      
      for (unsigned int mm = 0; mm < protein[i].itsE.size(); mm++) {
         if (protein[i].itsE[mm]->type != 0)                    //if it is a bending edge...
            protein[i].itsE[mm]->update_bending_forces(kb);
      }
   }
   
   
   //LJ forces calculation and update ES forces between subunits
   #pragma omp parallel for schedule(dynamic) default(shared) private(i, j, k)
   for (i = lowerBound; i <= upperBound; i++) {
      
      VECTOR3D r_vec = VECTOR3D(0, 0, 0);
      long double r2 = 0.0;
      double hbox = box.x / 2;
      
      if (updatePairlist == true) {
         fill(subunit_bead[i].itsN.begin(), subunit_bead[i].itsN.end(), -1);  //clear the pairlist to -1 (a number that cannot be bead index)
         unsigned int test = 0;
         for (unsigned int bead_j = 0; bead_j < subunit_bead.size(); bead_j++) {
            bool electrostatic = (i != bead_j && subunit_bead[i].q != 0 && subunit_bead[bead_j].q != 0);
            bool lj = subunit_bead[i].itsS[0]->id != subunit_bead[bead_j].itsS[0]->id;
            if (electrostatic || lj) {
               r_vec = subunit_bead[i].pos - subunit_bead[bead_j].pos;
               if (r_vec.x > hbox) r_vec.x -= box.x;
               else if (r_vec.x < -hbox) r_vec.x += box.x;
               if (r_vec.y > hbox) r_vec.y -= box.y;
               else if (r_vec.y < -hbox) r_vec.y += box.y;
               if (r_vec.z > hbox) r_vec.z -= box.z;
               else if (r_vec.z < -hbox) r_vec.z += box.z;
               r2 = r_vec.GetMagnitudeSquared();
               if ( r2 < (NListCutoff * NListCutoff) ) {
                  subunit_bead[i].itsN[test] = subunit_bead[bead_j].id;
                  test += 1;
               }
            } //if lj||es
         }// for j
         if(test > subunit_bead[0].itsN.size()) cout << endl << "ERROR! Neighborlist outgrew allocated vector size!" << endl;
      } //if
      
      VECTOR3D eForce = VECTOR3D(0, 0, 0);
      VECTOR3D ljForce = VECTOR3D(0, 0, 0);
      
    //  for (k = 0; k < subunit_bead[i].itsN.size(); k++) {
      for (unsigned int m = 0; m < subunit_bead[i].itsN.size(); m++){
         
         if (subunit_bead[i].itsN[m] == -1) break; // checking for -1 because this is the "empty" value, an index no bead can have.
         
            j = subunit_bead[i].itsN[m];
            
            //Add electrostatic cut offs here.
            bool electrostatic = (i != j && subunit_bead[i].q != 0 && subunit_bead[j].q != 0);
            bool lj = subunit_bead[i].itsS[0]->id != subunit_bead[j].itsS[0]->id;
            
            
            long double r = 0.0;
            
            
            if (electrostatic || lj) {
               
               r_vec = subunit_bead[i].pos - subunit_bead[j].pos;
               if (r_vec.x > hbox) r_vec.x -= box.x;
               else if (r_vec.x < -hbox) r_vec.x += box.x;
               if (r_vec.y > hbox) r_vec.y -= box.y;
               else if (r_vec.y < -hbox) r_vec.y += box.y;
               if (r_vec.z > hbox) r_vec.z -= box.z;
               else if (r_vec.z < -hbox) r_vec.z += box.z;
               r2 = r_vec.GetMagnitudeSquared();
               
            }
            
            if (electrostatic && r2 < (ecut_el * ecut_el) && ecut_el != 0) {
               r = sqrt(r2);
               eForce += r_vec ^ ((subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa * r)
               / (r2)) * (kappa + 1 / r));
            } else if (electrostatic && ecut_el == 0) {
               r = sqrt(r2);
               eForce += r_vec ^ ((subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa * r)
               / (r2)) * (kappa + 1 / r));
            }
            
            double sig1 = subunit_bead[i].sigma;
            double sig2 = subunit_bead[j].sigma;
            double shc = (sig1 +sig2)/2;
            
            if (r2 >= ecut*ecut && r2 >= ((1.12246205 * shc) * (1.12246205 * shc)))
               continue;
            
            if (lj) {
               
               bool lj_attractive = false;
               for (k = 0; k < lj_a[0].size(); k++) {
                  if (subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]) {
                     lj_attractive = true;
                  }
               }
               
              
               if (lj_attractive == true && r2 < (ecut*ecut)) {
                  double r6 = r2 * r2 * r2;
                  double sigma6 = sig1 * sig1 * sig1 * sig1 * sig1 * sig1;
                  ljForce += (r_vec ^ (48 * elj_att * ((sigma6 / r6) * ((sigma6 / r6) - 0.5) ) * (1 / (r2))));
               }
               else if (lj_attractive == false && r2 < ((1.12246205 * shc) * (1.12246205 * shc))) {
                  double r6 = r2 * r2 * r2;
                  double sigma6 = shc * shc * shc * shc * shc * shc;
                  double elj = 1;//subunit_bead[j].epsilon;
                  ljForce += (r_vec ^ (48 * elj * ((sigma6 / r6) * ((sigma6 / r6) - 0.5) ) * (1 / (r2))));
               } 
            } //if lj
        } //for m (j)
      forvec[i - lowerBound] = eForce + ljForce;
      } //for i
      
   //forvec broadcasting using all gather = gather + broadcast
   if (world.size() > 1) {
      all_gather(world, &forvec[0], forvec.size(), forvecGather);
   } else {
      for (i = lowerBound; i <= upperBound; i++)
         forvecGather[i] = forvec[i - lowerBound];
   }
      
   //cout << lowerBound << " " << upperBound << " "<< subunit_bead.size()<<endl;
   //Total force accumulation
   //for (i = 0; i < subunit_bead.size(); i++)
   //    subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + lj[i];
   
   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + forvecGather[i];
   
   
   forvec.clear();
   forvecGather.clear();
      
      
} //void fxn



 
 
 
