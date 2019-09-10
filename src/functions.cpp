//
// Created by lauren on 2/1/18.
//

#include<iostream>
#include <fstream>
#include <iomanip>
#include"functions.h"
#include "bead.h"
#include "subunit.h"
#include "edge.h"
#include "face.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

using namespace std;

void ProgressBar(double fraction_completed){
   int val = (int) (fraction_completed * 100);
   int lpad = (int) (fraction_completed * PBWIDTH);
   int rpad = PBWIDTH - lpad;
   printf("\r%3d%% |%.*s%*s|", val, lpad, PBSTR, rpad, "");
   fflush(stdout);
}

//part of thermostat
void update_chain_xi(unsigned int j, vector<THERMOSTAT> &bath, double dt, long double ke){ 
    if (bath[j].Q == 0)
        return;
// 	 double expfactor = exp(-0.25 * dt * bath[j + 1].xi);
// 	 
// 	 if (j != 0) {
//         bath[j].xi = bath[j].xi * expfactor * expfactor + 
//         0.5 * dt * (1.0 / bath[j].Q) * (bath[j - 1].Q * bath[j - 1].xi * bath[j - 1].xi - 
//         bath[j].dof * bath[j].kB * bath[j].T) * expfactor;
//     } else {
//         bath[j].xi = bath[j].xi * expfactor * expfactor + 
//         0.5 * dt * (1.0 / bath[j].Q) * (2 * ke - bath[j].dof * bath[j].kB * bath[j].T) * expfactor;
//     }
	 
    if (j != 0) {
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) + 0.5 * dt * (1.0 / bath[j].Q) *
                                                                    (bath[j - 1].Q * bath[j - 1].xi * bath[j - 1].xi -
                                                                     bath[j].dof * bath[j].kB * bath[j].T) *
                                                                    exp(-0.25 * dt * bath[j + 1].xi);
    } else {
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) +
                     0.5 * dt * (1.0 / bath[j].Q) * (2 * ke - bath[j].dof * bath[j].kB * bath[j].T) *
                     exp(-0.25 * dt * bath[j + 1].xi);
    }
    return;
}

void dress_up(vector<EDGE> &subunit_edge, vector<FACE> &subunit_face) {
   for (unsigned int i = 0; i < subunit_edge.size(); i++) {
      subunit_edge[i].update_length();
    }
   for (unsigned int i = 0; i < subunit_face.size(); i++) {
      subunit_face[i].update_area_normal();
   }
}

void update_pairlist(unsigned int i, vector<SUBUNIT> &protein, double NListCutoff, VECTOR3D box, VECTOR3D hbox) {
   VECTOR3D r_vec = VECTOR3D(0, 0, 0);
   long double r2 = 0.0;
   double NListCutoff2 = NListCutoff * NListCutoff;
   std::vector<unsigned int> NLindex( protein[i].itsB.size() );
   
   for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){
      fill(protein[i].itsB[ii]->itsN.begin(), protein[i].itsB[ii]->itsN.end(), -1); //clear the pairlist to -1 (a number that cannot be bead index)
      NLindex[ii] = 0;
   }
   
   for (unsigned int j = 0; j < protein.size(); j++) { // loop over subunits for ~(N/40)^2 scaling
      if (i != j){
         r_vec = protein[i].centerBead->pos - protein[j].centerBead->pos;  //Check center bead (determined in generate lattice)
         if (r_vec.x > hbox.x) r_vec.x -= box.x;
         else if (r_vec.x < -hbox.x) r_vec.x += box.x;
         if (r_vec.y > hbox.y) r_vec.y -= box.y;
         else if (r_vec.y < -hbox.y) r_vec.y += box.y;
         if (r_vec.z > hbox.z) r_vec.z -= box.z;
         else if (r_vec.z < -hbox.z) r_vec.z += box.z;
         r2 = r_vec.GetMagnitudeSquared();
      } else r2 = 0; //(i == j)
      if (r2 < (NListCutoff2 * 2.25) ) { //if the subunits are close together...
         for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){
            for (unsigned int jj = 0; jj < protein[j].itsB.size(); jj++) {
               bool electrostatic = (protein[i].itsB[ii]->id != protein[j].itsB[jj]->id && protein[i].itsB[ii]->q != 0 && protein[j].itsB[jj]->q != 0);
               bool lj = protein[i].id != protein[j].id;
               if (electrostatic || lj) {
                  if ( ii != protein[i].centerBead->id && jj != protein[j].centerBead->id) { //don't recalculate dist. center-to-center
                     r_vec = protein[i].itsB[ii]->pos - protein[j].itsB[jj]->pos;
                     if (r_vec.x > hbox.x) r_vec.x -= box.x;
                     else if (r_vec.x < -hbox.x) r_vec.x += box.x;
                     if (r_vec.y > hbox.y) r_vec.y -= box.y;
                     else if (r_vec.y < -hbox.y) r_vec.y += box.y;
                     if (r_vec.z > hbox.z) r_vec.z -= box.z;
                     else if (r_vec.z < -hbox.z) r_vec.z += box.z;
                     r2 = r_vec.GetMagnitudeSquared();
                  }
                  if ( r2 < (NListCutoff2) ) {
                     protein[i].itsB[ii]->itsN[NLindex[ii]] = protein[j].itsB[jj]->id; //add to pairlist
                     NLindex[ii] += 1;
                  }
               } //if lj||es
            } //for jj
            if(NLindex[ii] > protein[i].itsB[0]->itsN.size()) cout << endl << "ERROR! Neighborlist outgrew allocated vector size!" << endl;
         } //for ii
      } //if r2
   }// for j
} // update pairlist fxn

// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm) {
   string inPath = "sub_beads.traj.out";
   ifstream in(inPath.c_str(), ios::in);
   if (!in) {
      if (!in)
         if (world.rank() == 0)
            cout << "File could not be opened" << endl;
      return 0;
   }
   string dummy;
   double col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
   vector<double> ext, ke, pe;
   in >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
   while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11) {
   ext.push_back(col7);
   ke.push_back(col2);
   pe.push_back(col8);
   }

   double ext_mean = 0;
   for (unsigned int i = 0; i < ext.size(); i++)
      ext_mean += ext[i];
   ext_mean = ext_mean / ext.size();
   double ke_mean = 0;
   for (unsigned int i = 0; i < ke.size(); i++)
      ke_mean += ke[i];
   ke_mean = ke_mean / ke.size();

   double ext_sd = 0;
   for (unsigned int i = 0; i < ext.size(); i++)
      ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
   ext_sd = ext_sd / ext.size();
   ext_sd = sqrt(ext_sd);

   double ke_sd = 0;
   for (unsigned int i = 0; i < ke.size(); i++)
      ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
   ke_sd = ke_sd / ke.size();
   ke_sd = sqrt(ke_sd);

   double R = ext_sd / ke_sd;
//    if (world.rank() == 0)
//    {
   string outPath = "R.dat";
   ofstream out(outPath.c_str());
   out << "Sample size " << ext.size() << endl;
   out << "Sd: ext, kinetic energy and R" << endl;
   out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
//    }
   cout << endl << endl << "R is: " << R;
   return R;
}

