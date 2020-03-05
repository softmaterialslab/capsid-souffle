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
#include "oligomer.h"
#include <iterator>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

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
double compute_MD_trust_factor_R(int &hiteqm, bool &done, string directory) {
   string inPath = directory+"/energy.out";
   ifstream in(inPath.c_str(), ios::in);
   if (!in) {
      if (!in)
         cout << "ERR: FILE " << directory+"/energy.out" << " NOT OPENED. Check directory and/or filename." << endl;
      return 0;
   }
   string dummy;
   string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
   double col2_db, col7_db, col8_db, col5_db, col1_db;
   vector<double> ext, ke, pe, lj, time;
   in >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
   while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9) {
      try {  
         col1_db = boost::lexical_cast<double>(col1);
         col2_db = boost::lexical_cast<double>(col2); //use boost lexical_cast to look for extra headers from restart files
         col7_db = boost::lexical_cast<double>(col7);
         col8_db = boost::lexical_cast<double>(col8);
         col5_db = boost::lexical_cast<double>(col5);
      } catch (boost::bad_lexical_cast&) {
         cout << "Caught a restart header! This will artificially inflate R value." << endl;
         goto next;
      }
      in >> col10 >> col11;
      time.push_back(col1_db);
      ext.push_back(col7_db);
      ke.push_back(col2_db);
      pe.push_back(col8_db);
      lj.push_back(col5_db);
      next: ;
   }
   //determine when equilibrium is reached
   double x_sum, y_sum, xy_sum, x_m, y_m;
   double xx_sum, slope;
   int number_sections = 10;
   int section_size = floor(pe.size()/number_sections);
   vector<double> slope_vec;
   slope_vec.resize(number_sections);
   vector<double> stdev_vec;
   stdev_vec.resize(number_sections);
   cout << "Each energy section has " << section_size << " points." << endl;
   
   vector<RunningStat> RS;
   RS.resize(number_sections);
   double mean;
   double variance;
   double stdev;
   
   for (int i = 0; i < number_sections; i++){ //loop over sections
      x_sum = 0;
      y_sum = 0;
      xy_sum = 0;
      xx_sum = 0;
      for (int j = 0; j < section_size; j++){ //loop over points in sections
         x_sum += j;
         y_sum += pe[(i*section_size + j)];
         xy_sum += j * pe[(i*section_size + j)];
         xx_sum += j * j;

         RS[i].Push(pe[(i*section_size + j)]);
      }
      x_m = x_sum / double(section_size);
      y_m = y_sum / double(section_size);
      slope = ((xy_sum) - (double(section_size) * x_m * y_m) ) / (xx_sum - (double(section_size) * x_m * y_m) ); 
      //   cout << "For section " << i << " the slope is " << slope << endl;
      slope_vec[i] = slope;
      
      mean = RS[i].Mean();
      variance = RS[i].Variance();
      stdev = RS[i].StandardDeviation();
      stdev_vec[i] = stdev;
   }
   RunningStat rs;
   //find sections which have a slope +- 1e-6
   vector<bool> flat_vec;
   flat_vec.resize(number_sections);
   for (int i = 0; i < number_sections; i++) {
      //if (slope_vec[i] < 1e-6 && slope_vec[i] > -1e-6) {
      if (stdev_vec[i] < 0.02 && slope_vec[i] < 1e-4 && slope_vec[i] > -1e-4) {
         flat_vec[i] = true;
         //  cout << "section " << i << " is flat." << endl; 
      }
      else flat_vec[i] = false;
   }
   flat_vec[9] = true;
   //See how many of the final sections are at equilibrium
   int flat_sections = 0;
   for (int i = (number_sections - 1); i > -1; i--) {
      if (flat_vec[i] == true) {
         flat_sections += 1;
         hiteqm = i*section_size;
      }
      else break;
   }
   hiteqm = time[hiteqm];
   cout << "There are " << flat_sections << " ending consecutive flat sections" << endl;
   done = true;
   //Let the user know how long equilibrium has been reached (or if it needs restart)
   if (flat_sections > 1) {
      cout << 100.0*double(flat_sections)/double(number_sections) << "% of simulation is in equilibrium. Analyzing..." << endl;
      cout << "Equilibrium starts at timestep " << hiteqm << "." << endl;
   } else if (flat_sections ==1) {
      cout << "Simulation has not reached equilibrium! Analyzing last 10% anyway..." << endl;
     // done = false;
   } 
   //Compute R
   double ext_mean = 0;
   for (unsigned int i = 0; i < ext.size(); i++)
      ext_mean += ext[i];
   ext_mean = ext_mean / (ext.size() - 0 );
   double ke_mean = 0;
   for (unsigned int i = 0; i < ke.size(); i++)
      ke_mean += ke[i];
   ke_mean = ke_mean / ke.size();
   
   double ext_sd = 0;
   for (unsigned int i = 0; i < ext.size(); i++)
      ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
   ext_sd = ext_sd / (ext.size() - 0 );
   ext_sd = sqrt(ext_sd);
   
   double ke_sd = 0;
   for (unsigned int i = 0; i < ke.size(); i++)
      ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
   ke_sd = ke_sd / (ke.size() - 0 );
   ke_sd = sqrt(ke_sd);
   
   double R = ext_sd / ke_sd;
   //    if (world.rank() == 0)
   //    {
   ofstream out( (directory+"/analysis.rdat").c_str() );
   out << "Sample size " << ext.size() << endl;
   out << "Sd: ext, kinetic energy and R" << endl;
   out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
   //    }
   cout << endl << endl << "R is: " << R << endl << endl;
   
   return R;
}



vector<string> getFileNames(string directoryPath)
{
   namespace fs = boost::filesystem ;
   vector<string> names ;
   
   if ( fs::exists(directoryPath) )
   {
      fs::directory_iterator it(directoryPath) ;
      fs::directory_iterator end ;
      
      while ( it != end )
      {
         names.push_back(it->path().filename().string()) ;
         ++it ;
      }
   }
   
   return names ;
}




void filter(vector<string>& strings, string pattern)
{
   vector<string>::iterator pos = remove_if(strings.begin(), strings.end(), isRestart) ; 
   strings.erase(pos, strings.end()) ;
}



bool numeric_string_compare(const std::string& s1, const std::string& s2) {
   if (s1.length() < s2.length()) //You need this to sort integers in file name
      return true;
   if (s2.length() < s1.length())
      return false;
   else
      return (s1 < s2);
}



bool isRestart (string& s) {
   return s.find("restart") == std::string::npos ; 
}


