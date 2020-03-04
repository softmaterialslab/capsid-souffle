#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include "functions.h"
#include "analysis.h"
#include "oligomer.h"
#include "subunit.h"
#include "bead.h"

using namespace std;

int analyze_output(int argc, char *argv[]) {
   string directory ="outfiles";               //append bin to home directory
   string coords = directory+"/ovito.lammpstrj";
   string nrg = directory+"/energy.out";
   string analysis = directory+"/analysis";
   int number_of_particles = 0;
   
   ifstream crds;                                           //open coordinates file
   crds.open(coords.c_str());
   if (!crds) {                                             //check to make sure file is there
      cerr << "ERR: FILE " << coords << " NOT OPENED. Check directory and/or filename.";
      exit(1);
   }
   cout << endl;
   cout << "Mass spectrum file will be saved to " << directory+"/analysis.ms" << endl;
   cout << "Pair-correlation file will be saved to " << directory+"/analysis.gr" << endl;
   cout << "R value will be saved to " << directory+"/analysis.rdat" << endl;
   cout << "Sub-oligomer structural information will be saved to " << directory+"/analysis.so" << endl;
   
   string dummy ;                                           //dummy string for text in file
   unsigned int number_of_beads;                            //how many beads in the file
   int number_of_subunits = 0;
   //int number_of_timesteps = 15920;                       //how many timesteps in the file
   int bead_index;                                          //index of bead
   int type;                                                //type of bead
   long double x, y, z;                                     // x y z coordinates
   long double box_size;
   
   vector<OLIGOMER> oligomers_list;                         //create vector to hold oligomers for mass spectrum analysis
   vector<OLIGOMER> sub_oligomers_list;                     //create vector to hold sub_oligomers for morphology analysis
   vector<BEAD> subunit_bead;                               //Create particles, named subunit_bead
   vector<SUBUNIT> protein;                                 //create subunits, named protein
   
   ofstream msdata( (directory+"/analysis.ms").c_str() );
   ofstream sodata( (directory+"/analysis.so").c_str() );
   ofstream grdata( (directory+"/analysis.gr").c_str() );
   ofstream rdata( (directory+"/analysis.rdat").c_str() );
   ofstream grinfo( (directory+"/info.gr").c_str() );
   
   int mstime = -1;                                         //parameter for ms_bin filling      
   vector<int> massbins(protein.size());
   vector<vector<int> > ms_bin;
   int sotime = -1;                                         //parameter for so_bin filling      
   vector<int> sobins(protein.size()*6);
   vector<vector<int> > so_bin;
   vector<int> pair_bins;
   
   int hiteqm = 0;
   bool done;
   int a = 0;                                               //index
   double delta_g = 0;
   int gr_size = 100;                                       //bin number for g(r)
   int beadType = 5;                                        // which type of bead to analyze with pcf?
   int num_beadType = 0;                                    //number of beads in the box of beadType
   double rho;                                              //density of beadType in the box
   
   compute_MD_trust_factor_R( hiteqm, done, directory );
   
   if (done){
      while ( crds >> dummy ){                              //while it hasn't reached the end of the line yet
         int timestep;
         crds >> dummy >> timestep >> dummy >> dummy >> dummy >> dummy ;      //reading the header of the timestep
         crds >> number_of_beads;
         crds >> dummy >> dummy >> dummy >> dummy >> box_size >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy ;
         if (a==0){                                                           //set up the simulation box environment
            if (number_of_beads % 41 == 0) number_of_particles = 41;          //determine how many beads in the subunit
            if (number_of_beads % 43 == 0) number_of_particles = 43;
            if (number_of_beads % 78 == 0) number_of_particles = 78;
            subunit_bead.resize(number_of_beads);
            number_of_subunits = number_of_beads / number_of_particles;
            protein.resize(number_of_subunits);  
            for (unsigned int b = 0; b < number_of_beads; b++){
               subunit_bead[b].bx.x = box_size*2;
               subunit_bead[b].bx.y = box_size*2;
               subunit_bead[b].bx.z = box_size*2;
            }
            delta_g = subunit_bead[0].bx.x / (2*gr_size);
         }
         
         int myindex;
         pair_bins.resize(gr_size);
         
         for (unsigned int b = 0; b < number_of_beads; b++){                  //read bead info for that timestep
            crds >> bead_index >> type >> x >> y >> z >> dummy >> dummy;
            myindex = b / (int)number_of_particles;
            subunit_bead[b].id = bead_index-1;
            
            subunit_bead[b].type = type;
            subunit_bead[b].pos.x = x;
            subunit_bead[b].pos.y = y;
            subunit_bead[b].pos.z = z;
            if (a==0){
               protein[myindex].itsB.push_back(&subunit_bead[b]);
               protein[myindex].id = myindex;
               if (subunit_bead[b].type == beadType) num_beadType += 1;       //count how many beads are of beadType
            }
         }
         
         if (a==0) { //determine rho
            rho = num_beadType / (subunit_bead[0].bx.x * subunit_bead[0].bx.x * subunit_bead[0].bx.x); //number density ******** (ASSUMPTION sigma = 1nm!) *******
         }         
         
         if (timestep >= hiteqm && (a % 1 == 0)) {
            cout << "Analyzing position file " << a << " with timestep " << timestep << endl;
            ms_bin.push_back(vector<int>() );
            so_bin.push_back(vector<int>() );
                                                                           //Generate pair_correlation bin
            for (unsigned int i = 0; i < number_of_beads; i++){
               for (unsigned int j = i + 1; j < number_of_beads; j++) {
                  if (subunit_bead[i].type == beadType && subunit_bead[j].type == beadType){
                     double distance = dist(&subunit_bead[i], &subunit_bead[j]).GetMagnitude();
                     if (distance < (box_size)) {
                        int gr_index = floor(distance/delta_g);
                        pair_bins[gr_index] += 1;
                     }
                  }
               }
            }
                                                                           // Save pair data to a file
            for (unsigned int j = 0; j < pair_bins.size(); j++) {
               grdata << pair_bins[j] << setw(15);                         //print gr bin data to file
               pair_bins[j] = 0;                                           //Clear gr bin
            }
            grdata << endl;
            
            
            int index = -1;
            for (unsigned int i = 0; i < protein.size(); i++) {            //Create oligomers for mass spectrum analysis
               unsigned int oldsize = 0;
               
               if (protein[i].itsO.size() == 0) {                          //if the unit isn't already counted...
                  
                  oligomers_list.push_back(OLIGOMER(VECTOR3D(0, 0, 0)));   //create an oligomer for the unit
                  index += 1;
                  oligomers_list[index].itsS.push_back(&protein[i]);       //add unit to oligomer
                  oligomers_list[index].id = index;
                  protein[i].itsO.push_back(oligomers_list[index]);        //add oligomer to unit
                  while (oldsize < oligomers_list[index].itsS.size()) {    //while the oligomer is still growing...
                     int n = oldsize;
                     oldsize = (int) oligomers_list[index].itsS.size();    //see how much the oligomer has grown
                     for (unsigned int j = n; j < oldsize; j++) {                   //loop over the growth from last round
                        int g = oligomers_list[index].itsS[j]->id;
                        for (unsigned int k = i + 1; k < protein.size(); k++) { //look for new growth
                           if (protein[k].itsO.size() == 0) {              //if it isn't in an oligomer yet...
                              for (unsigned int m = 0;
                                   m < protein[g].itsB.size(); m++) {        //check to see if it is in this oligomer
                                      for (unsigned int n = 0; n < protein[k].itsB.size(); n++) {
                                         if (dist(protein[g].itsB[m], protein[k].itsB[n]).GetMagnitude() < 2) {
                                            //  if (protein[k].itsB[n]->id == 0 && protein[g].itsB[m]->id == 0){
                                            oligomers_list[index].itsS.push_back(&protein[k]);    //if it is attached, add it
                                            protein[k].itsO.push_back(oligomers_list[index]);     //mark subunit as bonded
                                            goto finish;
                                            //    }
                                         } //if
                                      } //for n
                                   } //for m
                           } //if
                           finish:;
                        } //for k
                     } //for j
                  } //while
               } //if
            } //for i
            
            mstime += 1;
            ms_bin[mstime].resize(protein.size());
            for (unsigned int i = 0; i < oligomers_list.size(); i++) {
               if (oligomers_list[i].itsS.size() >= 1) {
                  ms_bin[mstime][(oligomers_list[i].itsS.size() - 1)] += 1;//fill mass bins
               }
            }

            
            index = -1;
            for (unsigned int i = 0; i < oligomers_list.size(); i++) {              //Loop over oligomers
               for (unsigned int j = 0; j < oligomers_list[i].itsS.size(); j++){    //Loop over subunits within the oligomer
                  for (unsigned int jj = 0; jj < oligomers_list[i].itsS[j]->itsB.size(); jj++){ //Loop over beads within the subunit
                     unsigned int oldsize = 0;
                     if (oligomers_list[i].itsS[j]->itsB[jj]->itsSO.size() == 0 && oligomers_list[i].itsS[j]->itsB[jj]->type == 0){ // if bead not counted yet
                        sub_oligomers_list.push_back(OLIGOMER(VECTOR3D(0,0,0)));    //create a sub_oligomer for the bead
                        index += 1;
                        sub_oligomers_list[index].itsB.push_back(oligomers_list[i].itsS[j]->itsB[jj]); //add bead to oligomer
                        sub_oligomers_list[index].id = index;
                        oligomers_list[i].itsS[j]->itsB[jj]->itsSO.push_back(sub_oligomers_list[index]); //add oligomer to bead
                        while (oldsize < sub_oligomers_list[index].itsB.size()) {   //while the oligomer is still growing
                           int n = oldsize;
                           oldsize = (int) sub_oligomers_list[index].itsB.size();   //see how much the oligomer has grown
                           for (unsigned int growth = n; growth < oldsize; growth++) {                      //loop over the growth from the last round
                              int g = sub_oligomers_list[index].itsB[growth]->id;
                              for (unsigned int k = 0; k < oligomers_list[i].itsS.size(); k++) {     //Look for new growth
                                 for (unsigned int kk = 0; kk < oligomers_list[i].itsS[k]->itsB.size(); kk++) {   //Loop over other beads
                                    if (oligomers_list[i].itsS[k]->itsB[kk]->itsSO.size() == 0 && oligomers_list[i].itsS[k]->itsB[kk]->type == 0) {//if it isn't in an oligomer yet
                                       int newg = (int) oligomers_list[i].itsS[k]->itsB[kk]->id;
                                       if (dist(&subunit_bead[g], &subunit_bead[newg]).GetMagnitude() < 2) { //check to see if it is in this oligomer
                                          sub_oligomers_list[index].itsB.push_back(oligomers_list[i].itsS[k]->itsB[kk]); //if it is attached, add it
                                          oligomers_list[i].itsS[k]->itsB[kk]->itsSO.push_back(sub_oligomers_list[index]); //mark bead as bonded
                                          goto finish2;
                                       }
                                    }
                                    finish2:;
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
            
            sotime += 1;
            so_bin[sotime].resize(protein.size()*2);
            for (unsigned int i = 0; i < sub_oligomers_list.size(); i++) {
               if (sub_oligomers_list[i].itsB.size() >= 1) {
                  so_bin[sotime][((sub_oligomers_list[i].itsB.size() - 1)) / 3] += 1;//fill so bins
               }
            }
            
            for (unsigned int j = 0; j < so_bin[sotime].size(); j++) {
               sodata << so_bin[sotime][j] << setw(15);                    //print mass bin data to file
            }
            sodata << endl;
            
            cout << "Printed mass data from step " << a << " to file" << endl;
            
            
            
            for (unsigned int i = 0; i < protein.size(); i++) {            // clear oligomer pointers from subunit
               protein[i].itsO.clear();
               for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
                  protein[i].itsB[ii]->itsSO.clear();
               }
            }
            
            oligomers_list.erase(oligomers_list.begin(),oligomers_list.end());               //erases oligomer objects
            sub_oligomers_list.erase(sub_oligomers_list.begin(),sub_oligomers_list.end());               //erases oligomer objects
            
            
            
            
         } //if at equilibrium
         a++;
      } //while file
      
      //Print some extra info for pair correlation file
      grinfo << "Number of Subunits: " << number_of_subunits << endl;
      grinfo << "PCF bin size: " << delta_g << endl;
      grinfo << "PCF bin number: " << gr_size << endl;
      grinfo << "Box size: " << box_size*2 << " " << box_size*2 << " " << box_size*2 << endl;
      grinfo << "Number density: " << rho << endl;
      grinfo << "Bead Type: " << beadType << endl;
   }//if (done)
   
   if (!done) {
      remove( (directory+"/info.gr").c_str() );
   }
   
   vector<double> ms_histogram(protein.size());
                                                                        //Make histogram for mass spectrum
   for (unsigned int i = 0; i < ms_bin.size(); i++) {
      for (unsigned int j = 0; j < ms_bin[i].size(); j++) {
         ms_histogram[j] += ms_bin[i][j] * (j + 1);                     //sum and weight the bin
      }
   }
   for (unsigned int i = 0; i < ms_histogram.size(); i++) {
      ms_histogram[i] = ms_histogram[i]/(ms_bin.size());                //average the bin
   }
   while (!ms_histogram.empty() && ms_histogram[ms_histogram.size() -1] == 0) {
      ms_histogram.pop_back();                                          //remove trailing 0s
   }
   for (unsigned int i = 0; i < ms_histogram.size(); i++) {
      msdata << ms_histogram[i] << setw(15);                            //save msdata to file
   }
   msdata << endl;
   
   
   
   
   return 0;
}