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






 

void update_ES_forces_simplified(vector<BEAD>& subunit_bead, double lb, double ni, double qs)
 {
	 
	 for (unsigned int i=0; i < subunit_bead.size(); i++)
	 {
		 subunit_bead[i].eforce = VECTOR3D(0,0,0);
	 }
	 
	 for (unsigned int i = 0; i < subunit_bead.size(); i++)
	 {
		 for (unsigned int j = 0; j < subunit_bead.size(); j++)
		 {
			 if (subunit_bead[i].id != subunit_bead[j].id)
			 {
				 double kappa = sqrt (8 * 3.1416 * ni * lb * qs*qs);
				 VECTOR3D r_vec = dist( & subunit_bead[i] , & subunit_bead[j] );
				 long double r = r_vec.GetMagnitude();
				 VECTOR3D ff = r_vec ^ ( ( subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa * r )
				 / (r * r) ) * (kappa + 1/r) );
				 
				 subunit_bead[i].eforce += ff;
			 }
		 }
	 }
	 
 }
 
 
 
 
 
 
 
 
 
 void update_LJ_forces_simplified(vector<BEAD>& subunit_bead, double rcut, vector<vector<int> > lj_a){
	 for (unsigned int i = 0; i < subunit_bead.size(); i++){
		 subunit_bead[i].ljforce = VECTOR3D(0,0,0);
	 }
	 
	 double shc = 1.2;
	 double shc6 = shc * shc * shc * shc * shc * shc;
	 double shc12 = shc6 * shc6;
	 
	 for (unsigned int i = 0; i < subunit_bead.size(); i++){
		 for (unsigned int j = 0; j < subunit_bead.size(); j++){
			 if (subunit_bead[i].itsS[0]->id != subunit_bead[j].itsS[0]->id){
				 VECTOR3D r_vec = dist( &subunit_bead[i] , &subunit_bead[j] );
				 long double r = r_vec.GetMagnitude();
				 //double r2 = (r*r);
				 double r6 ;
				 double r12;
				 
				 double sig1 = subunit_bead[i].sigma;
				 double sig2 = subunit_bead[j].sigma;
				 double del = (sig1+sig2)/2 - shc;
				 bool lj_attractive = false;
				 for (unsigned int k = 0; k < lj_a[0].size(); k++){
					 if ( subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]  ){
						 lj_attractive = true;
					 }
				 }
				 if ( r < (del+1.12246205*shc) && lj_attractive == false){							//Repulsive
					 double elj = 1;//subunit_bead[j].epsilon;
					 r6 = (r-del) * (r-del) * (r-del) * (r-del) * (r-del) * (r-del);
					 r12 = r6 * r6;
					 subunit_bead[i].ljforce += (r_vec ^ (48 * elj * ((shc12 / r12) - 0.5 * (shc6 / r6)) * (1 / (r*(r-del))) ));
				 } else if ( r < ((del+1.12246205*shc)*rcut) && lj_attractive == true ){			//Attractive
					 r6 = (r-del) * (r-del) * (r-del) * (r-del) * (r-del) * (r-del);
					 r12 = r6 * r6;
					 double elj = 2;//subunit_bead[j].epsilon;
					 subunit_bead[i].ljforce += (r_vec ^ (48 * elj * ((shc12 / r12) - 0.5 * (shc6 / r6)) * (1 / (r*(r-del)))));
				 } else {
					 subunit_bead[i].ljforce += 0;
				 }
			 }
		 }
	 }
 }
 
 
 
 
