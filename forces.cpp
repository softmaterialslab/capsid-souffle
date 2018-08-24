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

void update_ES_forces_pairlist(vector<BEAD>& subunit_bead, double lb, double ni, double qs, vector<PAIR>& lj_pairlist){
	for (int i=0; i<subunit_bead.size(); i++){
		subunit_bead[i].eforce = VECTOR3D(0,0,0);
	}
	
	for (int i=0; i<lj_pairlist.size(); i++){
		VECTOR3D r_vec = dist( lj_pairlist[i].itsB[0] , lj_pairlist[i].itsB[1] );
		long double r = r_vec.GetMagnitude();
		double kappa = 8 * 3.1416 * ni * lb * qs*qs;
		VECTOR3D ff = r_vec ^ ( ( lj_pairlist[i].itsB[0]->q * lj_pairlist[i].itsB[1]->q * lb * exp(-kappa * r )
		/ (r * r) ) * (kappa + 1/r) );
		
		lj_pairlist[i].itsB[0]->eforce += ff;
		lj_pairlist[i].itsB[1]->eforce -= ff;
	}
	
}

void update_ES_forces_intra(vector<SUBUNIT>& protein, double lb, double ni, double qs){
	for (int i = 0; i < protein.size(); i++){
		for (int ii = 0; ii < protein[i].itsB.size(); ii++)							//intramolecular forces loop
		{
			for (int kk = ii + 1; kk < protein[i].itsB.size(); kk++)
			{
				double kappa = 8 * 3.1416 * ni * lb * qs*qs;
				VECTOR3D r_vec = dist( protein[i].itsB[ii] , protein[i].itsB[kk] );
				long double r = r_vec.GetMagnitude();
				VECTOR3D ff = r_vec ^ ( ( protein[i].itsB[ii]->q * protein[i].itsB[kk]->q * lb * exp(-kappa * r )
                          / (r*r) ) * (kappa + 1/r ) );

				protein[i].itsB[ii]->eforce += ff;
				protein[i].itsB[kk]->eforce -= ff;
			}
		}
	}
}

void update_ES_forces(vector<SUBUNIT>& protein, double lb, double ni, double qs){
    for (int i=0; i<protein.size(); i++){
		for (int ii=0; ii<protein[i].itsB.size(); ii++)
		{
			protein[i].itsB[ii]->eforce = VECTOR3D(0,0,0);
		}
    }
    
    for (int i = 0; i < protein.size() ; i++)									//intermolecular forces loop
	{
		for (int j = 0; j < protein.size(); j++)
		{
			if (protein[i].id != protein[j].id){
				for (int ii = 0; ii < protein[i].itsB.size(); ii++)
				{
					for (int jj = 0; jj < protein[j].itsB.size(); jj++)
					{
						double kappa = 8 * 3.1416 * ni * lb * qs*qs;
						VECTOR3D r_vec = dist( protein[i].itsB[ii] , protein[j].itsB[jj] );
						long double r = r_vec.GetMagnitude();
						VECTOR3D ff = r_vec ^ ( ( protein[i].itsB[ii]->q * protein[j].itsB[jj]->q * lb * exp(-kappa * r )
						/ (r * r) ) * (kappa + 1/r) );
						
						
						protein[i].itsB[ii]->eforce += ff;
						
					}
				}
			}
			
		}
		for (int ii = 0; ii < protein[i].itsB.size(); ii++)							//intramolecular forces loop
		{
			for (int kk = ii + 1; kk < protein[i].itsB.size(); kk++)
			{
				double kappa = 8 * 3.1416 * ni * lb * qs*qs;
				VECTOR3D r_vec = dist( protein[i].itsB[ii] , protein[i].itsB[kk] );
				long double r = r_vec.GetMagnitude();
				VECTOR3D ff = r_vec ^ ( ( protein[i].itsB[ii]->q * protein[i].itsB[kk]->q * lb * exp(-kappa * r )
                          / (r*r) ) * (kappa + 1/r ) );

				protein[i].itsB[ii]->eforce += ff;
				protein[i].itsB[kk]->eforce -= ff;
			}
		}
	}
	
	
 }
 
 void update_LJ_forces_pairlist(vector<SUBUNIT>& protein, double ecut, vector<PAIR>& lj_pairlist){
	 
	 for ( int i = 0; i<protein.size(); i++) {
		 for (int ii = 0; ii < protein[i].itsB.size(); ii++){
			 protein[i].itsB[ii]->ljforce = VECTOR3D(0,0,0);
		 }
	 }
	 
	 for (int i=0; i<lj_pairlist.size(); i++){
		 
		 VECTOR3D r_vec = dist( lj_pairlist[i].itsB[0] , lj_pairlist[i].itsB[1] );
		 long double r = r_vec.GetMagnitude();
		 //double r2 = (r*r);
		 double r6 ;
		 double r12;
		 double sigma6;
		 double shc = 1.2;
		 double elj = lj_pairlist[i].epsilon;
		 double sig1 = lj_pairlist[i].itsB[0]->sigma;
		 double sig2 = lj_pairlist[i].itsB[1]->sigma;
		 double del = (sig1+sig2)/2 - shc;
		 
		 if (lj_pairlist[i].type==0 && r < (del+1.12246205*shc)){
			 sigma6 = pow(shc,6);
			 double sigma12 = sigma6 * sigma6;
			 r6 = pow((r-del),6);
			 r12 = r6 * r6;
			 lj_pairlist[i].itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
			 lj_pairlist[i].itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
		 } else if (lj_pairlist[i].type==1 && r < ((del+1.12246205*shc)*ecut)){
			 r6 = pow((r-del),6);
			 r12 = r6 * r6;
			 sigma6 = pow(shc,6);
			 double sigma12 = sigma6 * sigma6;
			 lj_pairlist[i].itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
			 lj_pairlist[i].itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
		 } else {
			 lj_pairlist[i].itsB[0]->ljforce += 0;
			 lj_pairlist[i].itsB[1]->ljforce += 0;
		 }
	 }
 }

void update_LJ_forces(vector<SUBUNIT>& protein, double ecut, vector<PAIR>& lj_pairlist){

	for (int i = 0; i < protein.size(); i++){
		for (int ii = 0; ii < protein[i].itsB.size(); ii++){
			protein[i].itsB[ii]->ljforce = VECTOR3D(0,0,0);
		}
	}
	
	
		for (int i=0; i < protein.size(); i++)
		{
			for (int ii=0; ii<protein[i].itsB.size(); ii++)
				{
					for (int n = 0; n < protein[i].itsB[ii]->itsP.size(); n++)
						{
								VECTOR3D r_vec = dist( protein[i].itsB[ii] , protein[i].itsB[ii]->itsP[n]->itsB[1] );
								if ( r_vec.GetMagnitude() == 0) r_vec = dist( protein[i].itsB[ii] , protein[i].itsB[ii]->itsP[n]->itsB[0] );
								long double r = r_vec.GetMagnitude();
								//double r2 = (r*r);
								double r6 ;
								double r12;
								double sigma6;
								double shc = 1.2;
								double elj = protein[i].itsB[ii]->itsP[n]->epsilon;
								double sig1 = protein[i].itsB[ii]->itsP[n]->itsB[1]->sigma;
								double sig2 = protein[i].itsB[ii]->itsP[n]->itsB[0]->sigma;
								double del = (sig1+sig2)/2 - shc;

								if (protein[i].itsB[ii]->itsP[n]->type == 0 && r < (del+1.12246205*shc) ){
									sigma6 = pow(shc,6);
									double sigma12 = sigma6 * sigma6;
									r6 = pow((r-del),6);
									r12 = r6 * r6;
									protein[i].itsB[ii]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
								} else if (protein[i].itsB[ii]->itsP[n]->type == 1 && r < ((del+1.12246205*shc)*ecut) ){
									r6 = pow((r-del),6);
									r12 = r6 * r6;
									sigma6 = pow(shc,6);
									double sigma12 = sigma6 * sigma6;
									protein[i].itsB[ii]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
								} else {
									protein[i].itsB[ii]->ljforce += 0;
								}
							//}
						}
				}
		}
}

void update_ES_forces_simplified(vector<BEAD>& subunit_bead, double lb, double ni, double qs)
 {
	 
	 for (int i=0; i < subunit_bead.size(); i++)
	 {
		 subunit_bead[i].eforce = VECTOR3D(0,0,0);
	 }
	 
	 for (int i = 0; i < subunit_bead.size(); i++)
	 {
		 for (int j = 0; j < subunit_bead.size(); j++)
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
 
 void update_LJ_forces_simplified(vector<BEAD>& subunit_bead, double ecut, vector<vector<int> > lj_a){
	 for (int i = 0; i < subunit_bead.size(); i++){
		 subunit_bead[i].ljforce = VECTOR3D(0,0,0);
	 }
	 for (int i = 0; i < subunit_bead.size(); i++){
		 for (int j = 0; j < subunit_bead.size(); j++){
			 if (subunit_bead[i].itsS[0]->id != subunit_bead[j].itsS[0]->id){
				 VECTOR3D r_vec = dist( &subunit_bead[i] , &subunit_bead[j] );
				 long double r = r_vec.GetMagnitude();
				 //double r2 = (r*r);
				 double r6 ;
				 double r12;
				 double sigma6;
				 double shc = 1.2;
				 
				 double sig1 = subunit_bead[i].sigma;
				 double sig2 = subunit_bead[j].sigma;
				 double del = (sig1+sig2)/2 - shc;
				 bool lj_attractive = false;
				 for (int k = 0; k < lj_a[0].size(); k++){
					 if ( subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]  ){
						 lj_attractive = true;
					 }
				 }
				 if ( r < (del+1.12246205*shc) && lj_attractive == false){							//Attractive
					 sigma6 = pow(shc,6);
					 double sigma12 = sigma6 * sigma6;
					 double elj = 1;//subunit_bead[j].epsilon;
					 r6 = pow((r-del),6);
					 r12 = r6 * r6;
					 subunit_bead[i].ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
				 } else if ( r < ((del+1.12246205*shc)) && lj_attractive == true ){			//Repulsive
					 r6 = pow((r-del),6);
					 r12 = r6 * r6;
					 sigma6 = pow(shc,6);
					 double sigma12 = sigma6 * sigma6;
					 double elj = 2;//subunit_bead[j].epsilon;
					 subunit_bead[i].ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
				 } else {
					 subunit_bead[i].ljforce += 0;
				 }
			 }
		 }
	 }
 }
 
 
 
 
