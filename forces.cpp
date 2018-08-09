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



void update_ES_forces(vector<SUBUNIT>& protein, double lb, double ni, double qs){
    for (int i=0; i<protein.size(); i++){
		for (int ii=0; ii<protein[i].itsB.size(); ii++)
		{
			protein[i].itsB[ii]->eforce = VECTOR3D(0,0,0);
		}
    }
    
    for (int i = 0; i < protein.size()-1 ; i++)	//intermolecular forces loop
	{
		for (int j = i + 1; j < protein.size(); j++)
		{
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
					protein[j].itsB[jj]->eforce -= ff;
				}
			}
		}
	}
	for (int i = 0; i < protein.size(); i++)		//intramolecular forces loop
	{
		for (int ii = 0; ii < protein[i].itsB.size(); ii++)
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
		 double shc = 1;
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
			for (int n = 0; n < protein[i].itsB[ii]->itsP.size(); n++){
				protein[i].itsB[ii]->itsP[n]->lj_calculated = false;
			}
		}
	}
	
		for (int i=0; i < protein.size(); i++)
		{
			for (int ii=0; ii<protein[i].itsB.size(); ii++)
				{
					for (int n = 0; n < protein[i].itsB[ii]->itsP.size(); n++)
						{
							if (protein[i].itsB[ii]->itsP[n]->lj_calculated == false){

								VECTOR3D r_vec = dist( protein[i].itsB[ii]->itsP[n]->itsB[0] , protein[i].itsB[ii]->itsP[n]->itsB[1] );
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
                                    protein[i].itsB[ii]->itsP[n]->lj_calculated = true;
									protein[i].itsB[ii]->itsP[n]->itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
									protein[i].itsB[ii]->itsP[n]->itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
								} else if (protein[i].itsB[ii]->itsP[n]->type == 1 && r < ((del+1.12246205*shc)*ecut) ){
									r6 = pow((r-del),6);
									r12 = r6 * r6;
									sigma6 = pow(shc,6);
									double sigma12 = sigma6 * sigma6;
                                    protein[i].itsB[ii]->itsP[n]->lj_calculated = true;
									protein[i].itsB[ii]->itsP[n]->itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
									protein[i].itsB[ii]->itsP[n]->itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
								} else {
                                    protein[i].itsB[ii]->itsP[n]->lj_calculated = true;
									protein[i].itsB[ii]->itsP[n]->itsB[0]->ljforce += 0;
									protein[i].itsB[ii]->itsP[n]->itsB[1]->ljforce += 0;
								}
							}
						}
				}
		}
}



