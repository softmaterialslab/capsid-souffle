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


void update_ES_energies(vector<SUBUNIT>& protein, double lb, double ni, double qs)
{
	
	for (unsigned int i = 0; i < protein.size()-1 ; i++)							//intermolecular energies loop
	{
		for (unsigned int j = i + 1; j < protein.size(); j++)
		{
			for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
			{
				for (unsigned int jj = 0; jj < protein[j].itsB.size(); jj++)
				{
					double kappa = sqrt(8 * 3.1416 * ni * lb * qs*qs);
					VECTOR3D r_vec = dist( protein[i].itsB[ii] , protein[j].itsB[jj] );
					long double r = r_vec.GetMagnitude();
					
					protein[i].itsB[ii]->ce += (0.5 * protein[i].itsB[ii]->q * protein[j].itsB[jj]->q * lb * exp(-kappa*r) ) / (r);
					protein[j].itsB[jj]->ce += (0.5 * protein[i].itsB[ii]->q * protein[j].itsB[jj]->q * lb * exp(-kappa*r) ) / (r);
				}
			}
		}
	}
	for (unsigned int i = 0; i < protein.size(); i++)								//intramolecular energies loop
	{
		for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
		{
			for (unsigned int kk = ii + 1; kk < protein[i].itsB.size(); kk++)
			{
				double kappa = sqrt(8 * 3.1416 * ni * lb * qs*qs);
				VECTOR3D r_vec = dist( protein[i].itsB[ii] , protein[i].itsB[kk] );
				long double r = r_vec.GetMagnitude();
				
				protein[i].itsB[ii]->ce += (0.5 * protein[i].itsB[ii]->q * protein[i].itsB[kk]->q * lb * exp(-kappa*r) ) / (r);
				protein[i].itsB[kk]->ce += (0.5 * protein[i].itsB[ii]->q * protein[i].itsB[kk]->q * lb * exp(-kappa*r) ) / (r);
			}
		}
	}
	
}

void update_LJ_energies_pairlist(vector<SUBUNIT>& protein, double ecut, vector<PAIR>& lj_pairlist){
	
	for (unsigned int i = 0; i < protein.size(); i++){
		for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){
			protein[i].itsB[ii]->ne = 0;
		}
	}
	
	for (unsigned int i=0; i<lj_pairlist.size(); i++){
		
		VECTOR3D r_vec = dist( lj_pairlist[i].itsB[0] , lj_pairlist[i].itsB[1] );
		long double r = r_vec.GetMagnitude();
		//double r2 = r_vec.GetMagnitudeSquared();
		double r6 ;
		double sigma6;
		double elj = lj_pairlist[i].epsilon;
		double shc = 1;
		double sig1 = lj_pairlist[i].itsB[0]->sigma;
		double sig2 = lj_pairlist[i].itsB[1]->sigma;
		double del = (sig1+sig2)/2 - shc;
		
		if (lj_pairlist[i].type==0 && r < (del+1.12246205*shc)){
			sigma6 = pow(shc,6);
			r6 = pow((r-del),6);
			lj_pairlist[i].itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
			lj_pairlist[i].itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
		} else if (lj_pairlist[i].type==1 && r < ((del+1.12246205*shc)*ecut)){
			double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
			double ecut12 = ecut6 * ecut6;
			sigma6 = pow(shc,6);
			r6 = pow((r-del),6);
			lj_pairlist[i].itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
			(4 * elj * ((1 / ecut12) - (1 / ecut6))));
			lj_pairlist[i].itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
			(4 * elj * ((1 / ecut12) - (1 / ecut6))));
		} else {
			lj_pairlist[i].itsB[0]->ne += 0;
			lj_pairlist[i].itsB[1]->ne += 0;
		}
	}
}



void update_LJ_energies(vector<SUBUNIT>& protein, double ecut ){
	
	for (unsigned int i = 0; i < protein.size(); i++){
		for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){
			for (unsigned int n = 0; n < protein[i].itsB[ii]->itsP.size(); n++){
				protein[i].itsB[ii]->itsP[n]->lj_calculated = false;
			}
		}
	}

	for (unsigned int i=0; i < protein.size(); i++)
	{
		for (unsigned int ii=0; ii<protein[i].itsB.size(); ii++)
		{
			for (unsigned int n = 0; n < protein[i].itsB[ii]->itsP.size(); n++)
			{
				if (protein[i].itsB[ii]->itsP[n]->lj_calculated == false) {

					VECTOR3D r_vec = dist( protein[i].itsB[ii]->itsP[n]->itsB[0] , protein[i].itsB[ii]->itsP[n]->itsB[1] );
					long double r = r_vec.GetMagnitude();
					//double r2 = r_vec.GetMagnitudeSquared();
					double r6 ;
					double sigma6;
					double elj = protein[i].itsB[ii]->itsP[n]->epsilon;
					double shc = 1.2;
					double sig1 = protein[i].itsB[ii]->itsP[n]->itsB[0]->sigma;
					double sig2 = protein[i].itsB[ii]->itsP[n]->itsB[1]->sigma;
					double del = (sig1+sig2)/2 - shc;

					if (protein[i].itsB[ii]->itsP[n]->type == 0 && r < (del+1.12246205*shc)){
						sigma6 = pow(shc,6);
						r6 = pow((r-del),6);
						protein[i].itsB[ii]->itsP[n]->lj_calculated = true;
						protein[i].itsB[ii]->itsP[n]->itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
						protein[i].itsB[ii]->itsP[n]->itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
					} else if (protein[i].itsB[ii]->itsP[n]->type == 1 && r < ((del+1.12246205*shc)*ecut)){
						double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
						double ecut12 = ecut6 * ecut6;
						sigma6 = pow(shc,6);
						r6 = pow((r-del),6);
						protein[i].itsB[ii]->itsP[n]->lj_calculated = true;
						protein[i].itsB[ii]->itsP[n]->itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
															 (4 * elj * ((1 / ecut12) - (1 / ecut6))));
						protein[i].itsB[ii]->itsP[n]->itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
															 (4 * elj * ((1 / ecut12) - (1 / ecut6))));
					} else {
						protein[i].itsB[ii]->itsP[n]->lj_calculated = true;
						protein[i].itsB[ii]->itsP[n]->itsB[0]->ne += 0;
						protein[i].itsB[ii]->itsP[n]->itsB[1]->ne += 0;
					}
				}
			}
		}
	}
}

 

long double particle_kinetic_energy(vector <BEAD> &subunit_bead) {          //part of thermostat
	for (unsigned int i = 0; i < subunit_bead.size(); i++)
		subunit_bead[i].update_kinetic_energy();
	long double kinetic_energy = 0.0;
	for (unsigned int i = 0; i < subunit_bead.size(); i++)
		kinetic_energy += subunit_bead[i].ke;
	return kinetic_energy;
}

void update_LJ_energies_simplified(vector<BEAD>& subunit_bead, double ecut, vector<vector<int> > lj_a){
	 
	 for (unsigned int i = 0; i < subunit_bead.size(); i++){
		 for (unsigned int j = 0; j < subunit_bead.size(); j++){
			 if (subunit_bead[i].itsS[0]->id != subunit_bead[j].itsS[0]->id){
				 VECTOR3D r_vec = dist( &subunit_bead[i] , &subunit_bead[j] );
				 long double r = r_vec.GetMagnitude();
				 //double r2 = (r*r);
				 double r6 ;
				 //double r12;
				 double sigma6;
				 double shc = 1.2;
				 
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
					 sigma6 = pow(shc,6);
					 double elj = 1;//subunit_bead[j].epsilon;
					 r6 = pow((r-del),6);
					 subunit_bead[i].ne += 0.5*((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
				 } else if ( r < ((del+1.12246205*shc)*ecut) && lj_attractive == true ){			//Attractive
					 double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
					 double ecut12 = ecut6 * ecut6;
					 r6 = pow((r-del),6);
					 sigma6 = pow(shc,6);
					 double elj = 2;//subunit_bead[j].epsilon;
					 subunit_bead[i].ne += 0.5* (((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
															 (4 * elj * ((1 / ecut12) - (1 / ecut6)))));
				 } else {
					 subunit_bead[i].ne += 0;
				 }
			 }
		 }
	 }
 }
 
 void update_ES_energies_simplified(vector<BEAD>& subunit_bead, double lb, double ni, double qs)
 {
	 
	 for (unsigned int i = 0; i < subunit_bead.size(); i++)
	 {
		 for (unsigned int j = 0; j < subunit_bead.size(); j++)
		 {
			 if (subunit_bead[i].id != subunit_bead[j].id)
			 {
				 double kappa = sqrt (8 * 3.1416 * ni * lb * qs*qs);
				 VECTOR3D r_vec = dist( & subunit_bead[i] , & subunit_bead[j] );
				 long double r = r_vec.GetMagnitude();
				 subunit_bead[i].ce +=  (subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa*r) ) / (r);
				 
				 
			 }
		 }
	 }
	 
 }
