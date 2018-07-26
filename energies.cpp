#include<iostream>
#include <fstream>
#include <iomanip>
#include"energies.h"
#include "bead.h"
#include "LJpair.h"
#include "unit.h"
#include "edge.h"
#include "face.h"
#include "functions.h"


using namespace std;


void update_ES_energies(vector<UNIT>& protein, double lb, double ni, double qs)
{
	
	for (int i = 0; i < protein.size()-1 ; i++)							//intermolecular energies loop
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
					
					protein[i].itsB[ii]->ce += (0.5 * protein[i].itsB[ii]->q * protein[j].itsB[jj]->q * lb * exp(-kappa*r) ) / (r);
					protein[j].itsB[jj]->ce += (0.5 * protein[i].itsB[ii]->q * protein[j].itsB[jj]->q * lb * exp(-kappa*r) ) / (r);
				}
			}
		}
	}
	for (int i = 0; i < protein.size(); i++)								//intramolecular energies loop
	{
		for (int ii = 0; ii < protein[i].itsB.size(); ii++)
		{
			for (int kk = ii + 1; kk < protein[i].itsB.size(); kk++)
			{
				double kappa = 8 * 3.1416 * ni * lb * qs*qs;
				VECTOR3D r_vec = dist( protein[i].itsB[ii] , protein[i].itsB[kk] );
				long double r = r_vec.GetMagnitude();
				
				protein[i].itsB[ii]->ce += (0.5 * protein[i].itsB[ii]->q * protein[i].itsB[kk]->q * lb * exp(-kappa*r) ) / (r);
				protein[i].itsB[kk]->ce += (0.5 * protein[i].itsB[ii]->q * protein[i].itsB[kk]->q * lb * exp(-kappa*r) ) / (r);
			}
		}
	}
	
}


void update_LJ_energies(vector<BEAD>& sub_beads, double ecut, vector<PAIR>& sub_pairlist){
	
	for (int i=0; i<sub_pairlist.size(); i++){
		
		VECTOR3D r_vec = dist( sub_pairlist[i].itsB[0] , sub_pairlist[i].itsB[1] );
		long double r = r_vec.GetMagnitude();
		double r2 = r_vec.GetMagnitudeSquared();
		double r6 ;
		double sigma6;
		double elj = sub_pairlist[i].epsilon;
		double shc = 1.2;
		double sig1 = sub_pairlist[i].itsB[0]->sigma;
		double sig2 = sub_pairlist[i].itsB[1]->sigma;
		double del = (sig1+sig2)/2 - shc;
		
		if (sub_pairlist[i].type==0 && r < (del+1.12246205*shc)){
			sigma6 = pow(shc,6);
			r6 = pow((r-del),6);
			sub_pairlist[i].itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
			sub_pairlist[i].itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
		} else if (sub_pairlist[i].type==1 && r < ((del+1.12246205*shc)*ecut)){
			double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
			double ecut12 = ecut6 * ecut6;
			sigma6 = pow(shc,6);
			r6 = pow((r-del),6);
			sub_pairlist[i].itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
			(4 * elj * ((1 / ecut12) - (1 / ecut6))));
			sub_pairlist[i].itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
			(4 * elj * ((1 / ecut12) - (1 / ecut6))));
		} else {
			sub_pairlist[i].itsB[0]->ne += 0;
			sub_pairlist[i].itsB[1]->ne += 0;
		}
	}
}

long double particle_kinetic_energy(vector <BEAD> &sub_beads) {          //part of thermostat
	for (unsigned int i = 0; i < sub_beads.size(); i++)
		sub_beads[i].update_kinetic_energy();
	long double kinetic_energy = 0.0;
	for (unsigned int i = 0; i < sub_beads.size(); i++)
		kinetic_energy += sub_beads[i].ke;
	return kinetic_energy;
}
