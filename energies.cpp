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


void update_ES_energies(vector<UNIT>& garfield, double lb, double ni, double qs)
{
	
	for (int i = 0; i < garfield.size()-1 ; i++)							//intermolecular energies loop
	{
		for (int j = i + 1; j < garfield.size(); j++)
		{
			for (int ii = 0; ii < garfield[i].itsB.size(); ii++)
			{
				for (int jj = 0; jj < garfield[j].itsB.size(); jj++)
				{
					double kappa = 8 * 3.1416 * ni * lb * qs*qs;
					VECTOR3D r_vec = dist( garfield[i].itsB[ii] , garfield[j].itsB[jj] );
					long double r = r_vec.GetMagnitude();
					
					garfield[i].itsB[ii]->ce += (0.5 * garfield[i].itsB[ii]->q * garfield[j].itsB[jj]->q * lb * exp(-kappa*r) ) / (r);
					garfield[j].itsB[jj]->ce += (0.5 * garfield[i].itsB[ii]->q * garfield[j].itsB[jj]->q * lb * exp(-kappa*r) ) / (r);
				}
			}
		}
	}
	for (int i = 0; i < garfield.size(); i++)								//intramolecular energies loop
	{
		for (int ii = 0; ii < garfield[i].itsB.size(); ii++)
		{
			for (int kk = ii + 1; kk < garfield[i].itsB.size(); kk++)
			{
				double kappa = 8 * 3.1416 * ni * lb * qs*qs;
				VECTOR3D r_vec = dist( garfield[i].itsB[ii] , garfield[i].itsB[kk] );
				long double r = r_vec.GetMagnitude();
				
				garfield[i].itsB[ii]->ce += (0.5 * garfield[i].itsB[ii]->q * garfield[i].itsB[kk]->q * lb * exp(-kappa*r) ) / (r);
				garfield[i].itsB[kk]->ce += (0.5 * garfield[i].itsB[ii]->q * garfield[i].itsB[kk]->q * lb * exp(-kappa*r) ) / (r);
			}
		}
	}
	
}


void update_LJ_energies(vector<BEAD>& gary, double ecut, vector<PAIR>& gpair){
	
	for (int i=0; i<gpair.size(); i++){
		
		VECTOR3D r_vec = dist( gpair[i].itsB[0] , gpair[i].itsB[1] );
		long double r = r_vec.GetMagnitude();
		double r2 = r_vec.GetMagnitudeSquared();
		double r6 ;
		double sigma6;
		double elj = gpair[i].epsilon;
		double shc = 1.2;
		double sig1 = gpair[i].itsB[0]->sigma;
		double sig2 = gpair[i].itsB[1]->sigma;
		double del = (sig1+sig2)/2 - shc;
		
		if (gpair[i].type==0 && r < (del+1.12246205*shc)){
			sigma6 = pow(shc,6);
			r6 = pow((r-del),6);
			gpair[i].itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
			gpair[i].itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
		} else if (gpair[i].type==1 && r < ((del+1.12246205*shc)*ecut)){
			double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
			double ecut12 = ecut6 * ecut6;
			sigma6 = pow(shc,6);
			r6 = pow((r-del),6);
			gpair[i].itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
			(4 * elj * ((1 / ecut12) - (1 / ecut6))));
			gpair[i].itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) -
			(4 * elj * ((1 / ecut12) - (1 / ecut6))));
		} else {
			gpair[i].itsB[0]->ne += 0;
			gpair[i].itsB[1]->ne += 0;
		}
	}
}

long double particle_kinetic_energy(vector <BEAD> &gary) {          //part of thermostat
	for (unsigned int i = 0; i < gary.size(); i++)
		gary[i].update_kinetic_energy();
	long double kinetic_energy = 0.0;
	for (unsigned int i = 0; i < gary.size(); i++)
		kinetic_energy += gary[i].ke;
	return kinetic_energy;
}
