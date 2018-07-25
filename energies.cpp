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


void update_ES_energies(vector<UNIT>& garfield, double lb, double ni, double qs){
//     for (int i=0; i<gary.size()-1; i++)
//     {
//         for (int j=i+1; j<gary.size(); j++)
//         {
//             double kappa = 8 * 3.1416 * ni * lb * qs*qs;
//             VECTOR3D rij = dist( &gary[i] , &gary[j] );
// 
//             gary[i].ce += (0.5 * gary[i].q * gary[j].q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());
//             gary[j].ce += (0.5 * gary[i].q * gary[j].q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());
// 
//         }
//     }
    
    for (int i = 0; i < garfield.size()-1 ; i++)	//intermolecular energies loop
	{
		for (int j = i + 1; j < garfield.size(); j++)
		{
			for (int ii = 0; ii < garfield[i].itsB.size(); ii++)
			{
				for (int jj = 0; jj < garfield[j].itsB.size(); jj++)
				{
					double kappa = 8 * 3.1416 * ni * lb * qs*qs;
					VECTOR3D rij = dist( garfield[i].itsB[ii] , garfield[j].itsB[jj] );

					garfield[i].itsB[ii]->ce += (0.5 * garfield[i].itsB[ii]->q * garfield[j].itsB[jj]->q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());
					garfield[j].itsB[jj]->ce += (0.5 * garfield[i].itsB[ii]->q * garfield[j].itsB[jj]->q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());
				}
			}
		}
	}
	for (int i = 0; i < garfield.size(); i++)		//intramolecular energies loop
	{
		for (int j = 0; j < garfield[i].itsB.size(); j++)
		{
			for (int k = j + 1; k < garfield[i].itsB.size(); k++)
			{
				double kappa = 8 * 3.1416 * ni * lb * qs*qs;
				VECTOR3D rij = dist( garfield[i].itsB[j] , garfield[i].itsB[k] );

				garfield[i].itsB[j]->ce += (0.5 * garfield[i].itsB[j]->q * garfield[i].itsB[k]->q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());
				garfield[i].itsB[k]->ce += (0.5 * garfield[i].itsB[j]->q * garfield[i].itsB[k]->q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());
			}
		}
	}
    
}


void update_LJ_energies(std::vector<BEAD>& gary, double ecut, std::vector<PAIR>& gpair){

    for (int i=0; i<gpair.size(); i++){

        VECTOR3D rij = dist( gpair[i].itsB[0] , gpair[i].itsB[1] );
        double r2 = rij.GetMagnitudeSquared();
        double r6 ;
        double sigma6;
        double elj = gpair[i].epsilon;
        double shc = 1.2;
        double sig1 = gpair[i].itsB[0]->sigma;
        double sig2 = gpair[i].itsB[1]->sigma;
        double del = (sig1+sig2)/2 - shc;

        if (gpair[i].type==0 && rij.GetMagnitude() < (del+1.12246205*shc)){
            sigma6 = pow(shc,6);
            r6 = pow((rij.GetMagnitude()-del),6);
            gpair[i].itsB[0]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
            gpair[i].itsB[1]->ne += 0.5 * ((4 * elj * (sigma6 / r6) * ((sigma6 / r6) - 1)) + elj);
        } else if (gpair[i].type==1 && rij.GetMagnitude() < ((del+1.12246205*shc)*ecut)){
            double ecut6 = ecut * ecut * ecut * ecut * ecut * ecut;
            double ecut12 = ecut6 * ecut6;
            sigma6 = pow(shc,6);
            r6 = pow((rij.GetMagnitude()-del),6);
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
