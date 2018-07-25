#include<iostream>
#include <fstream>
#include <iomanip>
#include"forces.h"
#include "bead.h"
#include "LJpair.h"
#include "unit.h"
#include "edge.h"
#include "face.h"
#include "functions.h"


using namespace std;



void update_ES_forces(vector<UNIT>& garfield, double lb, double ni, double qs){
    for (int i=0; i<garfield.size(); i++){
		for (int j=0; j<garfield[i].itsB.size(); j++)
		{
			garfield[i].itsB[j]->eforce = VECTOR3D(0,0,0);
		}
    }
    
    for (int i = 0; i < garfield.size()-1 ; i++)	//intermolecular forces loop
	{
		for (int j = i + 1; j < garfield.size(); j++)
		{
			for (int ii = 0; ii < garfield[i].itsB.size(); ii++)
			{
				for (int jj = 0; jj < garfield[j].itsB.size(); jj++)
				{
					double kappa = 8 * 3.1416 * ni * lb * qs*qs;
					VECTOR3D rij = dist( garfield[i].itsB[ii] , garfield[j].itsB[jj] );
					VECTOR3D ff = rij ^ ( ( garfield[i].itsB[ii]->q * garfield[j].itsB[jj]->q * lb * exp(-kappa * rij.GetMagnitude() )
                          / rij.GetMagnitudeSquared() ) * (kappa + 1/rij.GetMagnitude() ) );

					garfield[i].itsB[ii]->eforce += ff;
					garfield[j].itsB[jj]->eforce -= ff;
				}
			}
		}
	}
	for (int i = 0; i < garfield.size(); i++)		//intramolecular forces loop
	{
		for (int j = 0; j < garfield[i].itsB.size(); j++)
		{
			for (int k = j + 1; k < garfield[i].itsB.size(); k++)
			{
				double kappa = 8 * 3.1416 * ni * lb * qs*qs;
				VECTOR3D rij = dist( garfield[i].itsB[j] , garfield[i].itsB[k] );
				VECTOR3D ff = rij ^ ( ( garfield[i].itsB[j]->q * garfield[i].itsB[k]->q * lb * exp(-kappa * rij.GetMagnitude() )
                          / rij.GetMagnitudeSquared() ) * (kappa + 1/rij.GetMagnitude() ) );

				garfield[i].itsB[j]->eforce += ff;
				garfield[i].itsB[k]->eforce -= ff;
			}
		}
	}
//     for (int i=0; i<gary.size()-1; i++)
//     {
//         for (int j=i+1; j<gary.size(); j++)
//         {
//             double kappa = 8 * 3.1416 * ni * lb * qs*qs;
//             VECTOR3D rij = dist( &gary[i] , &gary[j] );
//             VECTOR3D ff = rij ^ ( ( gary[i].q * gary[j].q * lb * exp(-kappa * rij.GetMagnitude() )
//                           / rij.GetMagnitudeSquared() ) * (kappa + 1/rij.GetMagnitude() ) );
// 
//             gary[i].eforce += ff;
//             gary[j].eforce -= ff;
//         }
//     }
 }

void update_LJ_forces(std::vector<BEAD>& gary, double ecut, std::vector<PAIR>& gpair){

    for ( int i = 0; i<gary.size(); i++) {
       gary[i].ljforce = VECTOR3D(0,0,0);                           //Clear force
    }

    for (int i=0; i<gpair.size(); i++){

        VECTOR3D rij = dist( gpair[i].itsB[0] , gpair[i].itsB[1] );
        //double r2 = rij.GetMagnitudeSquared();
        double r6 ;
        double r12;
        double sigma6;
        double shc = 1.2;
        double elj = gpair[i].epsilon;
        double sig1 = gpair[i].itsB[0]->sigma;
        double sig2 = gpair[i].itsB[1]->sigma;
        double del = (sig1+sig2)/2 - shc;

        if (gpair[i].type==0 && rij.GetMagnitude() < (del+1.12246205*shc)){
            sigma6 = pow(shc,6);
            double sigma12 = sigma6 * sigma6;
            r6 = pow((rij.GetMagnitude()-del),6);
            r12 = r6 * r6;
            gpair[i].itsB[0]->ljforce += (rij ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (rij.GetMagnitude()*(rij.GetMagnitude()-del))) ));
            gpair[i].itsB[1]->ljforce -= (rij ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (rij.GetMagnitude()*(rij.GetMagnitude()-del))) ));
        } else if (gpair[i].type==1 && rij.GetMagnitude() < ((del+1.12246205*shc)*ecut)){
            r6 = pow((rij.GetMagnitude()-del),6);
            r12 = r6 * r6;
            sigma6 = pow(shc,6);
            double sigma12 = sigma6 * sigma6;
            gpair[i].itsB[0]->ljforce += (rij ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (rij.GetMagnitude()*(rij.GetMagnitude()-del)))));
            gpair[i].itsB[1]->ljforce -= (rij ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (rij.GetMagnitude()*(rij.GetMagnitude()-del)))));
        } else {
            gpair[i].itsB[0]->ljforce += 0;
            gpair[i].itsB[1]->ljforce += 0;
        }
    }
}



