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
		for (int ii=0; ii<garfield[i].itsB.size(); ii++)
		{
			garfield[i].itsB[ii]->eforce = VECTOR3D(0,0,0);
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
					VECTOR3D r_vec = dist( garfield[i].itsB[ii] , garfield[j].itsB[jj] );
					long double r = r_vec.GetMagnitude();
					VECTOR3D ff = r_vec ^ ( ( garfield[i].itsB[ii]->q * garfield[j].itsB[jj]->q * lb * exp(-kappa * r )
                          / (r * r) ) * (kappa + 1/r) );


					garfield[i].itsB[ii]->eforce += ff;
					garfield[j].itsB[jj]->eforce -= ff;
				}
			}
		}
	}
	for (int i = 0; i < garfield.size(); i++)		//intramolecular forces loop
	{
		for (int ii = 0; ii < garfield[i].itsB.size(); ii++)
		{
			for (int kk = ii + 1; kk < garfield[i].itsB.size(); kk++)
			{
				double kappa = 8 * 3.1416 * ni * lb * qs*qs;
				VECTOR3D r_vec = dist( garfield[i].itsB[ii] , garfield[i].itsB[kk] );
				long double r = r_vec.GetMagnitude();
				VECTOR3D ff = r_vec ^ ( ( garfield[i].itsB[ii]->q * garfield[i].itsB[kk]->q * lb * exp(-kappa * r )
                          / (r*r) ) * (kappa + 1/r ) );

				garfield[i].itsB[ii]->eforce += ff;
				garfield[i].itsB[kk]->eforce -= ff;
			}
		}
	}
 }

void update_LJ_forces(vector<BEAD>& gary, double ecut, vector<PAIR>& gpair){

    for ( int i = 0; i<gary.size(); i++) {
       gary[i].ljforce = VECTOR3D(0,0,0);                           //Clear force
    }

    for (int i=0; i<gpair.size(); i++){

        VECTOR3D r_vec = dist( gpair[i].itsB[0] , gpair[i].itsB[1] );
		long double r = r_vec.GetMagnitude();
        //double r2 = (r*r);
        double r6 ;
        double r12;
        double sigma6;
        double shc = 1.2;
        double elj = gpair[i].epsilon;
        double sig1 = gpair[i].itsB[0]->sigma;
        double sig2 = gpair[i].itsB[1]->sigma;
        double del = (sig1+sig2)/2 - shc;

        if (gpair[i].type==0 && r < (del+1.12246205*shc)){
            sigma6 = pow(shc,6);
            double sigma12 = sigma6 * sigma6;
            r6 = pow((r-del),6);
            r12 = r6 * r6;
            gpair[i].itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
            gpair[i].itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
        } else if (gpair[i].type==1 && r < ((del+1.12246205*shc)*ecut)){
            r6 = pow((r-del),6);
            r12 = r6 * r6;
            sigma6 = pow(shc,6);
            double sigma12 = sigma6 * sigma6;
            gpair[i].itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
            gpair[i].itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
        } else {
            gpair[i].itsB[0]->ljforce += 0;
            gpair[i].itsB[1]->ljforce += 0;
        }
    }
}



