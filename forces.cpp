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



void update_ES_forces(vector<UNIT>& protein, double lb, double ni, double qs){
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

void update_LJ_forces(vector<BEAD>& sub_beads, double ecut, vector<PAIR>& sub_pairlist){

    for ( int i = 0; i<sub_beads.size(); i++) {
       sub_beads[i].ljforce = VECTOR3D(0,0,0);                           //Clear force
    }

    for (int i=0; i<sub_pairlist.size(); i++){

        VECTOR3D r_vec = dist( sub_pairlist[i].itsB[0] , sub_pairlist[i].itsB[1] );
		long double r = r_vec.GetMagnitude();
        //double r2 = (r*r);
        double r6 ;
        double r12;
        double sigma6;
        double shc = 1.2;
        double elj = sub_pairlist[i].epsilon;
        double sig1 = sub_pairlist[i].itsB[0]->sigma;
        double sig2 = sub_pairlist[i].itsB[1]->sigma;
        double del = (sig1+sig2)/2 - shc;

        if (sub_pairlist[i].type==0 && r < (del+1.12246205*shc)){
            sigma6 = pow(shc,6);
            double sigma12 = sigma6 * sigma6;
            r6 = pow((r-del),6);
            r12 = r6 * r6;
            sub_pairlist[i].itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
            sub_pairlist[i].itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del))) ));
        } else if (sub_pairlist[i].type==1 && r < ((del+1.12246205*shc)*ecut)){
            r6 = pow((r-del),6);
            r12 = r6 * r6;
            sigma6 = pow(shc,6);
            double sigma12 = sigma6 * sigma6;
            sub_pairlist[i].itsB[0]->ljforce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
            sub_pairlist[i].itsB[1]->ljforce -= (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) * (1 / (r*(r-del)))));
        } else {
            sub_pairlist[i].itsB[0]->ljforce += 0;
            sub_pairlist[i].itsB[1]->ljforce += 0;
        }
    }
}



