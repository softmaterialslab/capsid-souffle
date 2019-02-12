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

void forceCalculation(vector<SUBUNIT> &protein, double lb, double ni, double qs, vector<BEAD> &subunit_bead,
                      vector<PAIR> &lj_pairlist, double ecut, double ks, double bondlength, double kb,
                      vector<vector<int> > lj_a, double ecut_el, double kappa) {

    //ofstream forces("outfiles/forces.out", ios::app);

    //Common MPI Message objects
    vector<VECTOR3D> forvec(sizFVec, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> forvecGather(subunit_bead.size() + extraElements, VECTOR3D(0, 0, 0));

    //global variables
    unsigned int i, j, k;
    VECTOR3D box = subunit_bead[0].bx;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTRA MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//#pragma omp parallel for schedule(dynamic) default(shared) private(i)
    for (unsigned int i = 0; i < protein.size(); i++) {
        for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
            protein[i].itsB[ii]->update_stretching_force(ks, bondlength);
            protein[i].itsB[ii]->bforce = VECTOR3D(0, 0, 0);        //resetting bending force here
        }

        for (unsigned int mm = 0; mm < protein[i].itsE.size(); mm++) {
            if (protein[i].itsE[mm]->type != 0)                    //if it is a bending edge...
                protein[i].itsE[mm]->update_bending_forces(kb);
        }
    }
    

    //LJ forces calculation and update ES forces between subunits
#pragma omp parallel for schedule(dynamic) default(shared) private(i, j, k)
    for (i = lowerBound; i <= upperBound; i++) {
        VECTOR3D eForce = VECTOR3D(0, 0, 0);
        VECTOR3D ljForce = VECTOR3D(0, 0, 0);
        for (j = 0; j < subunit_bead.size(); j++) {

            //Add electrostatic cut offs here.
            bool electrostatic = (i != j && subunit_bead[i].q != 0 && subunit_bead[j].q != 0);
            bool lj = subunit_bead[i].itsS[0]->id != subunit_bead[j].itsS[0]->id;

            VECTOR3D r_vec = VECTOR3D(0, 0, 0);
            long double r = 0.0;

            if (electrostatic || lj) {

                r_vec = subunit_bead[i].pos - subunit_bead[j].pos;
                if (r_vec.x > box.x / 2) r_vec.x -= box.x;
                if (r_vec.x < -box.x / 2) r_vec.x += box.x;
                if (r_vec.y > box.y / 2) r_vec.y -= box.y;
                if (r_vec.y < -box.y / 2) r_vec.y += box.y;
                if (r_vec.z > box.z / 2) r_vec.z -= box.z;
                if (r_vec.z < -box.z / 2) r_vec.z += box.z;
                r = r_vec.GetMagnitude();

            }

            if (electrostatic && r < ecut_el && ecut_el != 0) {
                eForce += r_vec ^ ((subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa * r)
                                    / (r * r)) * (kappa + 1 / r));
	    } else if (electrostatic && ecut_el == 0) {
		eForce += r_vec ^ ((subunit_bead[i].q * subunit_bead[j].q * lb * exp(-kappa * r)
                                    / (r * r)) * (kappa + 1 / r));
	    }
	    
	    double sig1 = subunit_bead[i].sigma;
	    double sig2 = subunit_bead[j].sigma;
	    double shc = (sig1 +sig2)/2;
	    double del = (sig1 + sig2) / 2 - shc;
	
	    if (r >= ecut && r >= del + 1.12246205 * shc)
	      continue;
	    
            if (lj) {

                bool lj_attractive = false;
                for (k = 0; k < lj_a[0].size(); k++) {
                    if (subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]) {
                        lj_attractive = true;
                    }
                }
                
	    if (lj_attractive == false && subunit_bead[i].type != 5 && subunit_bead[j].type != 5 && subunit_bead[i].type != 6 && subunit_bead[j].type != 6 && r < (del + 1.12246205 * shc)) {
		    double r3 = r * r * r;
                    double r6 = r3 * r3;
                    double r12 = r6 * r6;
                    double elj = 1;//subunit_bead[j].epsilon;
                    ljForce += (r_vec ^ (48 * elj * ((1.0 / r12) - 0.5 * (1.0 / r6)) * (1 / (r * r))));
                }      
                else if (lj_attractive == true && r < (ecut)) {
                    double r3 = r * r * r;
                    double r6 = r3 * r3;
                    double r12 = r6 * r6;
                    double elj = 2;	//subunit_bead[j].epsilon;
                    ljForce += (r_vec ^ (48 * elj * ((1.0 / r12) - 0.5 * (1.0 / r6)) * (1 / (r * r))));
                }
                else if (lj_attractive == false && r < (del + 1.12246205 * shc)) {
                    double r3 = (r - del) * (r - del) * (r - del);
                    double r6 = r3 * r3;
                    double r12 = r6 * r6;
                    double sigma3 = shc * shc * shc;
                    double sigma6 = sigma3 * sigma3;
                    double sigma12 = sigma6 * sigma6;
                    double elj = 1;//subunit_bead[j].epsilon;
                    ljForce += (r_vec ^ (48 * elj * ((sigma12 / r12) - 0.5 * (sigma6 / r6)) *
                                         (1 / (r * (r - del)))));
                } 
            }
        }

        forvec[i - lowerBound] = eForce + ljForce;

    }

    //forvec broadcasting using all gather = gather + broadcast
    if (world.size() > 1) {
        all_gather(world, &forvec[0], forvec.size(), forvecGather);
    } else {
        for (i = lowerBound; i <= upperBound; i++)
            forvecGather[i] = forvec[i - lowerBound];
    }

    //cout << lowerBound << " " << upperBound << " "<< subunit_bead.size()<<endl;
    //Total force accumulation
    //for (i = 0; i < subunit_bead.size(); i++)
    //    subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + lj[i];

    for (unsigned int i = 0; i < subunit_bead.size(); i++)
        subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + forvecGather[i];


    forvec.clear();
    forvecGather.clear();


}



 
 
 
