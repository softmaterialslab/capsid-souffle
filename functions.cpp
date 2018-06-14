//
// Created by lauren on 2/1/18.
//

#include<iostream>
#include <fstream>
#include <iomanip>
#include"functions.h"
#include "bead.h"
#include "LJpair.h"
#include "unit.h"
#include "edge.h"
#include "face.h"


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

using namespace std;

void ProgressBar(double fraction_completed)
{
    int val = (int) (fraction_completed * 100);
    int lpad = (int) (fraction_completed * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% |%.*s%*s|", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
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



void update_ES_forces(vector<BEAD>& gary, double lb, double ni, double qs){
    for (int i=0; i<gary.size(); i++){
        gary[i].eforce = VECTOR3D(0,0,0);       //clearing forces
    }
    for (int i=0; i<gary.size()-1; i++)
    {
        for (int j=i+1; j<gary.size(); j++)
        {
            double kappa = 8 * 3.1416 * ni * lb * qs*qs;
            VECTOR3D rij = dist( &gary[i] , &gary[j] );
            VECTOR3D ff = rij ^ ( ( gary[i].q * gary[j].q * lb * exp(-kappa * rij.GetMagnitude() )
                          / rij.GetMagnitudeSquared() ) * (kappa + 1/rij.GetMagnitude() ) );

            gary[i].eforce += ff;
            gary[j].eforce -= ff;
        }
    }
}

void update_ES_energies(vector<BEAD>& gary, double lb, double ni, double qs){
    for (int i=0; i<gary.size()-1; i++)
    {
        for (int j=i+1; j<gary.size(); j++)
        {
            double kappa = 8 * 3.1416 * ni * lb * qs*qs;
            VECTOR3D rij = dist( &gary[i] , &gary[j] );

            gary[i].ce += (0.5 * gary[i].q * gary[j].q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());
            gary[j].ce += (0.5 * gary[i].q * gary[j].q * lb * exp(-kappa*rij.GetMagnitude()) ) / (rij.GetMagnitude());

        }
    }
}

VECTOR3D dist(BEAD* A, BEAD* B){                    //finds distance considering periodic boundaries.
    VECTOR3D r_vec; //= (A->pos - B->pos);
    r_vec.x = A->pos.x - B->pos.x;
    r_vec.y = A->pos.y - B->pos.y;
    r_vec.z = A->pos.z - B->pos.z;
    VECTOR3D box = A->bx;
    if (r_vec.x>box.x/2) r_vec.x -= box.x;
    if (r_vec.x<-box.x/2) r_vec.x += box.x;
    if (r_vec.y>box.y/2) r_vec.y -= box.y;
    if (r_vec.y<-box.y/2) r_vec.y += box.y;
    if (r_vec.z>box.z/2) r_vec.z -= box.z;
    if (r_vec.z<-box.z/2) r_vec.z += box.z;
    return r_vec;
}

void update_chain_xi(unsigned int j, vector<THERMOSTAT>& bath, double dt, long double ke)   //part of thermostat
{
    if (bath[j].Q == 0)
        return;
    if (j != 0) {
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) + 0.5 * dt * (1.0 / bath[j].Q) *
                                                                    (bath[j - 1].Q * bath[j - 1].xi * bath[j - 1].xi -
                                                                     bath[j].dof * bath[j].kB * bath[j].T) *
                                                                    exp(-0.25 * dt * bath[j + 1].xi);
    }else {
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) +
                     0.5 * dt * (1.0 / bath[j].Q) * (2 * ke - bath[j].dof * bath[j].kB * bath[j].T) *
                     exp(-0.25 * dt * bath[j + 1].xi);
    }
    return;
}

long double particle_kinetic_energy(vector <BEAD> &gary) {          //part of thermostat
    for (unsigned int i = 0; i < gary.size(); i++)
        gary[i].update_kinetic_energy();
    long double kinetic_energy = 0.0;
    for (unsigned int i = 0; i < gary.size(); i++)
        kinetic_energy += gary[i].ke;
    return kinetic_energy;
}


void dress_up(vector<EDGE> &gedge, vector<FACE> &gface){
    for (int i=0; i<gedge.size(); i++){
        gedge[i].update_length();
    }
    for (int i=0; i<gface.size();i++){
        gface[i].update_area_normal();
    }
}

// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm)
{
    string inPath= "gary.traj.out";
    ifstream in(inPath.c_str(), ios::in);
    if (!in)
    {
        if (!in)
            cout << "File could not be opened" << endl;
        return 0;
    }
    string dummy;
    double col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
    vector<double> ext, ke, pe;
    in >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >>col10 >> col11)
    {
        ext.push_back(col7);
        ke.push_back(col2);
        pe.push_back(col8);
    }

     double ext_mean = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_mean += ext[i];
    ext_mean = ext_mean / ext.size();
    double ke_mean = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_mean += ke[i];
    ke_mean = ke_mean / ke.size();

    double ext_sd = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
    ext_sd = ext_sd / ext.size();
    ext_sd = sqrt(ext_sd);

    double ke_sd = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
    ke_sd = ke_sd / ke.size();
    ke_sd = sqrt(ke_sd);

    double R = ext_sd / ke_sd;
//    if (world.rank() == 0)
//    {
        string outPath="R.dat";
        ofstream out (outPath.c_str());
        out << "Sample size " << ext.size() << endl;
        out << "Sd: ext, kinetic energy and R" << endl;
        out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
//    }
    cout << endl << endl << "R is: " << R;
    return R;
}

