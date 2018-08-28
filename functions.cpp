//
// Created by lauren on 2/1/18.
//

#include<iostream>
#include <fstream>
#include <iomanip>
#include"functions.h"
#include "bead.h"
#include "LJpair.h"
#include "subunit.h"
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









void dress_up(vector<EDGE> &subunit_edge, vector<FACE> &subunit_face){
    for (unsigned int i=0; i<subunit_edge.size(); i++){
        subunit_edge[i].update_length();
    }
    for (unsigned int i=0; i<subunit_face.size();i++){
        subunit_face[i].update_area_normal();
    }
}










// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm)
{
    string inPath= "sub_beads.traj.out";
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

