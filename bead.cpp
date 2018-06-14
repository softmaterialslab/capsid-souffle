//
// Created by lauren on 1/25/18.
//

#include "bead.h"
#include "edge.h"

void BEAD::update_stretching_energy(double ks, double vdwr)
{
    se = 0;
    for(int i=0;i< (itsE.size() );i++)      //Using 1.5 for ~5nm edge length in protein
    {
        se += ( 0.125 * ks * ( itsE[i]->len0-itsE[i]->length ) * ( itsE[i]->len0-itsE[i]->length ) );
    }
}

void BEAD::update_stretching_force(double ks, double vdwr)
{
    sforce = VECTOR3D(0,0,0);
    for(int i=0; i< itsE.size(); i++)
    {
        VECTOR3D r_vec = (pos - (itsE[i]->opposite(this)->pos));
        //VECTOR3D bx = itsT->size;
        if (r_vec.x > bx.x / 2) r_vec.x -= bx.x;
        if (r_vec.x < -bx.x / 2) r_vec.x += bx.x;
        if (r_vec.y > bx.y / 2) r_vec.y -= bx.y;
        if (r_vec.y < -bx.y / 2) r_vec.y += bx.y;
        if (r_vec.z > bx.z / 2) r_vec.z -= bx.z;
        if (r_vec.z < -bx.z / 2) r_vec.z += bx.z;

        sforce += ((r_vec) ^ ((itsE[i]->len0 - itsE[i]->length) * 0.5*ks / itsE[i]->length));
    }
}