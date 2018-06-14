//
// Created by lauren on 1/25/18.
//

#include "face.h"
#include "bead.h"
#include "edge.h"
#include "functions.h"



void FACE::update_normal()                                                  //updates the normal
{
    VECTOR3D first = dist(itsB[1],itsB[0]);
    VECTOR3D second = dist(itsB[2],itsB[1]);
    direction = first&second;                                                           //cross product
    normvec = ( first&second ) ^ ( 1/((first&second).GetMagnitude()) );
}

void FACE::update_area_normal(){                                            //updates both area and normal
    VECTOR3D first = dist(itsB[1],itsB[0]);
    VECTOR3D second = dist(itsB[2],itsB[1]);
    a = 0.5 * ((first&second).GetMagnitude());
    direction = first&second;                                                           //cross product
    normvec = ( first&second ) ^ ( 1/((first&second).GetMagnitude()) );
}

void FACE::update_area(){                                                   //updates the area
    VECTOR3D first = dist(itsB[1],itsB[0]);
    VECTOR3D second = dist(itsB[2],itsB[1]);
    a = 0.5 * ((first&second).GetMagnitude());
}

BEAD* FACE::across(EDGE* theE)                                              //finds the face across from an edge
{
    for(int i=0;i<itsB.size();i++)
    {
        if (this->itsB[i] != theE->itsB[0] && this->itsB[i] != theE->itsB[1])
        {
            return this->itsB[i];
        }
    }
    return NULL;
}
