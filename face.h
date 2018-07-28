//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_FACE_H
#define LEMONSOUFFLE_FACE_H

#include <vector>
#include "vector3d.h"



class BEAD;
class EDGE;
class SUBUNIT;

class FACE
{
public:
    //member variables
    int id;                                 //subunit id
    int type;
    VECTOR3D normvec;                       //normal vector
    VECTOR3D direction;                     //direction of vector
    double a;                          //area
    std::vector<BEAD*> itsB;                //particles making up the face
    std::vector<EDGE*> itsE;               //bending edges on the face

    //member fxns
    FACE(VECTOR3D initial=(0,0,0), int id_i=0)          //constructor
    {
        normvec=initial;
        id = id_i;
    }

    void update_normal();

    void update_area();

    void update_area_normal();

    BEAD* across(EDGE* theE);
};

#endif //LEMONSOUFFLE_FACE_H
