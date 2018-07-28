//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_UNIT_H
#define LEMONSOUFFLE_UNIT_H

#include <vector>
#include "vector3d.h"
#include "oligomer.h"


class BEAD;
class EDGE;
class FACE;
class BOX;
//class OLIGOMER;

class SUBUNIT
{
public:
    //member variables
    int id;                                 //subunit id
    std::vector<BEAD*> itsB;                //particles making up the subunit
    BOX* itsT;                              //it's box (tardis)
    std::vector<EDGE*> itsE;				//edges within the subunit
    std::vector<OLIGOMER> itsO;                         //its oligmer
    VECTOR3D comvec;                        //center of mass of the subunit
    VECTOR3D vsumvec;
    double com;                             // center of mass magnitude
    double vsum;                            // magnitude of average velocity

    //member fxns
    SUBUNIT(VECTOR3D initial=(0,0,0))       //constructor
    {
        vsumvec=initial;
    }

    void update_com();




};




#endif //LEMONSOUFFLE_UNIT_H
