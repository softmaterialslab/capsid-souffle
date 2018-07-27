//
// Created by lauren on 2/22/18.
//

#ifndef LEMONSOUFFLE_OLIGOMER_H
#define LEMONSOUFFLE_OLIGOMER_H
#include <vector>
#include "vector3d.h"

class SUBUNIT;

class OLIGOMER {
public:
    //data members
    unsigned int size;
    int id;
    double mass;
    std::vector<SUBUNIT*> itsS;            //Gary's parent subunit
    VECTOR3D dummy;

    OLIGOMER(VECTOR3D initial = (0, 0, 0))       //constructor
    {
        dummy = initial;
    }


    ~OLIGOMER() {                                //destructor

    }

};
#endif //LEMONSOUFFLE_OLIGOMER_H
