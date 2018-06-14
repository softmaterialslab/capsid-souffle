//
// Created by lauren on 4/12/18.
//

#ifndef LEMONSOUFFLE_LJPAIR_H
#define LEMONSOUFFLE_LJPAIR_H

#include <vector>
#include "vector3d.h"

class BEAD;

class PAIR{
public:
    //members
    int type;                                        //type classification, 0 for repulsive and 1 for attractive
    double epsilon;
    double sigma;
    std::vector<BEAD *> itsB;                        //vector of particles in the pair
    VECTOR3D dummy;                                 //dummy vector

    PAIR(VECTOR3D initial = (0, 0, 0))                  //constructor
    {
        dummy = initial;                              //dummy vector
    }
};



#endif //LEMONSOUFFLE_LJPAIR_H
