//
// Created by lauren on 2/8/18.
//

#ifndef LEMONSOUFFLE_BOX_H
#define LEMONSOUFFLE_BOX_H

#include <vector>
#include "vector3d.h"

class BEAD;
class UNIT;


class BOX {         //This class might be unnecessary. Absorb features into another class?
public:
    //member variables
    VECTOR3D size;                      // x y z size of the box
    std::vector<UNIT *> itsU;            // Units in the box
    std::vector<BEAD *> itsB;            // Beads in the box


};


#endif //LEMONSOUFFLE_BOX_H
