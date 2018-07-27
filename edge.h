//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_EDGE_H
#define LEMONSOUFFLE_EDGE_H

#include <vector>
#include "vector3d.h"

class BEAD;
class SUBUNIT;
class FACE;



class EDGE {
public:
    //member variables
    int id;                                         //identification
    int type;                                       //0=non-bending , 1=hard-bend , 2=soft-bend
    double length;                                  //actual length
    double len0;                                    //ideal length
    VECTOR3D lengthvec;
    std::vector<BEAD *> itsB;                        //vector of particles in the bond
    std::vector<FACE *> itsF;                       //vector of faces attached to the bond



    //member fxns
    EDGE(VECTOR3D initial = (0, 0, 0), int id_i = 0)                  //constructor
    {
        id = id_i;
        lengthvec = initial;
    }

    void update_length();

    void update_bending_forces(double Kb);

    void update_bending_energy(double Kb);

    BEAD *opposite(BEAD *theB);

    FACE *opposite(FACE *theF);

    VECTOR3D get_gradS(BEAD *wrt, VECTOR3D, VECTOR3D, VECTOR3D, VECTOR3D, VECTOR3D);

    VECTOR3D get_grad0(BEAD *wrt, VECTOR3D, VECTOR3D, VECTOR3D);

    VECTOR3D get_grad1(BEAD *wrt, VECTOR3D, VECTOR3D, VECTOR3D);
};

#endif //LEMONSOUFFLE_EDGE_H
