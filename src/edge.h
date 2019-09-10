//
// Created by lauren on 1/25/18.
//

#ifndef SOUFFLE_EDGE_H
#define SOUFFLE_EDGE_H

#include <vector>
#include "vector3d.h"

class BEAD;
class SUBUNIT;
class FACE;

class EDGE {
public:
//member variables
   int id;                                      //identification
   int type;                                    //bond type; 0=non-bending , 1=hard-bend , 2=soft-bend
   double length;                               //actual length
   double len0;                                 //ideal length
   VECTOR3D lengthvec;                          //actual length vector
   std::vector<BEAD *> itsB;                    //vector of beads in the bond
   std::vector<FACE *> itsF;                    //vector of faces attached to the bond

//member fxns
//constructor
   EDGE(VECTOR3D initial = VECTOR3D(0, 0, 0), int id_i = 0){                
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

#endif //SOUFFLE_EDGE_H
