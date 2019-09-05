//
// Created by lauren on 1/25/18.
//

#ifndef SOUFFLE_FACE_H
#define SOUFFLE_FACE_H

#include <vector>
#include "vector3d.h"

class BEAD;
class EDGE;
class SUBUNIT;

class FACE {
public:
//member variables
   int id;                                //id of the face
   int type;
   VECTOR3D normvec;                      //normal vector
   VECTOR3D areavector;              	//area vector
   double a;                          	//area
   std::vector<BEAD*> itsB;               //beads making up the face
   std::vector<EDGE*> itsE;               //bending edges on the face

//member fxns
//constructor
   FACE(VECTOR3D initial=VECTOR3D(0,0,0), int id_i=0){ 
      normvec=initial;
      id = id_i;
   }

   void update_area();

   void update_area_normal();

   BEAD* across(EDGE* theE);
};

#endif //SOUFFLE_FACE_H
