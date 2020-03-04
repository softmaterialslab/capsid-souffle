//
// Created by lauren on 1/25/18.
//

#ifndef SOUFFLE_UNIT_H
#define SOUFFLE_UNIT_H

#include <vector>
#include "vector3d.h"

class BEAD;
class EDGE;
class FACE;
class BOX;
class OLIGOMER;

class SUBUNIT {
public:
//member variables
   int id;                                      //subunit id
   std::vector<BEAD*> itsB;                     //particles making up the subunit
   std::vector<EDGE*> itsE;			      //edges within the subunit
   std::vector<FACE*> itsF;                     //faces within the subunit
   std::vector<OLIGOMER> itsO;                  //its oligmer
   VECTOR3D comvec;                             //center of mass of the subunit
   VECTOR3D vsumvec;
   double com;                                  // center of mass magnitude
   double vsum;                                 // magnitude of average velocity
   BEAD* centerBead;                            // id of center bead

   //member fxns
   SUBUNIT(VECTOR3D initial=VECTOR3D(0,0,0)){  //constructor
      vsumvec=initial;
   }
   void update_com();
};




#endif //SOUFFLE_UNIT_H
