//
// Created by lauren on 1/25/18.
//

#include "face.h"
#include "bead.h"
#include "edge.h"
#include "functions.h"

//updates both area and normal
void FACE::update_area_normal(){                                            
   VECTOR3D first = dist(itsB[1],itsB[0]);
   VECTOR3D second = dist(itsB[2],itsB[1]);
   areavector = (first&second) ^ 0.5;
   a = areavector.GetMagnitude();
   normvec = (areavector) ^ (1 / a);
}

//updates the area
void FACE::update_area(){                                                   
    VECTOR3D first = dist(itsB[1],itsB[0]);
    VECTOR3D second = dist(itsB[2],itsB[1]);
    a = 0.5 * ((first&second).GetMagnitude());
}

//finds the face across from an edge
BEAD* FACE::across(EDGE* theE){
   for(unsigned int i=0;i<itsB.size();i++){
      if (this->itsB[i] != theE->itsB[0] && this->itsB[i] != theE->itsB[1]){
         return this->itsB[i];
      }
   }
   return NULL;
}
