//
// Created by lauren on 1/25/18.
//

#include "bead.h"
#include "edge.h"
#include "functions.h"

void BEAD::update_stretching_energy(double ks){
   se = 0;
   for(unsigned int i = 0; i < (itsE.size()); i++){
      se += ( 0.25 * ks * ( itsE[i]->len0 - itsE[i]->length ) * ( itsE[i]->len0 - itsE[i]->length ) );
   }
}

// consider moving to EDGE class
void BEAD::update_stretching_force(double ks){
   sforce = VECTOR3D(0,0,0);
   for(unsigned int i = 0; i < itsE.size(); i++){
      VECTOR3D r_vec = dist(this, (itsE[i]->opposite(this)));
      sforce += ((r_vec) ^ ((itsE[i]->len0 - itsE[i]->length) * ks / itsE[i]->length));
   }
}
