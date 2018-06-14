//
// Created by lauren on 1/31/18.
//

#include "unit.h"
#include "bead.h"
#include <cstdlib>

void UNIT::update_com(){
    comvec = VECTOR3D(0,0,0);
    for(int i=0;i<itsB.size();i++){
        comvec += itsB[i]->pos;
        vsumvec = itsB[i]->vel;
    }
    com = comvec.GetMagnitude();
    vsum = vsumvec.GetMagnitude();
}

