//
// Created by lauren on 1/25/18.
//

#include <assert.h>
#include <iostream>
#include "edge.h"
#include "bead.h"
#include "face.h"
#include "functions.h"


using namespace std;


BEAD* EDGE::opposite(BEAD *theB) {          //Finds particle opposite to the one referred in a edge
    assert(itsB.size()==2);
    if (itsB[0]==theB){
        return itsB[1];
    }else if (itsB[1]==theB){
        return itsB[0];
    } return NULL;
}

FACE* EDGE::opposite(FACE* theF){           //Find face on the opposite side of the one referred over a edge
    assert(itsF.size()==2);
    if (itsF[0]==theF) {
        return itsF[1];
    } else if (itsF[1]==theF) {
        return itsF[0];
    } return NULL;
}


void EDGE::update_length() {                //update's the length of a edge (PBC)
	
    lengthvec = itsB[0]->pos - itsB[1]->pos;
	 
    VECTOR3D box = (itsB[0]->bx);
    VECTOR3D halfbox = (itsB[0]->bx) ^ 0.5;
    if (lengthvec.x > halfbox.x) lengthvec.x -= box.x;
    else if (lengthvec.x < -halfbox.x) lengthvec.x += box.x;
    if (lengthvec.y > halfbox.y) lengthvec.y -= box.y;
    else if (lengthvec.y < -halfbox.y) lengthvec.y += box.y;
    if (lengthvec.z > halfbox.z) lengthvec.z -= box.z;
    else if (lengthvec.z < -halfbox.z) lengthvec.z += box.z;
	 
    length = lengthvec.GetMagnitude();              //edgelength
//    VECTOR3D r_vec = dist(itsB[0], itsB[1]);
//    length = r_vec.GetMagnitude();              //edgelength
//    lengthvec = r_vec;                          //length vector (useful for some computations)
}



void EDGE::update_bending_energy(double Kb)
{
    vector<BEAD*> r;
    r.resize(4);                        //          r[2]
    r[0] = itsB[0];                     //         /    \ <---f[0]
    r[1] = itsB[1];                     //     r[0]------r[1]
    r[2] = itsF[0]->across(this);       //         \    / <---f[1]
    r[3] = itsF[1]->across(this);       //          r[3]

    for(int i=0; i<4; i++) {
        r[i]->be +=  Kb * 0.25 * ( 1 - (itsF[0]->normvec * itsF[1]->normvec));
    }
}

void EDGE::update_bending_forces(double Kb)
{
    vector<BEAD*> r;
    r.resize(4);                        //          r[2]
    r[0] = itsB[0];                     //         /    \ <---f[0]
    r[1] = itsB[1];                     //     r[0]------r[1]
    r[2] = itsF[0]->across(this);       //         \    / <---f[1]
    r[3] = itsF[1]->across(this);       //          r[3]


    double S = (itsF[0]->areavector * itsF[1]->areavector) * 4;
    VECTOR3D dist01 = lengthvec;
    VECTOR3D dist02 = dist(r[0],r[2]);
    VECTOR3D dist03 = dist(r[0],r[3]);
    VECTOR3D dist12 = dist(r[1],r[2]);
    VECTOR3D dist13 = dist(r[1],r[3]);

    for(int i=0; i<4; i++){

        r[i]->bforce += (get_gradS(r[i],dist13,dist01,dist12,dist02,dist03)
                    - (((get_grad0(r[i],dist01,dist12,dist02))^((S)/itsF[0]->a))
                     + ((get_grad1(r[i],dist13,dist01,dist03))^((S)/itsF[1]->a))))
                    ^ ( (Kb) / (itsF[0]->a * itsF[1]->a * 4) );
    }
}



VECTOR3D EDGE::get_gradS(BEAD* wrt, VECTOR3D r13, VECTOR3D r01, VECTOR3D r12, VECTOR3D r02, VECTOR3D r03){
                        //distances computed in update_bending_forces input into this fxn
    VECTOR3D result;
    
    if (wrt == itsB[0]) {//if this bead is itsB[0]...
       
        result = (r13 & ((r01^-1) & r12)) - (r12 & ((r01^-1) & (r13^-1))); 
        
    } else if (wrt == itsB[1]){
       
        result = (r03 & (r01 & r02)) - (r02 & (r01 & (r03^-1))); 
        
    } else {
        VECTOR3D r1x;
        for (unsigned int i=0; i < itsF[0]->itsB.size(); i++)
            if (itsF[0]->itsB[i] != itsB[0] && itsF[0]->itsB[i] != itsB[1] && itsF[0]->itsB[i] != wrt) {
                r1x = r12 ^-1;
            }
        for (unsigned int i=0; i < itsF[1]->itsB.size(); i++) {
            if (itsF[1]->itsB[i] != itsB[0] && itsF[1]->itsB[i] != itsB[1] && itsF[1]->itsB[i] != wrt) {
                r1x = r13 ^-1;
            }
        }
        result = ((r01^-1) & ((r01^-1) & r1x));
    } // else
    return result;
}

VECTOR3D EDGE::get_grad0(BEAD* wrt, VECTOR3D r01, VECTOR3D r12, VECTOR3D r02){

    VECTOR3D result;

    if (wrt == itsB[0]){ 
       
        result = (r12 & ( r02 & r01 )) ^ -1;
        
    } else if (wrt == itsB[1]){
       
        result = (r02 & ( r12 & r01 ));  
        
    } else if (wrt == itsF[0]->across(this)){
       
        result = (r01 & ( r02 & r01 )) ^ -1;
        
    } else if (wrt == itsF[1]->across(this)){
        result = VECTOR3D(0,0,0);
    }
    result = result ^ ( 0.25/(itsF[0]->a) );  
    
    return result;
}


VECTOR3D EDGE::get_grad1(BEAD* wrt, VECTOR3D r13, VECTOR3D r01, VECTOR3D r03){

    VECTOR3D result;

    if (wrt == itsB[0]){
       
        result = (r13 & ((r03^-1) & (r01^-1))) ^ -1;
        
    } else if (wrt == itsB[1]){
       
        result = (r03 & ((r13^-1) & r01)) ^ -1;

    } else if (wrt == itsF[1]->across(this)) {

        result = ((r01^-1) & ((r03^-1) & (r01^-1)));
        
    } else if (wrt == itsF[0]->across(this)) {
        result.x = 0;
        result.y = 0;
        result.z = 0;
    }
    result = result ^ ( 0.25/(itsF[1]->a) );
    
    return result;
}
