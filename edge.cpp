//
// Created by lauren on 1/25/18.
//

#include <assert.h>
#include <iostream>
#include "edge.h"
#include "bead.h"
#include "face.h"
#include "functions.h"





BEAD* EDGE::opposite(BEAD *theB) {          //Finds particle opposite to the one referred in a bond
    assert(itsB.size()==2);
    if (itsB[0]==theB){
        return itsB[1];
    }else if (itsB[1]==theB){
        return itsB[0];
    } return NULL;
}

FACE* EDGE::opposite(FACE* theF){           //Find face on the opposite side of the one referred over a bond
    assert(itsF.size()==2);
    if (itsF[0]==theF) {
        return itsF[1];
    } else if (itsF[1]==theF) {
        return itsF[0];
    } return NULL;
}


void EDGE::update_length() {                //update's the length of a bond (PBC)
    VECTOR3D r_vec = (itsB[0]->pos - itsB[1]->pos);
    VECTOR3D box = itsB[0]->bx;
    if (r_vec.x>box.x/2) r_vec.x -= box.x;
    if (r_vec.x<-box.x/2) r_vec.x += box.x;
    if (r_vec.y>box.y/2) r_vec.y -= box.y;
    if (r_vec.y<-box.y/2) r_vec.y += box.y;
    if (r_vec.z>box.z/2) r_vec.z -= box.z;
    if (r_vec.z<-box.z/2) r_vec.z += box.z;
    length = std::sqrt( r_vec * r_vec );              //bondlength of a garybond
    lengthvec = r_vec;                                //length vector (useful for some computations)
}



void EDGE::update_bending_energy(double Kb)
{
    std::vector<BEAD*> r;
    r.resize(4);                        //          r[2]
    r[0] = itsB[0];                     //         /    \ <---f[0]
    r[1] = itsB[1];                     //     r[0]------r[1]
    r[2] = itsF[0]->across(this);       //         \    / <---f[1]
    r[3] = itsF[1]->across(this);       //          r[3]


    for(int i=0; i<4; i++)
    {
        r[i]->be +=  Kb * 0.25 * ( 1 - (itsF[0]->normvec * itsF[1]->normvec));
    }

}

void EDGE::update_bending_forces(double Kb)
{
    std::vector<BEAD*> r;
    r.resize(4);                        //          r[2]
    r[0] = itsB[0];                     //         /    \ <---f[0]
    r[1] = itsB[1];                     //     r[0]------r[1]
    r[2] = itsF[0]->across(this);       //         \    / <---f[1]
    r[3] = itsF[1]->across(this);       //          r[3]


    double S = (itsF[0]->direction * itsF[1]->direction);
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

    if (wrt == itsB[0]) //if this bead is itsB[0]...
    {

        VECTOR3D dist23 = r13;      //dist## (left) and r(##) (right) may not be the same numbers. That's OK. dist## refers to the original
        VECTOR3D dist32 = r13^-1;   //orientation (shown in comments in update_bending_forces), whereas r## refers to the orientation
        VECTOR3D dist20 = r01^-1;   //with respect to the bead that was input. r## will be different for each if statement.
        VECTOR3D dist21 = r12;

        result.x = dist23.y * ( (dist20.x * dist21.y) - (dist20.y * dist21.x) )
                 - dist21.y * ( (dist20.x * dist32.y) - (dist20.y * dist32.x) )
                 + dist32.z * ( (dist20.z * dist21.x) - (dist20.x * dist21.z) )
                 + dist21.z * ( (dist20.z * dist32.x) - (dist20.x * dist32.z) );


        result.y = dist32.x * ( (dist20.x * dist21.y) - (dist20.y * dist21.x) )
                 + dist21.x * ( (dist20.x * dist32.y) - (dist20.y * dist32.x) )
                 + dist23.z * ( (dist20.y * dist21.z) - (dist20.z * dist21.y) )
                 - dist21.z * ( (dist20.y * dist32.z) - (dist20.z * dist32.y) );

        result.z = dist23.x * ( (dist20.z * dist21.x) - (dist20.x * dist21.z) )
                 - dist21.x * ( (dist20.z * dist32.x) - (dist20.x * dist32.z) )
                 + dist32.y * ( (dist20.y * dist21.z) - (dist20.z * dist21.y) )
                 + dist21.y * ( (dist20.y * dist32.z) - (dist20.z * dist32.y) );

    } else if (wrt == itsB[1]){
        VECTOR3D dist23 = r03;
        VECTOR3D dist32 = r03^-1;
        VECTOR3D dist20 = r01;
        VECTOR3D dist21 = r02;

        result.x = dist23.y * ( (dist20.x * dist21.y) - (dist20.y * dist21.x) )
                   - dist21.y * ( (dist20.x * dist32.y) - (dist20.y * dist32.x) )
                   + dist32.z * ( (dist20.z * dist21.x) - (dist20.x * dist21.z) )
                   + dist21.z * ( (dist20.z * dist32.x) - (dist20.x * dist32.z) );


        result.y = dist32.x * ( (dist20.x * dist21.y) - (dist20.y * dist21.x) )
                   + dist21.x * ( (dist20.x * dist32.y) - (dist20.y * dist32.x) )
                   + dist23.z * ( (dist20.y * dist21.z) - (dist20.z * dist21.y) )
                   - dist21.z * ( (dist20.y * dist32.z) - (dist20.z * dist32.y) );

        result.z = dist23.x * ( (dist20.z * dist21.x) - (dist20.x * dist21.z) )
                   - dist21.x * ( (dist20.z * dist32.x) - (dist20.x * dist32.z) )
                   + dist32.y * ( (dist20.y * dist21.z) - (dist20.z * dist21.y) )
                   + dist21.y * ( (dist20.y * dist32.z) - (dist20.z * dist32.y) );
    }
    else
    {

        VECTOR3D dist32;
        for (unsigned int i=0; i < itsF[0]->itsB.size(); i++)
            if (itsF[0]->itsB[i] != itsB[0] && itsF[0]->itsB[i] != itsB[1]
                && itsF[0]->itsB[i] != wrt) {
                dist32 = r12 ^-1;
            }
        for (unsigned int i=0; i < itsF[1]->itsB.size(); i++) {
            if (itsF[1]->itsB[i] != itsB[0] && itsF[1]->itsB[i] != itsB[1]
                && itsF[1]->itsB[i] != wrt) {
                dist32 = r13 ^-1;
            }
        }


        VECTOR3D dist20 = r01^-1;


        result.x = dist20.y * ( (dist20.x * dist32.y) - (dist20.y * dist32.x) )
                 - dist20.z * ( (dist20.z * dist32.x) - (dist20.x * dist32.z) );

        result.y = - dist20.x * ( (dist20.x * dist32.y) - (dist20.y * dist32.x) )
                   + dist20.z * ( (dist20.y * dist32.z) - (dist20.z * dist32.y) );

        result.z = dist20.x * ( (dist20.z * dist32.x) - (dist20.x * dist32.z) )
                 - dist20.y * ( (dist20.y * dist32.z) - (dist20.z * dist32.y) );

    }
    return result;
}

VECTOR3D EDGE::get_grad0(BEAD* wrt, VECTOR3D r01, VECTOR3D r12, VECTOR3D r02){

    VECTOR3D result;

    if (wrt == itsB[0])
    {

        VECTOR3D dist21 = r12;
        VECTOR3D dist10 = r02^-1;
        VECTOR3D dist20 = r01^-1;

        result.x = - dist21.y * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) )
                   + dist21.z * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) );
        result.y = dist21.x * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) )
                 - dist21.z * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) );
        result.z = - dist21.x * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) )
                   + dist21.y * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) );

    } else if (wrt == itsB[1]){
        VECTOR3D dist21 = r02;
        VECTOR3D dist10 = r12^-1;
        VECTOR3D dist20 = r01;


        result.x = - dist21.y * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) )
                   + dist21.z * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) );
        result.y = dist21.x * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) )
                   - dist21.z * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) );
        result.z = - dist21.x * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) )
                   + dist21.y * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) );
    }
    else if (wrt == itsF[0]->across(this))
    {


        VECTOR3D dist10 = r02^-1;
        VECTOR3D dist20 = r01^-1;

        result.x =   dist20.y * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) )
                   - dist20.z * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) );
        result.y = - dist20.x * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) )
                   + dist20.z * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) );
        result.z =   dist20.x * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) )
                   - dist20.y * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) );

    } else if (wrt == itsF[1]->across(this))
    {
        result.x = 0;
        result.y = 0;
        result.z = 0;
    }
    result = result ^ ( 0.25/(itsF[0]->a) );
    return result;
}


VECTOR3D EDGE::get_grad1(BEAD* wrt, VECTOR3D r13, VECTOR3D r01, VECTOR3D r03){

    VECTOR3D result;

    if (wrt == itsB[0])
    {

        VECTOR3D dist21 = r13;
        VECTOR3D dist10 = r03^-1;
        VECTOR3D dist20 = r01^-1;

        result.x = dist21.z * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) )
                 - dist21.y * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) );
        result.y =-dist21.z * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) )
                 + dist21.x * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) );
        result.z = dist21.y * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) )
                 - dist21.x * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) );
    } else if (wrt == itsB[1]){
        VECTOR3D dist21 = r03;
        VECTOR3D dist10 = r13^-1;
        VECTOR3D dist20 = r01;

        result.x = dist21.z * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) )
                   - dist21.y * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) );
        result.y =-dist21.z * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) )
                  + dist21.x * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) );
        result.z = dist21.y * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) )
                   - dist21.x * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) );
    }
    else if (wrt == itsF[1]->across(this))
    {

        VECTOR3D dist10 = r03^-1;
        VECTOR3D dist20 = r01^-1;

        result.x =-dist20.z * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) )
                 + dist20.y * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) );
        result.y = dist20.z * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) )
                 - dist20.x * ( (dist10.x * dist20.y) - (dist10.y * dist20.x) );
        result.z =-dist20.y * ( (dist10.y * dist20.z) - (dist10.z * dist20.y) )
                 + dist20.x * ( (dist10.z * dist20.x) - (dist10.x * dist20.z) );
    } else if (wrt == itsF[0]->across(this))
    {
        result.x = 0;
        result.y = 0;
        result.z = 0;
    }
    result = result ^ ( 0.25/(itsF[1]->a) );
    return result;
}
