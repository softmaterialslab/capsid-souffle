//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_FUNCTIONS_H
#define LEMONSOUFFLE_FUNCTIONS_H

#include "vector3d.h"
#include "bead.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>


//class BEAD;
class PAIR;
class SUBUNIT;
class EDGE;
class FACE;
class THERMOSTAT;

// display progress
void ProgressBar (double fraction_completed);

//VECTOR3D dist(BEAD* A, BEAD* B);

void update_chain_xi(unsigned int j, std::vector<THERMOSTAT>& bath, double dt, long double ke);

void dress_up(std::vector<EDGE> &subunit_edge, std::vector<FACE> &subunit_face);

double compute_MD_trust_factor_R(int hiteqm);


inline VECTOR3D dist(BEAD *A, BEAD *B) {   //finds distance considering periodic boundaries.
    VECTOR3D r_vec; //= (A->pos - B->pos);
    r_vec.x = A->pos.x - B->pos.x;
    r_vec.y = A->pos.y - B->pos.y;
    r_vec.z = A->pos.z - B->pos.z;
    VECTOR3D box = A->bx;
    if (r_vec.x > box.x / 2) r_vec.x -= box.x;
    if (r_vec.x < -box.x / 2) r_vec.x += box.x;
    if (r_vec.y > box.y / 2) r_vec.y -= box.y;
    if (r_vec.y < -box.y / 2) r_vec.y += box.y;
    if (r_vec.z > box.z / 2) r_vec.z -= box.z;
    if (r_vec.z < -box.z / 2) r_vec.z += box.z;
    return r_vec;
}


#endif //LEMONSOUFFLE_FUNCTIONS_H
