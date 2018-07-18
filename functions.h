//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_FUNCTIONS_H
#define LEMONSOUFFLE_FUNCTIONS_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>


class BEAD;
class PAIR;
class UNIT;
class EDGE;
class FACE;
class THERMOSTAT;

// display progress
void ProgressBar (double fraction_completed);

VECTOR3D dist(BEAD* A, BEAD* B);

void update_chain_xi(unsigned int j, std::vector<THERMOSTAT>& bath, double dt, long double ke);

void dress_up(std::vector<EDGE> &gedge, std::vector<FACE> &gface);

double compute_MD_trust_factor_R(int hiteqm);




#endif //LEMONSOUFFLE_FUNCTIONS_H
