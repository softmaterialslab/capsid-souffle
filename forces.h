#ifndef LEMONSOUFFLE_FORCES_H
#define LEMONSOUFFLE_FORCES_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>

class BEAD;
class PAIR;
class SUBUNIT;


void update_ES_forces(std::vector<SUBUNIT>& protein, double lb, double ni, double qs);

void update_LJ_forces(std::vector<SUBUNIT>& protein, double ecut, std::vector<PAIR>& lj_pairlist);









#endif //LEMONSOUFFLE_FORCES_H
