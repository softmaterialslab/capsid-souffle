#ifndef LEMONSOUFFLE_FORCES_H
#define LEMONSOUFFLE_FORCES_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>

class BEAD;
class PAIR;
class UNIT;


void update_ES_forces(std::vector<UNIT>& garfield, double lb, double ni, double qs);

void update_LJ_forces(std::vector<BEAD>& gary, double ecut, std::vector<PAIR>& gpair);









#endif //LEMONSOUFFLE_FORCES_H
