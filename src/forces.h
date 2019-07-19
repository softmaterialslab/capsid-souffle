#ifndef LEMONSOUFFLE_FORCES_H
#define LEMONSOUFFLE_FORCES_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>

class BEAD;

class PAIR;

class SUBUNIT;


void forceCalculation(std::vector<SUBUNIT> &protein, double lb, double ni, double qs, std::vector<BEAD> &subunit_bead,
                      std::vector<PAIR> &lj_pairlist, double ecut, double ks, double bondlength, double kb,
                      std::vector<std::vector<int> > lj_a, double ecut_el, double kappa, double elj_att, int a);



#endif //LEMONSOUFFLE_FORCES_H
