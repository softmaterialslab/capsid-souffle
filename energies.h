#ifndef LEMONSOUFFLE_ENERGIES_H
#define LEMONSOUFFLE_ENERGIES_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>

class BEAD;
class PAIR;
class SUBUNIT;


void update_LJ_energies(std::vector<BEAD>& subunit_bead, double ecut, std::vector<PAIR>& lj_pairlist);

void update_ES_energies(std::vector<SUBUNIT>& protein, double lb, double ni, double qs);

long double particle_kinetic_energy(std::vector <BEAD> &subunit_bead);









#endif //LEMONSOUFFLE_ENERGIES_H
