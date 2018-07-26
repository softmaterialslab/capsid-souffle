#ifndef LEMONSOUFFLE_ENERGIES_H
#define LEMONSOUFFLE_ENERGIES_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>

class BEAD;
class PAIR;
class UNIT;


void update_LJ_energies(std::vector<BEAD>& sub_beads, double ecut, std::vector<PAIR>& sub_pairlist);

void update_ES_energies(std::vector<UNIT>& protein, double lb, double ni, double qs);

long double particle_kinetic_energy(std::vector <BEAD> &sub_beads);









#endif //LEMONSOUFFLE_ENERGIES_H
