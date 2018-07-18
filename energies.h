#ifndef LEMONSOUFFLE_ENERGIES_H
#define LEMONSOUFFLE_ENERGIES_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>

class BEAD;
class PAIR;
class UNIT;


void update_LJ_energies(std::vector<BEAD>& gary, double ecut, std::vector<PAIR>& gpair);

void update_ES_energies(std::vector<BEAD>& gary, double lb, double ni, double qs);

long double particle_kinetic_energy(std::vector <BEAD> &rgary);









#endif //LEMONSOUFFLE_ENERGIES_H
