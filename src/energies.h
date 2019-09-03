#ifndef LEMONSOUFFLE_ENERGIES_H
#define LEMONSOUFFLE_ENERGIES_H

#include "vector3d.h"
//#include "thermostat.h"
#include <vector>
//#include <cmath>

class BEAD;
class PAIR;
class SUBUNIT;


void update_LJ_energies_simplified(std::vector<BEAD>& subunit_bead, double ecut, std::vector<std::vector<int> > lj_a, double elj_att);

void update_ES_energies_simplified(std::vector<BEAD>& subunit_bead, double lb, double ni, double qs, double ecut_el, double kappa);

void update_LJ_ES_energies_simplified(std::vector<BEAD>& subunit_bead, double ecut, std::vector<std::vector<int> > lj_a, double elj_att, double lb, double ni, double qs, double ecut_el, double kappa);

long double particle_kinetic_energy(std::vector <BEAD> &subunit_bead);









#endif //LEMONSOUFFLE_ENERGIES_H
