#ifndef SOUFFLE_ENERGIES_H
#define SOUFFLE_ENERGIES_H

#include "vector3d.h"
#include <vector>

class BEAD;
class SUBUNIT;

void update_LJ_ES_energies_simplified(std::vector<BEAD>& subunit_bead, double ecut, std::vector<std::vector<int> > lj_a, double elj_att, double lb, double ni, double qs, double ecut_el, double kappa);

long double particle_kinetic_energy(std::vector <BEAD> &subunit_bead);

void update_LJ_ES_energies_bigbeads(std::vector<BEAD>& big_bead, std::vector<BEAD>& subunit_bead, double ecut, double lb, double ni, double qs, double ecut_el, double kappa, double shc);

#endif //SOUFFLE_ENERGIES_H
