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

void update_LJ_forces(std::vector<BEAD>& gary, double ecut, std::vector<PAIR>& gpair);

void update_LJ_energies(std::vector<BEAD>& gary, double ecut, std::vector<PAIR>& gpair);

void update_ES_forces(std::vector<BEAD>& gary, double lb, double ni, double qs);

void update_ES_energies(std::vector<BEAD>& gary, double lb, double ni, double qs);

VECTOR3D dist(BEAD* A, BEAD* B);

void update_chain_xi(unsigned int j, std::vector<THERMOSTAT>& bath, double dt, long double ke);

long double particle_kinetic_energy(std::vector <BEAD> &rgary);

void dress_up(std::vector<EDGE> &gedge, std::vector<FACE> &gface);

double compute_MD_trust_factor_R(int hiteqm);




#endif //LEMONSOUFFLE_FUNCTIONS_H
