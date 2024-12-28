#ifndef SOUFFLE_FORCES_H
#define SOUFFLE_FORCES_H

#include "vector3d.h"
#include <vector>

class BEAD;
class SUBUNIT;
class EDGE;
class FACE;


void forceCalculation_long(std::vector<SUBUNIT> &protein, double lb, double ni, double qs, std::vector<BEAD> &subunit_bead,
                      double ecut, double ks, double kb, std::vector<std::vector<int> > lj_a, double ecut_el, 
                      double kappa, double elj_att, bool updatePairlist, double NListCutoff);


void forceCalculation_short(std::vector<SUBUNIT> &protein, std::vector<EDGE> &subunit_edge, std::vector<FACE> &subunit_face, double ks, double kb);

#endif //SOUFFLE_FORCES_H
