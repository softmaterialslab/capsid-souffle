//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_INITIALIZE_H
#define LEMONSOUFFLE_INITIALIZE_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include "bead.h"           //stores particle class (subunit_bead)
#include "edge.h"
#include "subunit.h"
#include "face.h"
#include "LJpair.h"
#include "oligomer.h"
#include "vector3d.h"           //stores VECTOR3D class


void initialize_outputfile(std::ofstream & reftraj, std::ofstream & refofile);

std::vector<std::vector<int> > generate_lattice (double capsomere_concentration ,unsigned int number_capsomeres, std::string file_name, double  &bondlength,  double &SIsigma,  double &SImass, double &SItime, std::vector<BEAD> &subunit_bead, std::vector<EDGE> &subunit_edge, std::vector<SUBUNIT> &protein, std::vector<FACE> &subunit_face);

void initialize_bead_velocities(std::vector<SUBUNIT> &protein, std::vector<BEAD> &subunit_bead, double T);

void initialize_constant_bead_velocities(std::vector<SUBUNIT> &protein, std::vector<BEAD> &subunit_bead, double T);


#endif //LEMONSOUFFLE_INITIALIZE_H
