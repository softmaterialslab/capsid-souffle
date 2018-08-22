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

std::vector<std::vector<int> > initialize_system(std::vector<BEAD> & subunit_bead,std::vector<EDGE> & subunit_edge,std::vector<SUBUNIT> & protein, \
                        std::vector<FACE> & subunit_face, VECTOR3D bxsz, std::vector<PAIR> & lj_pairlist);

void initialize_outputfile(std::ofstream & reftraj, std::ofstream & refofile);

void generate_lattice (double capsomere_concentration ,unsigned int number_capsomeres, std::string file_name, double  &bondlength,  double &SIsigma,  double &SImass, double &SItime);

void initialize_bead_velocities(std::vector<SUBUNIT> &protein, std::vector<BEAD> &subunit_bead, double T);

void initialize_constant_bead_velocities(std::vector<SUBUNIT> &protein, std::vector<BEAD> &subunit_bead, double T);


#endif //LEMONSOUFFLE_INITIALIZE_H
