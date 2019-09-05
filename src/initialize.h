//
// Created by lauren on 1/25/18.
//

#ifndef SOUFFLE_INITIALIZE_H
#define SOUFFLE_INITIALIZE_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include "bead.h"           
#include "edge.h"
#include "subunit.h"
#include "face.h"
#include "vector3d.h"          
#include "functions.h"


void initialize_outputfile(std::ofstream & reftraj, std::ofstream & refofile);

std::vector<std::vector<int> > generate_lattice (double capsomere_concentration ,unsigned int number_capsomeres, std::string file_name, double  &bondlength,  
                                                 double &SIsigma,  double &SImass, std::vector<BEAD> &subunit_bead, std::vector<EDGE> &subunit_edge, 
                                                 std::vector<SUBUNIT> &protein, std::vector<FACE> &subunit_face, bool restartFile, int &restartStep);

void initialize_bead_velocities(std::vector<SUBUNIT> &protein, std::vector<BEAD> &subunit_bead, double T);

void initialize_constant_bead_velocities(std::vector<SUBUNIT> &protein, std::vector<BEAD> &subunit_bead, double T);


#endif //SOUFFLE_INITIALIZE_H
