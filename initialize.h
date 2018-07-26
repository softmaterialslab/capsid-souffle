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
#include "bead.h"           //stores particle class (sub_beads)
#include "edge.h"
#include "unit.h"
#include "face.h"
#include "box.h"
#include "LJpair.h"
#include "oligomer.h"
#include "vector3d.h"           //stores VECTOR3D class

void initialize_system(std::vector<BEAD> & sub_beads,std::vector<EDGE> & sub_edges,std::vector<UNIT> & protein, \
                        std::vector<FACE> & sub_faces, VECTOR3D bxsz, BOX & tardis, std::vector<PAIR> & sub_pairlist);

void initialize_outputfile(std::ofstream & reftraj, std::ofstream & refofile);

void allonsy (double capconc,unsigned int numden, std::string file_name);

void initialize_bead_velocities(std::vector<UNIT> &protein, std::vector<BEAD> &sub_beads, double T);

void initialize_constant_bead_velocities(std::vector<UNIT> &protein, std::vector<BEAD> &sub_beads, double T);


#endif //LEMONSOUFFLE_INITIALIZE_H
