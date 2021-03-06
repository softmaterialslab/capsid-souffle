//
// Created by lauren on 1/25/18.
//

#ifndef SOUFFLE_FUNCTIONS_H
#define SOUFFLE_FUNCTIONS_H

#include <iostream>
#include <string>
#include <sstream>
#include "vector3d.h"
#include "bead.h"
#include <vector>
#include <boost/filesystem.hpp>
#include <algorithm>
#include "RunningStat.h"

class SUBUNIT;
class EDGE;
class FACE;
class THERMOSTAT;
class OLIGOMER;

// display progress
void ProgressBar (double fraction_completed);

void update_chain_xi(unsigned int j, std::vector<THERMOSTAT>& bath, double dt, long double ke);

void dress_up(std::vector<EDGE> &subunit_edge, std::vector<FACE> &subunit_face);

void update_pairlist(unsigned int i, std::vector<SUBUNIT> &protein, double NListCutoff, VECTOR3D box, VECTOR3D hbox);

double compute_MD_trust_factor_R(int &hiteqm, bool &done, std::string directory) ;

std::vector<std::string> getFileNames(std::string path = ".") ;

void filter(std::vector<std::string>& strings, std::string pattern) ;

bool numeric_string_compare(const std::string& s1, const std::string& s2) ;

//finds distance considering periodic boundaries.
inline VECTOR3D dist(BEAD *A, BEAD *B) {   
   VECTOR3D r_vec; //= (A->pos - B->pos);
   r_vec.x = A->pos.x - B->pos.x;
   r_vec.y = A->pos.y - B->pos.y;
   r_vec.z = A->pos.z - B->pos.z;
   VECTOR3D box = A->bx;
   VECTOR3D hbox = A->hbx;
   if (r_vec.x > hbox.x) r_vec.x -= box.x;
   else if (r_vec.x < -hbox.x) r_vec.x += box.x;
   if (r_vec.y > hbox.y) r_vec.y -= box.y;
   else if (r_vec.y < -hbox.y) r_vec.y += box.y;
   if (r_vec.z > hbox.z) r_vec.z -= box.z;
   else if (r_vec.z < -hbox.z) r_vec.z += box.z;
   return r_vec;
}

bool isRestart (std::string& s);

#endif //SOUFFLE_FUNCTIONS_H
