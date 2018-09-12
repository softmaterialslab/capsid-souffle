// This is header file for the UTILITY class.
// This file includes standard library files and gsl functions that are utilized in the code
// This file also has useful constant parameters for the problem at hand

#ifndef _UTILITY_H
#define _UTILITY_H

//OPENMP
#include <omp.h>
//BOOST MPI
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/program_options.hpp>

using namespace boost::program_options;

namespace mpi = boost::mpi;

extern mpi::environment env;
extern mpi::communicator world;

//MPI boundary parameters
extern unsigned int lowerBound;
extern unsigned int upperBound;
extern unsigned int sizFVec;
extern unsigned int extraElements;

#endif
