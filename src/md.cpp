//
// Created by lauren on 6/7/18.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "initialize.h"
#include "functions.h"
#include "energies.h"
#include "forces.h"
#include "md.h"

using namespace std;
using namespace boost::program_options;

//MPI boundary parameters
unsigned int lowerBound;
unsigned int upperBound;
unsigned int sizFVec;
unsigned int extraElements;

mpi::environment env;
mpi::communicator world;

int run_simulation(int argc, char *argv[]) {
   
   bool brownian = true;
                                                            //parameters from user
   char response;
   string file_name;
   double capsomere_concentration, ks, kb, number_capsomeres, ecut_c, elj_att;	// capsomere hamiltonian					
   double salt_concentration, temperature;	// environmental or control parameters					
   double totaltime, delta_t, fric_zeta, chain_length_real, NListCutoff;	// computational parameters
   bool verbose, restartFile;
   int buildFrequency ;
	
	double qs = 1;                                           //salt valency
   double T = 1;                                            //set temperature (reduced units)
   double Q = 10;                                           //nose hoover mass (reduced units)
	
	double const Avagadro = 6.022e23; // mol^-1					//useful constants
   double const Boltzmann = 1.3806e-23; // m2kg/s2K
   //double const e0 = 8.854187e-12; // C2/Nm2
   //double const q_electron = 1.602e-19; // C
   double const Pi = 3.14159;
                                                            //Progress bar paras
   double percentage = 0, percentagePre = -1;
                                                            // Get input values from the user
   options_description desc("Usage:\nrandom_mesh <options>");
   desc.add_options()
   ("help,h", "print usage message")
   ("Dynamics Response,D", value<char>(&response)->default_value('m'),
   "To run brownian dynamics (overdamped langevin) enter 'b'. Otherwise, to run molecular dynamics with nose' hoover thermostat enter 'm'. [b/m]")
   ("filename,f", value<string>(&file_name)->default_value("41part"), "Filename?")
   ("capsomere concentration,C", value<double>(&capsomere_concentration)->default_value(75),
   "capsomere concentration (micromolar)") // box size adjusteded for micromolar conc.
   ("salt concentration,c", value<double>(&salt_concentration)->default_value(200),
   "salt concentration (millimolar)") // electrostatic variables adjusted for millimolar conc.
   ("number of subunits,S", value<double>(&number_capsomeres)->default_value(8),
   "number of subunits")
   ("stretching constant,s", value<double>(&ks)->default_value(50), "stretching constant (KbT)")
   ("bending constant,b", value<double>(&kb)->default_value(20), "bending constant (KbT)")
   ("total time,T", value<double>(&totaltime)->default_value(100), "total time (MD steps)") // # of steps is total time / timestep
   ("timestep,t", value<double>(&delta_t)->default_value(0.001), "timestep (MD steps)")
   ("friction coefficient,r", value<double>(&fric_zeta)->default_value(1),
   "friction coefficient (reduced unit)") //used in brownian
   ("chain length,q", value<double>(&chain_length_real)->default_value(5), "nose hoover chain length") //used in MD
   ("temperature,K", value<double>(&temperature)->default_value(298), "temperature (Kelvin)")
   ("ecut_c,e", value<double>(&ecut_c)->default_value(20), "electrostatics cutoff coefficient, input 0 for no cutoff")
   ("Restart bool,R", value<bool>(&restartFile)->default_value(false), "restartFile true: initializes from a restart file in outfiles/")
   ("verbose,V", value<bool>(&verbose)->default_value(true), "verbose true: provides detailed output")
   ("lennard jones well depth,E", value<double>(&elj_att)->default_value(2), "lennard jones well depth")
   ("Neighbor list build frequency,B", value<int>(&buildFrequency)->default_value(20), "Neighbor list build frequency")
   ("Neighbor list cutoff,L", value<double>(&NListCutoff)->default_value(4.0), "Neighbor list cutoff");
   
   variables_map vm;
   store(parse_command_line(argc, argv, desc), vm);
   notify(vm);
   if (vm.count("help")) {
         std::cout << desc << "\n";
         return 0;
   }
   
   ofstream traj("outfiles/energy.out", ios_base::app);              //setting up file outputs
   ofstream ofile("outfiles/ovito.lammpstrj", ios_base::app);
   ofstream sysdata("outfiles/model.parameters.out", ios::out);
   ofstream restart;
   int restartStep;
   initialize_outputfile(traj, ofile);
   
   if (response == 'b') {                    			//set flag for brownian vs. molecular dynamics
      brownian = true;
      if (world.rank() == 0)
         sysdata << "Running brownian dynamics." << endl;
   } else if (response == 'm') {
      brownian = false;
      if (world.rank() == 0)
         sysdata << "Running molecular dynamics with Nose Hoover thermostat." << endl;
   }
																		//System specific paramters (modelled for HBV)
   double bondlength;     
   double SImass;                             			//SI value for a single bead (kg)
   double SIsigma;                                    // (nm)
   double SItime;                                     // (seconds)
   
   vector<BEAD> subunit_bead;                              //Create particles, named subunit_bead
   vector<EDGE> subunit_edge;                              //create edges between subunit_bead's
   vector<SUBUNIT> protein;                                //create subunits, named protein
   vector<FACE> subunit_face;                              //create faces of subunit_bead's on protein
   vector<PAIR> lj_pairlist;                               //create vector to hold LJ pairings
   vector<THERMOSTAT> real_bath;                           //vector of thermostats
    
   vector<vector<int> > lj_a;
   lj_a = generate_lattice(capsomere_concentration, number_capsomeres, file_name, bondlength, SIsigma, SImass, 
                        subunit_bead, subunit_edge, protein, subunit_face, restartFile, restartStep);     //Setting up the input file (uses user specified file to generate lattice)
                        
   double SIenergy =  temperature * Boltzmann; 		// Joules
   SItime = sqrt(SIsigma*SIsigma*SImass/SIenergy); 	//seconds
                        
   if (world.rank() == 0) {
      sysdata << "Simulation will run for " << totaltime * SItime / (1e-9) << " nanoseconds with a "
      << delta_t * SItime / (1e-12) << " picosecond timestep." << endl;
      sysdata << "Capsomere concentration: " << capsomere_concentration << " micromolar" << endl;
      sysdata << "Salt concentration: " << salt_concentration << " millimolar" << endl;
      sysdata << "Stretching constant: " << ks << " KbT" << endl;
      sysdata << "Bending constant: " << kb << " KbT" << endl;
      sysdata << "Bondlength between beads is " << bondlength << " LJ reduced units, which is "
      << bondlength * SIsigma / (1e-9) << " nanometers." << endl;
      sysdata << "Mass of a bead is " << SImass << " kg." << endl;
      sysdata << "Diameter of a bead is " << SIsigma / (1e-9) << " nanometers." << endl;
      sysdata << "Total number of subunits is " << number_capsomeres << endl;
      sysdata << "Temperature is " << temperature << " K" << endl;
     }
     
	// LJ features
   double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * Avagadro)),1.0 / 3.0);    //calculating box size, prefactor of 1000 used to combine units
   VECTOR3D box_size = VECTOR3D(box_x, box_x, box_x);
   double ecut = 2.5 * (SIsigma / 1e-9);	// Lennard-Jones cut-off distance
	
	// Electrostatic features
   double lb = (0.701e-9) / SIsigma;   // e^2 / (4 pi Er E0 Kb T) ; value for T = 298 K.
   //number density (1/sigma*^3)
   double ni = salt_concentration * Avagadro * SIsigma * SIsigma * SIsigma;
   //electrostatics parameter
   double kappa = sqrt(8 * Pi * ni * lb * qs * qs);
   double screen = 1 / kappa;							 	 // 
   double ecut_el = screen * ecut_c;			// screening length times a constant so that electrostatics is cutoff at approximately 0.015
   
   if (world.rank() == 0) {
      sysdata << "Box length is " << box_x * SIsigma / (1e-9) << " nanometers." << endl;
      sysdata << "Screening length is " << screen << " nanometers." << endl;
     }

   if (brownian == false) {             //for molecular, set up the nose hoover thermostat
      if (chain_length_real == 1)
         real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0, 1)));
      else {
         real_bath.push_back((THERMOSTAT(Q, T, 3 * subunit_bead.size(), 0, 0, 0, 1)));
         while (real_bath.size() != chain_length_real - 1)
            real_bath.push_back((THERMOSTAT(1, T, 1, 0, 0, 0, 1)));
         real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0, 1)));
        // final bath is dummy bath (dummy bath always has zero mass)
		}
   }

   dress_up(subunit_edge, subunit_face);	// Calculate initial forces

   //MPI Boundary calculation for ions
   unsigned int rangeIons = subunit_bead.size() / world.size() + 1.5;
   lowerBound = world.rank() * rangeIons;
   upperBound = (world.rank() + 1) * rangeIons - 1;
   extraElements = world.size() * rangeIons - subunit_bead.size();
   sizFVec = upperBound - lowerBound + 1;
   if (world.rank() == world.size() - 1) {
      upperBound = subunit_bead.size() - 1;
      sizFVec = upperBound - lowerBound + 1 + extraElements;
   }
   if (world.size() == 1) {
      lowerBound = 0;
      upperBound = subunit_bead.size() - 1;
   }

   int numOfNodes = world.size();
   if (world.rank() == 0) {
   #pragma omp parallel default(shared)
   {
      if (omp_get_thread_num() == 0) {
         printf("The app comes with MPI and OpenMP (Hybrid) parallelization)\n");
         printf("Number of MPI processes used %d\n", numOfNodes);
         printf("Number of OpenMP threads per MPI process %d\n", omp_get_num_threads());
         printf("Make sure that number of beads is greater than %d\n", omp_get_num_threads() * numOfNodes);
         }
      }
   }
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*									Initial Force Calculation								*/
   //////////////////////////////////////////////////////////////////////////////////////////////////////////

   bool updatePairlist = true;

   int NListVectorSize;
                                                      //calculating max number of neighbors (conservative estimate assuming 100% packing efficiency)
   NListVectorSize = ceil( (NListCutoff + 0.5) * (NListCutoff + 0.5) * (NListCutoff + 0.5) ) ;

   for (unsigned int i = 0; i < subunit_bead.size(); i ++) {
      subunit_bead[i].itsN.assign(NListVectorSize, -1);                       //Making "empty" pairlist (fill with -1)
     //subunit_bead[i].itsN.assign(subunit_bead.size(), -1);                    
   }

   forceCalculation(protein, lb, ni, qs, subunit_bead, lj_pairlist, ecut, ks, bondlength, kb, lj_a, ecut_el, kappa, elj_att, updatePairlist, NListCutoff);

   double senergy = 0;                                                        //blank all the energy metrics
   double kenergy = 0;
   double benergy = 0;
   double tenergy = 0;
   double ljenergy = 0;
   double cenergy = 0;
   double tpenergy = 0;
   double tkenergy = 0;

   if (restartFile == false) {                                               //assign random velocities based on initial temperature
      initialize_bead_velocities(protein, subunit_bead, T);
     //initialize_constant_bead_velocities(protein, subunit_bead, T);
     }

                                                                           //thermostat variables for nose hoover
   double particle_ke = particle_kinetic_energy(subunit_bead);
   double expfac_real;                 //= exp(-0.5 * delta_t * real_bath[0].xi);
                        
                                                                           //setting up random seed for brownian
   gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
   unsigned long int Seed = 23410981;
   gsl_rng_set(r, Seed);

                        
                        
   /*                  ___                        __      __      ___
   *      /|   /|      |   \             |        /  \    /  \    |   \
   *     / |  / |      |    \            |       |    |  |    |   |    |
   *    /  | /  | ---- |     |           |       |    |  |    |   |___/
   *   /   |/   |      |    /            |       |    |  |    |   |
   *  /         |      |___/             |_____   \__/    \__/    |                       */

   int loopStart;
   if (restartFile == false) loopStart = 0;
   if (restartFile == true) loopStart = restartStep;
                        
   for (unsigned int a = loopStart; a < (totaltime / delta_t); a++) {        // BEGIN MD LOOP
   
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*								VELOCITY VERLET															*/
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                                 
      if (brownian == false) {   //FOR MOLECULAR DYNAMICS
         for (int i = real_bath.size() - 1; i > -1; i--)                    //thermostat update
            update_chain_xi(i, real_bath, delta_t, particle_ke);
         for (unsigned int i = 0; i < real_bath.size(); i++)
            real_bath[i].update_eta(delta_t);
                                 
         expfac_real = exp(-0.25 * delta_t * real_bath[0].xi);
                                 
         for (unsigned int i = 0; i < protein.size(); i++) {               //velocity verlet loop
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);  //update velocity half step
               protein[i].itsB[ii]->update_position(delta_t);                                   //update position full step
            }
         } // for i
      } else {            // FOR BROWNIAN DYNAMICS
         for (unsigned int i = 0; i < subunit_bead.size(); i++) {
            subunit_bead[i].vel.x +=
            (subunit_bead[i].vel.x * (-0.5 * fric_zeta * delta_t)) +
            (subunit_bead[i].tforce.x * (0.5 * delta_t / subunit_bead[i].m)) +
            sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) * (gsl_rng_uniform(r) - 0.5);
            subunit_bead[i].vel.y +=
            (subunit_bead[i].vel.y * (-0.5 * fric_zeta * delta_t)) +
            (subunit_bead[i].tforce.y * (0.5 * delta_t / subunit_bead[i].m)) +
            sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) * (gsl_rng_uniform(r) - 0.5);
            subunit_bead[i].vel.z +=
            (subunit_bead[i].vel.z * (-0.5 * fric_zeta * delta_t)) +
            (subunit_bead[i].tforce.z * (0.5 * delta_t / subunit_bead[i].m)) +
            sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) * (gsl_rng_uniform(r) - 0.5);
         }
         for (unsigned int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->update_position(delta_t);  //update position full step
            }  // for ii
         }  // for i
      }  // else
                        
      dress_up(subunit_edge, subunit_face);                              //update edge and face properties
                        
                        
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*                                                                      UPDATE PAIRLIST                                                                                          */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////  
                        
      // VECTOR3D r_vec = VECTOR3D(0, 0, 0);
      // long double r2 = 0.0;
      // VECTOR3D box = subunit_bead[0].bx;
      // double hbox = box.x / 2;
      updatePairlist = false;

      if ( a % buildFrequency == 0) {
         updatePairlist = true;
      }

      // if (updatePairlist == true) {
      //    for (unsigned int i = 0; i < subunit_bead.size(); i++) {
      //       fill(subunit_bead[i].itsN.begin(), subunit_bead[i].itsN.end(), -1);  //clear the pairlist to -1 (a number that cannot be bead index)
      //       int test = 0;
      //       for (unsigned int j = 0; j < subunit_bead.size(); j++) {
      //          r_vec = subunit_bead[i].pos - subunit_bead[j].pos;
      //          if (r_vec.x > hbox) r_vec.x -= box.x;
      //          else if (r_vec.x < -hbox) r_vec.x += box.x;
      //          if (r_vec.y > hbox) r_vec.y -= box.y;
      //          else if (r_vec.y < -hbox) r_vec.y += box.y;
      //          if (r_vec.z > hbox) r_vec.z -= box.z;
      //          else if (r_vec.z < -hbox) r_vec.z += box.z;
      //          r2 = r_vec.GetMagnitudeSquared();
      //          if (i != j && r2 < (4*4)) {
      //             subunit_bead[i].itsN[test] = subunit_bead[j].id;
      //             test += 1;
      //           //  cout << "here";
      //          }
      //       }// for j
      //      if(test > subunit_bead[0].itsN.size() ) cout << "ERROR! NEIGHBORLIST OUTGREW ALLOCATED VECTOR SIZE!" << endl;
      //    } // for i
      // } //if
                        
                        
                        
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*									MD LOOP FORCES												*/
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
                        
                        
      forceCalculation(protein, lb, ni, qs, subunit_bead, lj_pairlist, ecut, ks, bondlength, kb, lj_a, ecut_el, kappa, elj_att, updatePairlist, NListCutoff);
                        
                        
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*								VELOCITY VERLET															*/
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
                        
      if (brownian == false) {         //FOR MOLECULAR DYNAMICS
         for (unsigned int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               //update velocity the other half step
               protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);
            }
         }

         particle_ke = particle_kinetic_energy(subunit_bead);
                              // forward update of Nose-Hoover chain
         for (unsigned int i = 0; i < real_bath.size(); i++)
            real_bath[i].update_eta(delta_t);
         for (unsigned int i = 0; i < real_bath.size(); i++)
            update_chain_xi(i, real_bath, delta_t, particle_ke);
      } else {        //FOR BROWNIAN DYNAMICS
         for (unsigned int i = 0; i < subunit_bead.size(); i++) {
            subunit_bead[i].vel.x += (subunit_bead[i].vel.x * (-0.5 * fric_zeta * delta_t)) +
            (subunit_bead[i].tforce.x * (0.5 * delta_t / subunit_bead[i].m)) +
            sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) *
            (gsl_rng_uniform(r) - 0.5);
            subunit_bead[i].vel.y += (subunit_bead[i].vel.y * (-0.5 * fric_zeta * delta_t)) +
            (subunit_bead[i].tforce.y * (0.5 * delta_t / subunit_bead[i].m)) +
            sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) *
            (gsl_rng_uniform(r) - 0.5);
            subunit_bead[i].vel.z += (subunit_bead[i].vel.z * (-0.5 * fric_zeta * delta_t)) +
            (subunit_bead[i].tforce.z * (0.5 * delta_t / subunit_bead[i].m)) +
            sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) *
            (gsl_rng_uniform(r) - 0.5);
         }  
      }  // else

      /*      __                 __                        ____     ________     ____
      *     /  \    |\    |    /  \    |       \   /     /    \        |       /    \
      *    |____|   | \   |   |____|   |        \ /      \_            |       \_
      *    |    |   |  \  |   |    |   |         |         \__         |         \__
      *    |    |   |   \ |   |    |   |         |            \        |            \
      *    |    |   |    \|   |    |   |_____    |       \____/    ____|____   \____/                                */

      if (a % 1000 == 0) {                                           //analysis loop
   
         //////////////////////////////////////////////////////////////////////////////////////////////////////////
         /*                                                              MAKING RESTART FILE                                                                                                                */
         //////////////////////////////////////////////////////////////////////////////////////////////////////////         
   
         if (a % 1000 == 0 && world.rank() == 0) {
            restart.open("outfiles/forcescatter.out", ofstream::out | ofstream::trunc);
          //  restart << "Velocities & Positions for " << a << endl;
            for (unsigned int i = 0; i < subunit_bead.size(); i++) {
              // restart << i << "  " << subunit_bead[i].vel.x << setw(25) << setprecision(12) << subunit_bead[i].vel.y << setw(25) << setprecision(12) << subunit_bead[i].vel.z  << setw(25) << setprecision(12)
                //       << subunit_bead[i].pos.x << setw(25) << setprecision(12) << subunit_bead[i].pos.y << setw(25) << setprecision(12) << subunit_bead[i].pos.z  << setw(25) << setprecision(12) << endl;
               restart << subunit_bead[i].tforce.x << setw(25) << setprecision(12) << subunit_bead[i].tforce.y << setw(25) << setprecision(12) << subunit_bead[i].tforce.z  << setw(25) << setprecision(12)
                       << endl;
            }
            restart.close();
         }
   
   
         //////////////////////////////////////////////////////////////////////////////////////////////////////////
         /*								ANALYZE ENERGIES														*/
         //////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   
         for (unsigned int i = 0; i < protein.size(); i++) {       //blanking out energies here
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->be = 0;
               protein[i].itsB[ii]->ne = 0;
               protein[i].itsB[ii]->ce = 0;
            }
         }
   
   
         for (unsigned int i = 0; i < protein.size(); i++) {       // Intramolecular Energies
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->update_stretching_energy(ks, bondlength);
               protein[i].itsB[ii]->update_kinetic_energy();
            }
            for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++) {
               if (protein[i].itsE[kk]->type != 0)            //if it is a bending edge...
                  protein[i].itsE[kk]->update_bending_energy(kb);
            }
         }
                                                               //Intermolecular Energies
         update_ES_energies_simplified(subunit_bead, lb, ni, qs, ecut_el, kappa);
         update_LJ_energies_simplified(subunit_bead, ecut, lj_a, elj_att);
   
         for (unsigned int i = 0; i < protein.size(); i++) {       //blanking out energies here
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               senergy += protein[i].itsB[ii]->se;                          //sum up total energies
               kenergy += protein[i].itsB[ii]->ke;
               ljenergy += protein[i].itsB[ii]->ne;
               benergy += protein[i].itsB[ii]->be;
               cenergy += protein[i].itsB[ii]->ce;
            }
         }
   
   
         for (unsigned int i = 0; i < real_bath.size(); i++) {        //thermostat energies
            real_bath[i].potential_energy();
            real_bath[i].kinetic_energy();
         }
         for (unsigned int i = 0; i < real_bath.size(); i++) {
            tpenergy += real_bath[i].pe;                        //sum up total energies
            tkenergy += real_bath[i].ke;
         }
   
         tenergy = senergy + kenergy + ljenergy + benergy + cenergy + tpenergy +
         tkenergy;      //print info to files for data analysis
         //////////////////////////////////////////////////////////////////////////////////////////////////////////
         /*								STORE ENERGY INFO TO FILE												*/
         //////////////////////////////////////////////////////////////////////////////////////////////////////////
   
         if (world.rank() == 0) {
            ofile << "ITEM: TIMESTEP" << endl << a << endl << "ITEM: NUMBER OF ATOMS" << endl << subunit_bead.size()
            << endl
            << "ITEM: BOX BOUNDS" << endl << -box_size.x / 2 << setw(15) << box_size.x / 2 << endl << -box_size.y / 2
            << setw(15)
            << box_size.y / 2 << endl << -box_size.z / 2 << setw(15) \
            << box_size.z / 2 << endl << "ITEM: ATOMS index type x y z b charge" << endl;
            
            
            traj << a * delta_t << setw(15) << kenergy / subunit_bead.size() << setw(15)
            << senergy / subunit_bead.size() << setw(15) <<
            benergy / subunit_bead.size() << setw(15) << ljenergy / subunit_bead.size() << setw(15)
            << cenergy / subunit_bead.size()
            << setw(15) << tenergy / subunit_bead.size()
            << setw(15) << (benergy + senergy + ljenergy + cenergy) / subunit_bead.size() << setw(15)
            << kenergy * 2 / (3 * subunit_bead.size()) << setw(15) << tpenergy / subunit_bead.size()
            << setw(15)
            << tkenergy / subunit_bead.size() << endl;
         }
         if (world.rank() == 0) {
            for (unsigned int b = 0; b < subunit_bead.size(); b++) {
               ofile << b + 1 << setw(15) << subunit_bead[b].type << setw(15) << subunit_bead[b].pos.x << setw(15)
               << subunit_bead[b].pos.y << setw(15) << subunit_bead[b].pos.z << setw(15)
               << subunit_bead[b].be << setw(15) << subunit_bead[b].q << endl;
            }
         }
         senergy = 0;                            //blanking out energies
         kenergy = 0;
         benergy = 0;
         tenergy = 0;
         ljenergy = 0;
         cenergy = 0;
         tpenergy = 0;
         tkenergy = 0;
      } // end of energy analysis loop


      if (world.rank() == 0) {
         percentage = roundf(a / (totaltime / delta_t) * 100 * 10) / 10;
         //percentage output
         if (percentage != percentagePre) {
            double fraction_completed = percentage / 100;
            ProgressBar(fraction_completed);
            percentagePre = percentage;
         }
      }
   } //time loop end


   //  compute_MD_trust_factor_R(1);                   //computes R
   gsl_rng_free (r);
   
   return 0;
} //end simulation fxn


