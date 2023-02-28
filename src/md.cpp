//
// Created by lauren on 6/7/18.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>
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
   double capsomere_concentration, ks, kb, number_capsomeres, ecut_c, elj_att;            // capsomere hamiltonian					
   double salt_concentration, temperature;	                                          // environmental or control parameters					
   double computationSteps, totaltime, delta_t, fric_zeta, chain_length_real, NListCutoff_c, NListCutoff, damp;// computational parameters
   bool verbose, restartFile, clusters;
   int buildFrequency, moviefreq, writefreq, restartfreq;
	
   double qs = 1;                                           //salt valency
   double T;                                                //set temperature (reduced units)
   double Q;                                           //nose hoover mass (reduced units)
	
   double const Avagadro = 6.022e23; // mol^-1		      //useful constants
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
   ("filename,f", value<string>(&file_name)->default_value("41part_c"), "Filename?")
   ("capsomere concentration,C", value<double>(&capsomere_concentration)->default_value(200),
   "capsomere concentration (micromolar)")                  // box size adjusteded for micromolar conc.
   ("salt concentration,c", value<double>(&salt_concentration)->default_value(500),
   "salt concentration (millimolar)")                       // electrostatic variables adjusted for millimolar conc.
   ("number of subunits,S", value<double>(&number_capsomeres)->default_value(64),
   "number of subunits")
   ("stretching constant,s", value<double>(&ks)->default_value(50), "stretching constant (KbT)")
   ("bending constant,b", value<double>(&kb)->default_value(20), "bending constant (KbT)")
   ("total time,T", value<double>(&computationSteps)->default_value(2500), "total time (computational steps)") // # of steps is total time / timestep
   ("timestep,t", value<double>(&delta_t)->default_value(0.004), "timestep (MD steps)")
   ("friction coefficient,r", value<double>(&fric_zeta)->default_value(1),
   "friction coefficient (reduced unit)")                   //used in brownian
   ("damping coefficient,d", value<double>(&damp)->default_value(100),"damping coefficient (unit of LJ time)")                   //used in brownian
   ("chain length,q", value<double>(&chain_length_real)->default_value(5), "nose hoover chain length") //used in MD
   ("thermal mass,Q", value<double>(&Q)->default_value(10), "nose hoover mass") //used in MD
   ("temperature,K", value<double>(&temperature)->default_value(298), "temperature (Kelvin)")
   ("ecut_c,e", value<double>(&ecut_c)->default_value(12), "electrostatics cutoff coefficient, input 0 for no cutoff")
   ("Restart bool,R", value<bool>(&restartFile)->default_value(true), "restartFile true: initializes from a restart file in outfiles/")
   ("Chunks bool,X", value<bool>(&clusters)->default_value(false), "clusters true: initializes from preformed clusters/")
   ("verbose,V", value<bool>(&verbose)->default_value(true), "verbose true: provides detailed output")
   ("lennard jones well depth,E", value<double>(&elj_att)->default_value(2), "lennard jones well depth")
   ("Neighbor list build frequency,B", value<int>(&buildFrequency)->default_value(20), "Neighbor list build frequency")
   ("Neighbor list cutoff,L", value<double>(&NListCutoff_c)->default_value(7.5), "Neighbor list cutoff (x + es & lj cutoff)")
   ("moviefreq,M", value<int>(&moviefreq)->default_value(1000), "The frequency of shooting the movie")
   ("writefreq,W", value<int>(&writefreq)->default_value(1000),"frequency of dumping energy file")
   ("restartfreq,w", value<int>(&restartfreq)->default_value(1000000), "The frequency of making restart files");
   
   variables_map vm;
   store(parse_command_line(argc, argv, desc), vm);
   notify(vm);
   if (vm.count("help")) {
         std::cout << desc << "\n";
         return 0;
   }
   
   ofstream traj("outfiles/energy.out", ios_base::app);     //setting up file outputs
   ofstream ofile("outfiles/ovito.lammpstrj", ios_base::app);
   ofstream sysdata("outfiles/model.parameters.out", ios::out);
   ofstream restart;
   string restartFilename;
   int restartStep;
   initialize_outputfile(traj, ofile);
   
   if (response == 'b') {                    		      //set flag for brownian vs. molecular dynamics
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
   double SImass;//SI value for a single bead (kg)
   double SIsigma;// (nm)
   double SItime;// (seconds)
   unsigned int cluster_size;
   
   vector<BEAD> subunit_bead;                               //Create particles, named subunit_bead
   vector<EDGE> subunit_edge;                               //create edges between subunit_bead's
   vector<SUBUNIT> protein;                                 //create subunits, named protein
   vector<FACE> subunit_face;                               //create faces of subunit_bead's on protein
   vector<THERMOSTAT> real_bath;                            //vector of thermostats
    
   vector<vector<int> > lj_a;
   lj_a = generate_lattice(capsomere_concentration, number_capsomeres, file_name, bondlength, SIsigma, SImass, 
                           subunit_bead, subunit_edge, protein, subunit_face, restartFile, restartStep, clusters, cluster_size);     //Setting up the input file (uses user specified file to generate lattice)
                        
   double SIenergy =  temperature * Boltzmann;// Joules
   SItime = sqrt(SIsigma*SIsigma*SImass/SIenergy);//seconds
   totaltime = computationSteps * delta_t;
   T = temperature / double(298);
   double UnitEnergy = SIenergy * SItime * SItime / (SImass * SIsigma * (1e-9));
   //cout << "Unit energy is " << UnitEnergy << endl;
                        
   if (world.rank() == 0) {                                 // making model.parameters file
      sysdata << "Simulation will run for " << totaltime * SItime / (1e-9) << " nanoseconds with a "
      << delta_t * SItime / (1e-12) << " picosecond timestep." << endl;
      sysdata << "Capsomere concentration: " << capsomere_concentration << " micromolar" << endl;
      sysdata << "Salt concentration: " << salt_concentration << " millimolar" << endl;
      sysdata << "Stretching constant: " << ks << " KbT" << "( " << (ks * SImass / (SItime * SItime)) << " N/m)" << endl;
      sysdata << "Bending constant: " << kb << " KbT" << "( " << (kb * SImass * SIsigma * SIsigma / (SItime * SItime)) << " J)" << endl;;
      sysdata << "Bondlength between beads is " << bondlength << " LJ reduced units, which is "
      << bondlength * SIsigma / (1e-9) << " nanometers." << endl;
      sysdata << "Mass of a bead is " << SImass << " kg." << endl;
      sysdata << "Diameter of a bead is " << SIsigma / (1e-9) << " nanometers." << endl;
      sysdata << "Total number of subunits is " << number_capsomeres << endl;
      sysdata << "Temperature is " << temperature << " K" << " Which is " << T << " in reduced units" << endl;
      sysdata << "Attractive LJ paramerter is " << elj_att << " Which is " << elj_att * (SIenergy * Avagadro / (1000 * UnitEnergy) ) << " kJ/mol." << endl;
   }
     
	// LJ features
   double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * Avagadro)),1.0 / 3.0); //calculating box size, prefactor of 1000 used to combine units
   VECTOR3D box_size = VECTOR3D(box_x, box_x, box_x);
   double ecut = 2.5 * (SIsigma / 1e-9);	                  // Lennard-Jones cut-off distance
	
	// Electrostatic features
   double lb = (0.701e-9) / SIsigma;   // e^2 / (4 pi Er E0 Kb T) ; value for T = 298 K.
   //number density (1/sigma*^3)
   double ni = salt_concentration * Avagadro * SIsigma * SIsigma * SIsigma;
   //electrostatics parameters
   double kappa = sqrt(8 * Pi * ni * lb * qs * qs);
   double screen = 1 / kappa;	
   double ecut_el = screen * ecut_c;			      // screening length times a constant so that electrostatics is cutoff at approximately 0.015
   if (ecut_el > ecut) NListCutoff = ecut_el + NListCutoff_c;
       else NListCutoff = ecut + NListCutoff_c;
   
   if (world.rank() == 0) {
      sysdata << "Box length is " << box_x * SIsigma / (1e-9) << " nanometers." << endl;
      sysdata << "Screening length is " << screen << " nanometers." << endl;
      sysdata << "electrostatic cutoff is " << ecut_el * SIsigma / (1e-9) << " nanometers." << endl;
      sysdata << "Neighborlist cutoff is " << NListCutoff * SIsigma / (1e-9) << " nanometers." << endl;
     }

   if (brownian == false) {                                 //for molecular, set up the nose hoover thermostat
      if (chain_length_real == 1)
         real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0, 1)));
      else {
         real_bath.push_back((THERMOSTAT(Q, T, 3 * subunit_bead.size(), 0, 0, 0, 1)));
         while (real_bath.size() != chain_length_real - 1)  // final bath is dummy bath (dummy bath always has zero mass)
            real_bath.push_back((THERMOSTAT(1, T, 1, 0, 0, 0, 1)));
         real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0, 1)));
		} //else
   }

   //MPI Boundary calculation for ions
   unsigned int rangeIons = (protein.size() + world.size() - 1) / (1.0 * world.size());
   lowerBound = world.rank() * rangeIons;
   upperBound = (world.rank() + 1) * rangeIons - 1;
  // extraElements = world.size() * rangeIons - protein.size();
   sizFVec = rangeIons; //upperBound - lowerBound + 1;
   if (world.rank() == world.size() - 1) {
      upperBound = protein.size() - 1;
    //  sizFVec = upperBound - lowerBound + 1;// + extraElements;
   }
//    if (world.size() == 1) {
//       lowerBound = 0;
//       upperBound = protein.size() - 1;
//    }

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
   /*									Initial Force Calculation				     */
   //////////////////////////////////////////////////////////////////////////////////////////////////////////

   bool updatePairlist = true;
   unsigned int NListVectorSize;//calculating max number of neighbors (conservative estimate assuming 100% packing efficiency)
   NListVectorSize = ceil( (NListCutoff + 0.5) * (NListCutoff + 0.5) * (NListCutoff + 0.5) ) ;
   if (NListVectorSize > subunit_bead.size()) NListVectorSize = subunit_bead.size();
   if (world.rank() == 0) {
      sysdata << "Neighborlist has a maximum of " << NListVectorSize << " neighbors per bead." << endl;
     }

   for (unsigned int i = 0; i < subunit_bead.size(); i ++) {
      subunit_bead[i].itsN.assign(NListVectorSize, -1);                       //Making "empty" pairlist (fill with -1)                    
   }

   forceCalculation(protein, lb, ni, qs, subunit_bead, ecut, ks, kb, lj_a, ecut_el, kappa, elj_att, updatePairlist, NListCutoff);

   double senergy = 0;                                                        //blank all the energy metrics
   double kenergy = 0;
   double benergy = 0;
   double tenergy = 0;
   double ljenergy = 0;
   double cenergy = 0;
   double tpenergy = 0;
   double tkenergy = 0;

   if (restartFile == false) {                                                //assign random velocities based on initial temperature
      initialize_bead_velocities(protein, subunit_bead, T, clusters, cluster_size);
     //initialize_constant_bead_velocities(protein, subunit_bead, T);
   }

   double particle_ke = particle_kinetic_energy(subunit_bead);                //thermostat variables for nose hoover
   double expfac_real = 0;


    /*
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);                               //setting up random seed for brownian
    unsigned long int Seed = 23410981;
    gsl_rng_set(r, Seed);
    if (world.rank() == 0) {
       ofile << "ITEM: TIMESTEP" << endl << 0 << endl << "ITEM: NUMBER OF ATOMS" << endl << subunit_bead.size()
       << endl
       << "ITEM: BOX BOUNDS" << endl << -box_size.x / 2 << setw(15) << box_size.x / 2 << endl << -box_size.y / 2
       << setw(15)
       << box_size.y / 2 << endl << -box_size.z / 2 << setw(15) \
       << box_size.z / 2 << endl << "ITEM: ATOMS index type x y z b charge" << endl;
    }
    if (world.rank() == 0) {
       for (unsigned int b = 0; b < subunit_bead.size(); b++) {
          ofile << b + 1 << setw(15) << subunit_bead[b].type << setw(15) << subunit_bead[b].pos.x << setw(15)
          << subunit_bead[b].pos.y << setw(15) << subunit_bead[b].pos.z << setw(15)
          << subunit_bead[b].be << setw(15) << subunit_bead[b].q << endl;
       }
    }*/
   
   

   /*                  ___                        __      __      ___
   *      /|   /|      |   \             |        /  \    /  \    |   \
   *     / |  / |      |    \            |       |    |  |    |   |    |
   *    /  | /  | ---- |     |           |       |    |  |    |   |___/
   *   /   |/   |      |    /            |       |    |  |    |   |
   *  /         |      |___/             |_____   \__/    \__/    |                       */

   int loopStart;
   if (restartFile == false) loopStart = 0;
   if (restartFile == true) loopStart = restartStep;


   //set up random distribution for brownian
   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator (seed);
   std::normal_distribution<double> distribution (0.0,1.0);
   std::uniform_real_distribution<> distr(-0.5,0.5);
                        
   for (unsigned int a = loopStart; a < ((totaltime / delta_t)+1); a++) {        // BEGIN MD LOOP
   
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*								VELOCITY VERLET		                                */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      //  Brownian Dynamics Equations
      //  r1 = -m / damp;
      //  r2 = sqrt((24*m*kB)/(damp*dt));
      //  frandom = r2 * uniform(-0.5, 0.5);
      //  fdrag = r1 * v
      // force = force + frandom + fdrag

      
      if (brownian == false) {                                                //FOR MOLECULAR DYNAMICS
         for (int i = real_bath.size() - 1; i > -1; i--)                      //thermostat update
            update_chain_xi(i, real_bath, delta_t, particle_ke);
         for (unsigned int i = 0; i < real_bath.size(); i++)
            real_bath[i].update_eta(delta_t);
			
         expfac_real = exp(-0.5 * delta_t * real_bath[0].xi);

                                 
         for (unsigned int i = 0; i < protein.size(); i++) {                  //velocity verlet loop
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);  //update velocity half step
               protein[i].itsB[ii]->update_position(delta_t);                                   //update position full step
            }
         } // for i
      } else {                                                                    // FOR BROWNIAN DYNAMICS
         for (int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){
               protein[i].itsB[ii]->compute_fdrag(damp);
            }
          }
         for (int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++){
               protein[i].itsB[ii]->update_velocity(delta_t);
               protein[i].itsB[ii]->update_position(delta_t);
            }
         }

      }  // else
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*                                                                      UPDATE PAIRLIST                 */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////  
      updatePairlist = false;
      if ( a % buildFrequency == 0) {
         updatePairlist = true;
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*									MD LOOP FORCES						  */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
      forceCalculation(protein, lb, ni, qs, subunit_bead, ecut, ks, kb, lj_a, ecut_el, kappa, elj_att, updatePairlist, NListCutoff);
            
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*								VELOCITY VERLET							  */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
      if (brownian == false) {                                                //FOR MOLECULAR DYNAMICS
         for (unsigned int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);  //update velocity the other half step
            }
         }
         particle_ke = particle_kinetic_energy(subunit_bead);
         for (unsigned int i = 0; i < real_bath.size(); i++)                  // forward update of Nose-Hoover chain
            real_bath[i].update_eta(delta_t);
         for (unsigned int i = 0; i < real_bath.size(); i++)
            update_chain_xi(i, real_bath, delta_t, particle_ke);
      } else {                                                                //FOR BROWNIAN DYNAMICS
         for (int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->compute_fran(delta_t, damp, distr(generator));
               //protein[i].itsB[ii]->update_tforce();
               protein[i].itsB[ii]->update_velocity(delta_t);
            }
         }

      }  // else

      /*     __                 __                        ____     ________     ____
      *     /  \    |\    |    /  \    |       \   /     /    \        |       /    \
      *    |____|   | \   |   |____|   |        \ /      \_            |       \_
      *    |    |   |  \  |   |    |   |         |         \__         |         \__
      *    |    |   |   \ |   |    |   |         |            \        |            \
      *    |    |   |    \|   |    |   |_____    |       \____/    ____|____   \____/                                */
   
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*                                                              MAKING RESTART FILE                                                                                                                */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////         
      if (a % restartfreq == 0 && world.rank() == 0) {
         stringstream step;
         step << a;
         restartFilename = "outfiles/restart_" + step.str();
         restartFilename += ".out";
         restart.open(restartFilename.c_str(), ios::out);
         restart << "Velocities & Positions for " << a << endl;
         for (unsigned int i = 0; i < subunit_bead.size(); i++) {
             restart << i << "  " << subunit_bead[i].vel.x << setw(25) << setprecision(12) << subunit_bead[i].vel.y << setw(25) << setprecision(12) << subunit_bead[i].vel.z  << setw(25) << setprecision(12) << subunit_bead[i].pos.x << setw(25) << setprecision(12) << subunit_bead[i].pos.y << setw(25) << setprecision(12) << subunit_bead[i].pos.z  << setw(25) << setprecision(12) << endl;
       //     restart << subunit_bead[i].tforce.x << setw(25) << setprecision(12) << subunit_bead[i].tforce.y << setw(25) << setprecision(12) << subunit_bead[i].tforce.z  << setw(25) << setprecision(12)
                //     << endl;
         }
         restart.close();
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*								ANALYZE ENERGIES							     */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
         
      if (a % writefreq == 0 && world.rank() == 0) {
   
         dress_up(subunit_edge, subunit_face);                     //update edge and face properties
   
         for (unsigned int i = 0; i < protein.size(); i++) {       //blanking out energies here
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->be = 0;
               protein[i].itsB[ii]->ne = 0;
               protein[i].itsB[ii]->ce = 0;
            }
         }
   
         for (unsigned int i = 0; i < protein.size(); i++) {       // Intramolecular Energies
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               protein[i].itsB[ii]->update_stretching_energy(ks);
               protein[i].itsB[ii]->update_kinetic_energy();
            }
            for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++) {
               if (protein[i].itsE[kk]->type != 0)                //if it is a bending edge...
                  protein[i].itsE[kk]->update_bending_energy(kb);
            }
         }
                                                                  //Intermolecular Energies
         update_LJ_ES_energies_simplified(subunit_bead, ecut, lj_a, elj_att, lb, ni, qs, ecut_el, kappa);
   
         for (unsigned int i = 0; i < protein.size(); i++) {      //blanking out energies here
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               senergy += protein[i].itsB[ii]->se;                //sum up total energies
               kenergy += protein[i].itsB[ii]->ke;
               ljenergy += protein[i].itsB[ii]->ne;
               benergy += protein[i].itsB[ii]->be;
               cenergy += protein[i].itsB[ii]->ce;
            }
         }
   
         for (unsigned int i = 0; i < real_bath.size(); i++) {    //thermostat energies
            real_bath[i].potential_energy();
            real_bath[i].kinetic_energy();
         }
         for (unsigned int i = 0; i < real_bath.size(); i++) {
            tpenergy += real_bath[i].pe;                          //sum up total energies
            tkenergy += real_bath[i].ke;
         }
   
         tenergy = senergy + kenergy + ljenergy + benergy + cenergy + tpenergy +
         tkenergy;                                                
         
         if (world.rank() == 0) {                                 //print info to files for data analysis
            traj << a << setw(15) << kenergy / subunit_bead.size() << setw(15)
                 << senergy / subunit_bead.size() << setw(15) << benergy / subunit_bead.size() 
                 << setw(15) << ljenergy / subunit_bead.size() << setw(15) << cenergy / subunit_bead.size()
                 << setw(15) << tenergy / subunit_bead.size() << setw(15) 
                 << (benergy + senergy + ljenergy + cenergy) / subunit_bead.size() << setw(15)
                 << kenergy * 2 / (3 * subunit_bead.size()) << setw(15) << tpenergy / subunit_bead.size()
                 << setw(15) << tkenergy / subunit_bead.size() << endl;
         }
         senergy = 0;                                             //blanking out energies
         kenergy = 0;
         benergy = 0;
         tenergy = 0;
         ljenergy = 0;
         cenergy = 0;
         tpenergy = 0;
         tkenergy = 0;
      } // energy analysis loop
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*								STORE POSITION INFO TO FILE					     */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////
      if (a % moviefreq == 0 && world.rank() == 0) {
         if (world.rank() == 0) {
            ofile << "ITEM: TIMESTEP" << endl << a << endl << "ITEM: NUMBER OF ATOMS" << endl << subunit_bead.size()
            << endl
            << "ITEM: BOX BOUNDS" << endl << -box_size.x / 2 << setw(15) << box_size.x / 2 << endl << -box_size.y / 2
            << setw(15)
            << box_size.y / 2 << endl << -box_size.z / 2 << setw(15) \
            << box_size.z / 2 << endl << "ITEM: ATOMS index type x y z b charge" << endl;
         }
         if (world.rank() == 0) {
            for (unsigned int b = 0; b < subunit_bead.size(); b++) {
               ofile << b + 1 << setw(15) << subunit_bead[b].type << setw(15) << subunit_bead[b].pos.x << setw(15)
               << subunit_bead[b].pos.y << setw(15) << subunit_bead[b].pos.z << setw(15)
               << subunit_bead[b].be << setw(15) << subunit_bead[b].q << endl;
            }
         }
      } //position print loop
      
      if (world.rank() == 0) {                                    // PRINT PROGRESS BAR
         percentage = roundf(a / (totaltime / delta_t) * 100 * 10) / 10;
         if (percentage != percentagePre) {
            double fraction_completed = percentage / 100;
            ProgressBar(fraction_completed);
            percentagePre = percentage;
         }
      }
   } //time loop end

//gsl_rng_free (r);                            //free gsl memory from brownian random variables
   
   return 0;
} //end simulation fxn


/* old brownian function
 *          double c2;
         for (unsigned int i = 0; i < subunit_bead.size(); i++) {
            c2 = 1 + (delta_t * 0.5 * fric_zeta / subunit_bead[i].m);
            subunit_bead[i].noise.x = gsl_ran_gaussian(r,1) * sqrt(2 * UnitEnergy * fric_zeta * delta_t);                  //determine noise term
            subunit_bead[i].noise.y = gsl_ran_gaussian(r,1) * sqrt(2 * UnitEnergy * fric_zeta * delta_t);
            subunit_bead[i].noise.z = gsl_ran_gaussian(r,1) * sqrt(2 * UnitEnergy * fric_zeta * delta_t);
            subunit_bead[i].oldtforce = subunit_bead[i].tforce;
            //subunit_bead[i].vel += (subunit_bead[i].tforce ^ (0.5 * delta_t / subunit_bead[i].m)) + (subunit_bead[i].pos ^ (fric_zeta / subunit_bead[i].m)) + subunit_bead[i].noise;
           // cout << "velocities are " << subunit_bead[i].vel.x << " , " << subunit_bead[i].vel.y << " , " << subunit_bead[i].vel.z << endl;
         }
         for (unsigned int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
               c2 = 1 + (delta_t * 0.5 * fric_zeta / protein[i].itsB[ii]->m);
               protein[i].itsB[ii]->update_position_brownian(delta_t, c2, fric_zeta);                                   //update position full step
            }  // for ii
         }  // for i


         for (unsigned int i = 0; i < subunit_bead.size(); i++) {
            double c2 = 1 + (delta_t * 0.5 * fric_zeta / subunit_bead[i].m);
            //subunit_bead[i].vel += (subunit_bead[i].tforce ^ (0.5 * delta_t / subunit_bead[i].m)) + (subunit_bead[i].pos ^ (fric_zeta / subunit_bead[i].m)) + subunit_bead[i].noise;
             VECTOR3D r_vec; //= (A->pos - B->pos);
             r_vec.x = subunit_bead[i].oldpos.x - subunit_bead[i].pos.x;
             r_vec.y = subunit_bead[i].oldpos.y - subunit_bead[i].pos.y;
             r_vec.z = subunit_bead[i].oldpos.z - subunit_bead[i].pos.z;
             VECTOR3D box = subunit_bead[i].bx;
             VECTOR3D hbox = subunit_bead[i].hbx;
             if (r_vec.x > hbox.x) r_vec.x -= box.x;
             else if (r_vec.x < -hbox.x) r_vec.x += box.x;
             if (r_vec.y > hbox.y) r_vec.y -= box.y;
             else if (r_vec.y < -hbox.y) r_vec.y += box.y;
             if (r_vec.z > hbox.z) r_vec.z -= box.z;
             else if (r_vec.z < -hbox.z) r_vec.z += box.z;
             subunit_bead[i].vel += ((subunit_bead[i].oldtforce + subunit_bead[i].tforce) ^ (0.5 * delta_t / subunit_bead[i].m)) + (subunit_bead[i].noise ^ (1/subunit_bead[i].m)) + (r_vec ^ (fric_zeta / subunit_bead[i].m));
         }

*/
