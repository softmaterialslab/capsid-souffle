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
    double capsomere_concentration, salt_concentration, ks, kb, number_capsomeres, totaltime, delta_t, fric_zeta, chain_length_real, temperature, ecut_c, elj_att;					
    bool verbose, restartFile;

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
	    ("lennard jones well depth,E", value<double>(&elj_att)->default_value(2), "lennard jones well depth");


    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }
    
    
    ofstream traj("outfiles/energy.out", ios_base::app);              //setting up file outputs
    ofstream ofile("outfiles/ovito.lammpstrj", ios_base::app);
  //  ofstream msdata("outfiles/ms.out", ios_base::app);
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

    double bondlength;//= pow(2, 0.166666666667);       	//System specific paramters (modelled for HBV)
    double SImass;                             			//SI value for a single bead (kg)
    double SIsigma;                                             // (nm)
    double SItime;                                              // (seconds)
    double const Avagadro = 6.022e23; // mol^-1			//useful constants
    double const Boltzmann = 1.3806e-23; // m2kg/s2K
    //double const e0 = 8.854187e-12; // c2/Nm2
    //double const q_electron = 1.602e-19; // C
    double const Pi = 3.14159;

    vector<BEAD> subunit_bead;                              //Create particles, named subunit_bead
    vector<EDGE> subunit_edge;                              //create edges between subunit_bead's
    vector<SUBUNIT> protein;                                //create subunits, named protein
    vector<FACE> subunit_face;                              //create faces of subunit_bead's on protein
    vector<PAIR> lj_pairlist;                               //create vector to hold LJ pairings
    vector<THERMOSTAT> real_bath;                           //vector of thermostats

    
    double qs = 1;                                           //salt valency
    double T = 1;                                            //set temperature (reduced units)
    double Q = 10;                                           //nose hoover mass (reduced units)


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
    double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * Avagadro)),
                       1.0 / 3.0);    							//calculating box size, prefactor of 1000 used to combine units
    VECTOR3D bxsz = VECTOR3D(box_x, box_x, box_x);
    double lb = (0.701e-9) / SIsigma; // at 300 K only!!!                            	// e^2 / (4 pi Er E0 Kb T)
    double ni = salt_concentration * Avagadro * SIsigma * SIsigma * SIsigma;      	//number density (1/sigma*^3)
    double kappa = sqrt(8 * Pi * ni * lb * qs * qs);					//electrostatics parameter
    double screen = 1 / kappa;							 	// 
    double ecut_el = screen * ecut_c;							// screening length times a constant so that electrostatics is cutoff at approximately 0.015
    double ecut = 2.5 * (SIsigma / 1e-9);                                               //lennard jones cut-off distance
    unsigned int M = floor( bxsz.x / ecut );
    double cell_x = bxsz.x / M;                                                          //Cell width for neighborlist
    int head [M][M][M];
    if (world.rank() == 0) {
       cout << "GOT HERE 1!" << endl;
    }
    vector<int> list(subunit_bead.size());
    
    /*
    int i_x, j_y, k_z; //index in x y and z directions
    int icell; //index of cell
    double r_cut = ecut*2;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Set up cell-linked neighborlist
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (unsigned int i = 0; i < M; i++) {
       for (unsigned int j = 0; j < M; j++) {
          for (unsigned int k = 0; k < M; k++) {
             head[i][j][k] = -1;
         }
      }
   }
   for (unsigned int i = 0; i < list.size(); i++) {
      list[i] = -1;
   }
    for (unsigned int i = 0; i < subunit_bead.size(); i++) { //put beads into cells
       i_x = floor(subunit_bead[i].pos.x / cell_x);
       j_y = floor(subunit_bead[i].pos.y / cell_x);
       k_z = floor(subunit_bead[i].pos.z / cell_x);
       //icell = ( (i_x * M * M) + (j_y * M) + k_z );
       list[i] = head[i_x][j_y][k_z];
       head[i_x][j_y][k_z] = i;
   }
   
   if (world.rank() == 0) {
      cout << "GOT HERE!" << endl;
   }
   
   int i_y, i_z;
   int j;
   VECTOR3D r_vec;
   double r_dist;
   VECTOR3D box = subunit_bead[0].bx;
   
   for (unsigned int i = 0; i < subunit_bead.size(); i++){ // for all particles...
      for (int n = -1; n < 2; n++){ //look in all of the neighboring cells (27)
         for (int m = -1; m < 2; m++){
            for(int k = -1; k < 2; k++){
               i_x = subunit_bead[i].cell_id.x + n;
               i_y = subunit_bead[i].cell_id.y + m;
               i_z = subunit_bead[i].cell_id.z + k;
               if (i_x < 0) i_x = M;    //account for periodic boundaries
               if (i_x > M) i_x = 0;
               if (i_y < 0) i_y = M;
               if (i_y > M) i_y = 0;
               if (i_z < 0) i_z = M;
               if (i_z > M) i_z = 0;
               j = head[i_x][i_y][i_z];
               while (j != -1) {        //scan through all the beads in the cell
                  if (i > j) {         //check distance to see if it is within r_cut
                     r_vec = subunit_bead[i].pos - subunit_bead[j].pos;
                     if (r_vec.x > box.x / 2) r_vec.x -= box.x;
                     if (r_vec.x < -box.x / 2) r_vec.x += box.x;
                     if (r_vec.y > box.y / 2) r_vec.y -= box.y;
                     if (r_vec.y < -box.y / 2) r_vec.y += box.y;
                     if (r_vec.z > box.z / 2) r_vec.z -= box.z;
                     if (r_vec.z < -box.z / 2) r_vec.z += box.z;
                     r_dist = r_vec.GetMagnitude();
                     
                     if (r_dist < r_cut) {
                        subunit_bead[i].itsN.push_back(&subunit_bead[j]);
                        subunit_bead[j].itsN.push_back(&subunit_bead[i]);
                     }
                  }
                  j = list[j];
               }
               
               
            }
         }
      }
   } */

    
    

    
    if (world.rank() == 0) {
        sysdata << "Box length is " << box_x * SIsigma / (1e-9) << " nanometers." << endl;
        sysdata << "Screening length is " << screen << " nanometers." << endl;
    }


    if (brownian == false)                                    //for molecular, set up the nose hoover thermostat
    {
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

    dress_up(subunit_edge,
             subunit_face);                                               		// Calculate initial forces



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
                printf("Make sure that number of beads is greater than %d\n",
                       omp_get_num_threads() * numOfNodes);
            }
        }
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									Initial Force Calculation												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

int test = 0;

   for (unsigned int i = 0; i < subunit_bead.size(); i ++) {
      subunit_bead[i].itsN.assign(subunit_bead.size(), 1);                      //To turn pairlist on, comment this line and uncomment the next one. Uncomment lines 50-67 in forces.cpp
      //subunit_bead[i].itsN.assign(subunit_bead.size(), 0); 
   }
  
    forceCalculation(protein, lb, ni, qs, subunit_bead, lj_pairlist, ecut, ks, bondlength, kb, lj_a, ecut_el, kappa, elj_att, test);


    double senergy = 0;                                                //blank all the energy metrics
    double kenergy = 0;
    double benergy = 0;
    double tenergy = 0;
    double ljenergy = 0;
    double cenergy = 0;
    double tpenergy = 0;
    double tkenergy = 0;

    if (restartFile == false) {
    initialize_bead_velocities(protein, subunit_bead, T);        //assign random velocities based on initial temperature
//    initialize_constant_bead_velocities(protein, subunit_bead, T);
    }

    double particle_ke = particle_kinetic_energy(subunit_bead);     //thermostat variables for nose hoover
    double expfac_real;//= exp(-0.5 * delta_t * real_bath[0].xi);

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);                    //setting up random seed for brownian
    unsigned long int Seed = 23410981;
    gsl_rng_set(r, Seed);



/*                  ___                        __      __      ___
      /|   /|      |   \             |        /  \    /  \    |   \
     / |  / |      |    \            |       |    |  |    |   |    |
    /  | /  | ---- |     |           |       |    |  |    |   |___/
   /   |/   |      |    /            |       |    |  |    |   |
  /         |      |___/             |_____   \__/    \__/    |                       */

int loopStart;
if (restartFile == false) loopStart = 0;
if (restartFile == true) loopStart = restartStep;

    for (unsigned int a = loopStart; a < (totaltime / delta_t); a++)         // BEGIN MD LOOP
    {

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								VELOCITY VERLET															*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (brownian == false)    //FOR MOLECULAR DYNAMICS
        {
            for (int i = real_bath.size() - 1; i > -1; i--)                    //thermostat update
                update_chain_xi(i, real_bath, delta_t, particle_ke);
            for (unsigned int i = 0; i < real_bath.size(); i++)
                real_bath[i].update_eta(delta_t);

            expfac_real = exp(-0.5 * delta_t * real_bath[0].xi);

            for (unsigned int i = 0; i < protein.size(); i++)                //velocity verlet loop
            {
                for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
                    protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);  //update velocity half step
                    protein[i].itsB[ii]->update_position(delta_t);                                   //update position full step
                }
            }
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
                }
            }
        }

        dress_up(subunit_edge, subunit_face);                              //update edge and face properties
        
        
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                                                      UPDATE PAIRLIST                                                                                          */
//////////////////////////////////////////////////////////////////////////////////////////////////////////  
// if (a % 20) {
// 
// for (unsigned int ii = 0; ii < subunit_bead.size(); ii++) {
//    subunit_bead[ii].itsN.clear();                                //clear the pairlist
// }
// 
// 
// //       for (unsigned int i = 0; i < protein.size(); i++) {
// //            for (unsigned int j = i + 1; j < protein.size(); j++) {
// //               if (dist(protein[i].itsB[8], protein[j].itsB[8]).GetMagnitudeSquared() < (9.7766*9.7766) ) { // If it is within the cutoff...
// //                  for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
// //                     for (unsigned int jj = 0; jj < protein[j].itsB.size(); jj++) {
// //                        protein[i].itsB[ii]->itsN.push_back(protein[j].itsB[jj]); //Add it to the pairlist
// //                        protein[j].itsB[jj]->itsN.push_back(protein[i].itsB[ii]);
// //                   }
// //                }
// //             }
// //          }
// //       }
// 
//          for (unsigned int ii = 0; ii < subunit_bead.size(); ii++) {
//             for (unsigned int jj = ii + 1; jj < subunit_bead.size(); jj++) {
//                if (dist(&subunit_bead[ii], &subunit_bead[jj]).GetMagnitudeSquared() < (6.5*6.5) ) {
//                   subunit_bead[ii].itsN.push_back(&subunit_bead[jj]);
//                   subunit_bead[jj].itsN.push_back(&subunit_bead[ii]);
//                }
//             }
//          }
//         
//         
//         
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									MD LOOP FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

        forceCalculation(protein, lb, ni, qs, subunit_bead, lj_pairlist, ecut, ks, bondlength, kb, lj_a, ecut_el, kappa, elj_att, a);


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								VELOCITY VERLET															*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////


        if (brownian == false) {        //FOR MOLECULAR DYNAMICS
            for (unsigned int i = 0; i < protein.size(); i++) {
                for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
                    protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);  //update velocity the other half step
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
        }

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

         if (a % 1000 == 0) {
            restart.open("outfiles/restart.out", ofstream::out | ofstream::trunc);
            
            restart << "Velocities & Positions for " << a << endl;
            for (unsigned int i = 0; i < subunit_bead.size(); i++) {
               restart << i << "  " << subunit_bead[i].vel.x << setw(25) << setprecision(12) << subunit_bead[i].vel.y << setw(25) << setprecision(12) << subunit_bead[i].vel.z  << setw(25) << setprecision(12)
                       << subunit_bead[i].pos.x << setw(25) << setprecision(12) << subunit_bead[i].pos.y << setw(25) << setprecision(12) << subunit_bead[i].pos.z  << setw(25) << setprecision(12) << endl;
            }
            restart.close();
         }

           
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								ANALYZE ENERGIES														*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////


            for (unsigned int i = 0; i < protein.size(); i++)        //blanking out energies here
            {
                for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
                    protein[i].itsB[ii]->be = 0;
                    protein[i].itsB[ii]->ne = 0;
                    protein[i].itsB[ii]->ce = 0;
                }
            }


            for (unsigned int i = 0; i < protein.size(); i++)        // Intramolecular Energies
            {
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

            for (unsigned int i = 0; i < protein.size(); i++)        //blanking out energies here
            {
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
                      << "ITEM: BOX BOUNDS" << endl << -bxsz.x / 2 << setw(15) << bxsz.x / 2 << endl << -bxsz.y / 2
                      << setw(15)
                      << bxsz.y / 2 << endl << -bxsz.z / 2 << setw(15) \
 << bxsz.z / 2 << endl << "ITEM: ATOMS index type x y z b charge" << endl;


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
}


