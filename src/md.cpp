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
    double capsomere_concentration, salt_concentration, ks, kb;
    double totaltime;                                        //total time in MD units
    double delta_t;;                                //time step in MD units
    double fric_zeta;
    bool verbose;

    // Get input values from the user
    options_description desc("Usage:\nrandom_mesh <options>");
    desc.add_options()
            ("help,h", "print usage message")
            ("Dynamics Response,D", value<char>(&response)->default_value('m'),
             "To run brownian dynamics (overdamped langevin) enter 'b'. Otherwise, to run molecular dynamics with nose' hoover thermostat enter 'm'. [b/m]")
            ("filename,f", value<string>(&file_name)->default_value("41part"), "Filename?")
            ("capsomere concentration,C", value<double>(&capsomere_concentration)->default_value(931),
             "capsomere concentration (micromolar)")
            ("salt concentration,c", value<double>(&salt_concentration)->default_value(0),
             "salt concentration (millimolar)")
            ("stretching constant,s", value<double>(&ks)->default_value(100), "stretching constant (KbT)")
            ("bending constant,b", value<double>(&kb)->default_value(20), "bending constant (KbT)")
            ("total time,T", value<double>(&totaltime)->default_value(100), "total time (MD steps)")
            ("timestep,t", value<double>(&delta_t)->default_value(0.002), "timestep (MD steps)")
            ("friction coefficient,r", value<double>(&fric_zeta)->default_value(1),
             "friction coefficient (reduced unit)")
            ("verbose,V", value<bool>(&verbose)->default_value(true), "verbose true: provides detailed output");


    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }


    ofstream traj("outfiles/energy.out", ios::out);              //setting up file outputs
    ofstream ofile("outfiles/ovito.lammpstrj", ios::out);
    ofstream msdata("outfiles/ms.out", ios::out);
    ofstream sysdata("outfiles/model.parameters.out", ios::out);

    initialize_outputfile(traj, ofile);


    if (response == 'b') {                    //set flag for brownian vs. molecular dynamics
        brownian = true;
        if (world.rank() == 0)
            sysdata << "Running brownian dynamics." << endl;
    } else if (response == 'm') {
        brownian = false;
        if (world.rank() == 0)
            sysdata << "Running molecular dynamics with Nose Hoover thermostat." << endl;
    }

    double bondlength;//= pow(2, 0.166666666667);       //System specific paramters (modelled for HBV)
    double SImass;//= 6.18e-24; //kg               		//SI value for a single bead
    double SIsigma;//= 1.67e-9; //m
    double SItime;//= 7.06e-11; //s
    double const Avagadro = 6.022e23; // mol^-1				//useful constants
    //int filenumber = 100000;                                //used in pair correlation file generation

    vector<BEAD> subunit_bead;                              //Create particles, named subunit_bead
    vector<EDGE> subunit_edge;                              //create edges between subunit_bead's
    vector<SUBUNIT> protein;                                //create subunits, named protein
    vector<FACE> subunit_face;                              //create faces of subunit_bead's on protein
    vector<PAIR> lj_pairlist;                               //create vector to hold LJ pairings
    vector<OLIGOMER> oligomers_list;                        //create vector to hold oligomers for mass spectrum analysis
    vector<THERMOSTAT> real_bath;                           //vector of thermostats

    double ecut = 2.5;                                        //lennard jones cut-off distance
    double qs = 1;                                            //salt valency
    double number_capsomeres = 8;                            //number of subunits in the box
    double T = 1;                                            //set temperature (reduced units)
    double chain_length_real = 5;                           //nose hoover chain length
    double Q = 1;                                           //nose hoover mass (reduced units)


    vector<vector<int> > lj_a;
    lj_a = generate_lattice(capsomere_concentration, number_capsomeres, file_name, bondlength, SIsigma, SImass,
                            SItime, subunit_bead, subunit_edge, protein,
                            subunit_face);     //Setting up the input file (uses user specified file to generate lattice)

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

    }
    double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * Avagadro)),
                       1.0 / 3.0);    //calculating box size
    VECTOR3D bxsz = VECTOR3D(box_x, box_x, box_x);

    if (world.rank() == 0) {
        sysdata << "Box length is " << box_x * SIsigma / (1e-9) << " nanometers." << endl;
        sysdata << "Screening length is " << 8.7785 / sqrt(salt_concentration) << " nanometers." << endl;
    }


    //user-derived parameters (not edittable)
    double lb = (0.76e-9) / SIsigma;                                        // e^2 / (4 pi Er E0 Kb T)
    double ni = salt_concentration * Avagadro * SIsigma * SIsigma * SIsigma;       //number density (1/sigma*^3)
    int count = 0;                                            //used in mass spectrum analysis
    int mstime = -1;                                        //parameter for ms_bin filling	
    vector<int> massbins(protein.size());
    vector<vector<int> > ms_bin(totaltime / (delta_t * 1000), vector<int>(protein.size()));


    if (brownian == false)                                    //for molecular, set up the nose hoover thermostat
    {
        if (chain_length_real == 1)
            real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0, 1)));
        else {
            real_bath.push_back((THERMOSTAT(Q, T, 3 * subunit_bead.size(), 0, 0, 0, 1)));
            while (real_bath.size() != chain_length_real - 1)
                real_bath.push_back((THERMOSTAT(Q / (3 * subunit_bead.size()), T, 1, 0, 0, 0, 1)));
            real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0, 1)));
            // final bath is dummy bath (dummy bath always has zero mass)
        }
    }

    dress_up(subunit_edge,
             subunit_face);                                                                    // Calculate initial forces



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


    forceCalculation(protein, lb, ni, qs, subunit_bead, lj_pairlist, ecut, ks, bondlength, kb, lj_a);


    double senergy = 0;                                                //blank all the energy metrics
    double kenergy = 0;
    double benergy = 0;
    double tenergy = 0;
    double ljenergy = 0;
    double cenergy = 0;
    double tpenergy = 0;
    double tkenergy = 0;

    //initialize_bead_velocities(protein, subunit_bead, T);        //assign random velocities based on initial temperature
    initialize_constant_bead_velocities(protein, subunit_bead, T);

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

    for (unsigned int a = 0; a < (totaltime / delta_t); a++)         // BEGIN MD LOOP
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
                    protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0],
                                                               expfac_real);  //update velocity half step
                    protein[i].itsB[ii]->update_position(
                            delta_t);                                      //update position full step
                }
            }
        } else {            // FOR BROWNIAN DYNAMICS
            for (unsigned int i = 0; i < subunit_bead.size(); i++) {
                //subunit_bead[i].brownian_update_velocity(delta_t, fric_zeta);  //update velocity half step

                //subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + subunit_bead[i].ljforce +
                //                         subunit_bead[i].eforce;

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
/*									INTRA MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

        for (unsigned int i = 0; i < protein.size(); i++) {
            for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
                protein[i].itsB[ii]->update_stretching_force(ks, bondlength);
                protein[i].itsB[ii]->bforce = VECTOR3D(0, 0, 0);            //zeroing bending force here
            }

            for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++) {
                if (protein[i].itsE[kk]->type != 0)                        //if it is a bending edge...
                    protein[i].itsE[kk]->update_bending_forces(kb);
            }
        }


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									MD LOOP FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

        forceCalculation(protein, lb, ni, qs, subunit_bead, lj_pairlist, ecut, ks, bondlength, kb, lj_a);


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								VELOCITY VERLET															*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////


        if (brownian == false) {        //FOR MOLECULAR DYNAMICS
            for (unsigned int i = 0; i < protein.size(); i++) {
                for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++) {
                    protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0],
                                                               expfac_real);  //update velocity the other half step
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
                //subunit_bead[i].brownian_update_velocity(delta_t, fric_zeta);    				//update velocity the other half step

                //subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + subunit_bead[i].ljforce +
                //                         subunit_bead[i].eforce;

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

        if (a % 100 == 0) {                                           //analysis loop

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
            //update_ES_energies(protein, lb, ni, qs);
            update_ES_energies_simplified(subunit_bead, lb, ni, qs);

            //update_LJ_energies(protein, ecut);
            update_LJ_energies_simplified(subunit_bead, ecut, lj_a);

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

                count += 1;

            }
            senergy = 0;                            //blanking out energies
            kenergy = 0;
            benergy = 0;
            tenergy = 0;
            ljenergy = 0;
            cenergy = 0;
            tpenergy = 0;
            tkenergy = 0;
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								RADIAL DISTRIBUTION FUNCTION											*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (a % 1000 == 0) {
//            if (a > 100000) {
//                char filename[100];                                     //pair-correlation fxn file generation
//                filenumber += 1000;
//                sprintf(filename, "data.coords.all.%d", filenumber);
//                ofstream pairout(filename, ios::out);
//                for (int i = 0; i < subunit_bead.size(); i++) {
//                    if (subunit_bead[i].type == 3) {
//                        pairout << subunit_bead[i].id << setw(15) << subunit_bead[i].type << setw(15) << subunit_bead[i].m << setw(15)
//                                << subunit_bead[i].pos.x
//                                << setw(15) << subunit_bead[i].pos.y << setw(15) << subunit_bead[i].pos.z << endl;
//                    }
//                }
//            }

            int index = -1;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								MASS SPECTRUM ANALYSIS													*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
            for (unsigned int i = 0; i < protein.size(); i++)           //Create oligomers for mass spectrum analysis
            {
                double oldsize = 0;

                if (protein[i].itsO.size() == 0) {                       //if the unit isn't already counted...

                    oligomers_list.push_back(OLIGOMER(VECTOR3D(0, 0, 0)));     //create an oligomer for the unit
                    index += 1;
                    oligomers_list[index].itsS.push_back(&protein[i]);      //add unit to oligomer
                    oligomers_list[index].id = index;
                    protein[i].itsO.push_back(oligomers_list[index]);               //add oligomer to unit
                    while (oldsize < oligomers_list[index].itsS.size()) {    //while the oligomer is still growing...
                        int n = oldsize;
                        oldsize = oligomers_list[index].itsS.size();           //see how much the oligomer has grown
                        for (int j = n; j < oldsize; j++) {             //loop over the growth from last round
                            int g = oligomers_list[index].itsS[j]->id;
                            for (unsigned int k = i + 1; k < protein.size(); k++) { //look for new growth
                                if (protein[k].itsO.size() == 0) {  //if it isn't in an oligomer yet...
                                    for (unsigned int m = 0;
                                         m < protein[g].itsB.size(); m++) { //check to see if it is in this oligomer
                                        for (unsigned int n = 0; n < protein[k].itsB.size(); n++) {
                                            if (dist(protein[g].itsB[m], protein[k].itsB[n]).GetMagnitude() < 1.5) {
                                                oligomers_list[index].itsS.push_back(
                                                        &protein[k]);   //if it is attached, add it
                                                protein[k].itsO.push_back(
                                                        oligomers_list[index]);   //mark subunit as bonded
                                                goto finish;
                                            }
                                        }
                                    }


                                }
                                finish:;
                            }
                        }

                    }
                }
            }
            mstime += 1;
            for (unsigned int i = 0; i < oligomers_list.size(); i++) {
                if (oligomers_list[i].itsS.size() >= 1) {
                    ms_bin[mstime][(oligomers_list[i].itsS.size() - 1)] += 1;          //fill mass bins
                }
            }
            if (world.rank() == 0){
                for (unsigned int j = 0; j < ms_bin[mstime].size(); j++) {
                    msdata << ms_bin[mstime][j] << setw(15);                //print mass bin data to file
                }
                msdata << endl;
            }
            for (unsigned int i = 0; i < protein.size(); i++) {             // clear oligomer pointers from subunit
                protein[i].itsO.clear();
            }

            oligomers_list.erase(oligomers_list.begin(),
                                 oligomers_list.end());                 //erases oligomer objects*/


        }//end of energy analysis loop


        double fraction_completed = ((a + 1) / (totaltime / delta_t));   //progress bar
        ProgressBar(fraction_completed);
    } //time loop end



    compute_MD_trust_factor_R(1);                   //computes R


    return 0;
}


