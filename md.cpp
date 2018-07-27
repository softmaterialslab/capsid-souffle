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


void run_molecular () {

    double const bondlength = pow(2, 0.166666666667);       //System specific paramters (modelled for HBV)
    double const SImass = 6.18e-24; //kg               //SI value for a single bead
    double const SIsigma = 1.67e-9; //m
    double const SItime = 7.06e-11; //s
    double const Avagadro = 6.022e23; // mol^-1
    int filenumber = 100000;

    vector<BEAD> subunit_bead;                                      //Create particles, named subunit_bead 
    vector<EDGE> subunit_edge;                                     //create edges between subunit_bead's
    vector<SUBUNIT> protein;                                  //create subunits, named protein
    vector<FACE> subunit_face;                                     //create faces of subunit_bead's on protein
    vector<PAIR> lj_pairlist;                                     //create vector to hold LJ pairings
    vector<OLIGOMER> oligomers_list;                                 //setting up other objects used in analysis
    vector<THERMOSTAT> real_bath;                           //vector of thermostats
    

    double capsomere_concentration, salt_concentration, ks, kb ;
    string file_name; // = "41part_c";
	double totaltime ;//= 100;				//total time in MD units
	double ecut = 2.5;					//lennard jones cut-off distance
	double qs = 1;						//salt valency
	double delta_t ;//= .002;					//time step in MD units
	double number_capsomeres = 8;					//number of subunits in the box
	double chain_length_real = 5;		//nose hoover chain length
	double Q = 1;						//nose hoover mass (reduced units)
	double T = 1;						//set temperature (reduced units)
   
    cout << "Filename?" << endl;								 cin >> file_name;
    cout << "capsomere concentration (micromolar):" << endl;     cin >> capsomere_concentration;
    cout << "salt concentration (millimolar):" << endl;          cin >> salt_concentration;
    cout << "stretching constant (KbT):" << endl;          		 cin >> ks;
    cout << "bending constant (KbT):" << endl;             		 cin >> kb;
	cout << "total time (MD steps):" << endl;					 cin >> totaltime;
	cout << "timestep (MD steps):" << endl;						 cin >> delta_t;
        



    generate_lattice(capsomere_concentration, number_capsomeres, file_name);     //Setting up the input file (uses user specified file to generate lattice)

    double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * 6.022e23)), 1.0 / 3.0);    //calculating box size
    VECTOR3D bxsz = VECTOR3D(box_x, box_x, box_x);

    initialize_system(subunit_bead, subunit_edge, protein, subunit_face, bxsz, lj_pairlist);


                                                                        //user-derived parameters (not edittable)
    double lb = 0.416;                                // e^2 / (4 pi Er E0 Kb T)
    double ni = salt_concentration * Avagadro * SIsigma * SIsigma * SIsigma;       //number density (1/sigma*^3)
    int count = 0;                                    //used in oligomer counting later on...
    vector<int> massbins(protein.size());
    vector<vector<int> > ms_bin(totaltime / (delta_t * 1000), vector<int>(protein.size()));
	
	

    cout << endl << "Simulation will run for " << totaltime * SItime / (1e-9) << " nanoseconds with a "
         << delta_t * SItime / (1e-12) << " picosecond timestep." << endl;

    ofstream traj("energy.out", ios::out);              //setting up file outputs
    ofstream ofile("ovito.lammpstrj", ios::out);
    ofstream msdata("ms.out", ios::out);
    initialize_outputfile(traj, ofile);

    /*                                                    NOSE-HOOVER CHAIN THERMOSTAT                                  */
    if (chain_length_real == 1)
        real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0, 1)));
    else
    {
        real_bath.push_back((THERMOSTAT(Q, T, 3 * subunit_bead.size(), 0, 0, 0, 1)));
        while (real_bath.size() != chain_length_real - 1)
            real_bath.push_back((THERMOSTAT(Q / (3 * subunit_bead.size()), T, 1, 0, 0, 0,1)));
        real_bath.push_back((THERMOSTAT(0, T, 3 * subunit_bead.size(), 0.0, 0, 0,1)));
// final bath is dummy bath (dummy bath always has zero mass)
    }



/*                  ___                        __      __      ___
      /|   /|      |   \             |        /  \    /  \    |   \
     / |  / |      |    \            |       |    |  |    |   |    |
    /  | /  | ---- |     |           |       |    |  |    |   |___/
   /   |/   |      |    /            |       |    |  |    |   |
  /         |      |___/             |_____   \__/    \__/    |                       */

// Calculate initial forces


    dress_up(subunit_edge, subunit_face);
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTRA MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    for (unsigned int i = 0; i < protein.size(); i++) 
	{
		for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
		{
			protein[i].itsB[ii]->update_stretching_force(ks, bondlength);
			protein[i].itsB[ii]->bforce = VECTOR3D(0,0,0);		//resetting bending force here
		}
		
		for (unsigned int m = 0; m < protein[i].itsE.size(); m++)
		{
			if (protein[i].itsE[m]->type != 0)					//if it is a bending edge...
				protein[i].itsE[m]->update_bending_forces(kb);
		}
	}


 //////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTER MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////   

    update_ES_forces(protein, lb, ni, qs);		//ALSO INTRA-MOLECULAR
	
    update_LJ_forces(subunit_bead, ecut, lj_pairlist);

    double senergy = 0;                            //blank all the energy metrics
    double kenergy = 0;
    double benergy = 0;
    double tenergy = 0;
    double ljenergy = 0;
    double cenergy = 0;
    double tpenergy = 0;
    double tkenergy = 0;

    //initialize_bead_velocities(protein, subunit_bead, T);            //assign random velocities based on initial temperature
	initialize_constant_bead_velocities(protein,subunit_bead,T);

    double particle_ke = particle_kinetic_energy(subunit_bead);     //thermostat variables
    double expfac_real;

    int mstime = -1;                                          //parameter for ms_bin filling

    cout << endl << endl << "Calculating..." << endl << endl;
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								BEGIN MD LOOP															*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    for (unsigned int a = 0; a < (totaltime / delta_t); a++)         // BEGIN MD LOOP
    {
                                                       
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								THERMOSTAT UPDATE														*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////																
        for (int i = real_bath.size() - 1; i > -1; i--)
            update_chain_xi(i, real_bath, delta_t, particle_ke);
        for (unsigned int i = 0; i < real_bath.size(); i++)
            real_bath[i].update_eta(delta_t);

        expfac_real = exp(-0.5 * delta_t * real_bath[0].xi);
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								VELOCITY VERLET															*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

		for (unsigned int i = 0; i < protein.size(); i++)
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
					{
						protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);  //update velocity half step
						protein[i].itsB[ii]->update_position(delta_t);									  //update position full step
					}
			}

        dress_up(subunit_edge, subunit_face);                              //update edge and face properties
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTRA MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

for (unsigned int i = 0; i < protein.size(); i++) 
{
	for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
	{
		protein[i].itsB[ii]->update_stretching_force(ks, bondlength);
		protein[i].itsB[ii]->bforce = VECTOR3D(0,0,0);			//zeroing bending force here
	}
	
	for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++)
	{
		if (protein[i].itsE[kk]->type != 0)						//if it is a bending edge...
			protein[i].itsE[kk]->update_bending_forces(kb);
	}
}
        
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTER MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        update_ES_forces(protein, lb, ni, qs);		//ALSO INTRAMOLECULAR

        update_LJ_forces(subunit_bead, ecut, lj_pairlist);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								VELOCITY VERLET															*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////        

        for (unsigned int i = 0; i < protein.size(); i++)
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
					{
						protein[i].itsB[ii]->therm_update_velocity(delta_t, real_bath[0], expfac_real);  //update velocity the other half step
					}
			}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								THERMOSTAT UPDATE														*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        particle_ke = particle_kinetic_energy(subunit_bead);
																				// forward update of Nose-Hoover chain
        for (unsigned int i = 0; i < real_bath.size(); i++)
            real_bath[i].update_eta(delta_t);
        for (unsigned int i = 0; i < real_bath.size(); i++)
            update_chain_xi(i, real_bath, delta_t, particle_ke);


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

         
            for (unsigned int i = 0; i < protein.size(); i++) 		//blanking out energies here 
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
				{
					protein[i].itsB[ii]->be = 0;
					protein[i].itsB[ii]->ne = 0;
					protein[i].itsB[ii]->ce = 0;
				}
			}
			
            
            for (unsigned int i = 0; i < protein.size(); i++) 		// Intramolecular Energies
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
				{
					protein[i].itsB[ii]->update_stretching_energy(ks, bondlength);
					protein[i].itsB[ii]->update_kinetic_energy();
				}
				for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++)
				{
					if (protein[i].itsE[kk]->type != 0)			//if it is a bending edge...
						protein[i].itsE[kk]->update_bending_energy(kb);
				}
			}
																	//Intermolecular Energies
            update_ES_energies(protein, lb, ni, qs);

            update_LJ_energies(subunit_bead, ecut, lj_pairlist);

           for (unsigned int i = 0; i < protein.size(); i++) 		//blanking out energies here 
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
				{
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
                tpenergy += real_bath[i].pe;                    	//sum up total energies
                tkenergy += real_bath[i].ke;
            }

            tenergy = senergy + kenergy + ljenergy + benergy + cenergy + tpenergy +
                      tkenergy;      //print info to files for data analysis
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								STORE ENERGY INFO TO FILE												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

            ofile << "ITEM: TIMESTEP" << endl << a << endl << "ITEM: NUMBER OF ATOMS" << endl << subunit_bead.size() << endl
                  << "ITEM: BOX BOUNDS" << endl << -bxsz.x / 2 << setw(15) << bxsz.x / 2 << endl << -bxsz.y / 2
                  << setw(15)
                  << bxsz.y / 2 << endl << -bxsz.z / 2 << setw(15) \
				  << bxsz.z / 2 << endl << "ITEM: ATOMS index type x y z b charge" << endl;

            traj << a * delta_t << setw(15) << kenergy / subunit_bead.size() << setw(15) << senergy / subunit_bead.size() << setw(15) <<
                 benergy / subunit_bead.size() << setw(15) << ljenergy / subunit_bead.size() << setw(15) << cenergy / subunit_bead.size()
                 << setw(15) << tenergy / subunit_bead.size()
                 << setw(15) << (benergy + senergy + ljenergy + cenergy) / subunit_bead.size() << setw(15)
                 << kenergy * 2 / (3 * subunit_bead.size()) << setw(15) << tpenergy / subunit_bead.size() << setw(15)
                 << tkenergy / subunit_bead.size() << endl;

            for (unsigned int b = 0; b < subunit_bead.size(); b++) {
                ofile << b + 1 << setw(15) << subunit_bead[b].type << setw(15) << subunit_bead[b].pos.x << setw(15) << subunit_bead[b].pos.y \
					  << setw(15) << subunit_bead[b].pos.z << setw(15) << subunit_bead[b].be << setw(15) << subunit_bead[b].q << endl;

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
                            for (int k = i + 1; k < protein.size(); k++) { //look for new growth
                                if (protein[k].itsO.size() == 0) {  //if it isn't in an oligomer yet...
                                    for (int m = 0;
                                         m < protein[g].itsB.size(); m++) { //check to see if it is in this oligomer
                                        for (int n = 0; n < protein[k].itsB.size(); n++) {
                                            if (dist(protein[g].itsB[m], protein[k].itsB[n]).GetMagnitude() < 1.5) {
                                                oligomers_list[index].itsS.push_back(&protein[k]);   //if it is attached, add it
                                                protein[k].itsO.push_back(oligomers_list[index]);   //mark subunit as bonded
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
            for (int i = 0; i < oligomers_list.size(); i++) {
                if (oligomers_list[i].itsS.size() >= 1) {
                    ms_bin[mstime][(oligomers_list[i].itsS.size() - 1)] += 1;          //fill mass bins
                }
            }

            for (int j = 0; j < ms_bin[mstime].size(); j++) {
                msdata << ms_bin[mstime][j] << setw(15);                //print mass bin data to file
            }
            msdata << endl;

            for (int i = 0; i < protein.size(); i++) {             // clear oligomer pointers from subunit
                protein[i].itsO.clear();
            }

            oligomers_list.erase(oligomers_list.begin(),oligomers_list.end());                 //erases oligomer objects


        }//end of energy analysis loop


        double fraction_completed = ((a + 1) / (totaltime / delta_t));   //progress bar
        ProgressBar(fraction_completed);
    } //time loop end



    compute_MD_trust_factor_R(1);                   //computes R


}


























void run_brownian(){

    double const bondlength = pow(2, 0.166666666667);       //System specific paramters (modelled for HBV)
    double const SImass = 6.18e-24; //kg               //SI value for a single bead
    double const SIsigma = 1.67e-9; //m
    double const SItime = 7.06e-11; //s
    double const Avagadro = 6.022e23; // mol^-1
    int filenumber = 100000;

    vector<BEAD> subunit_bead;                                      //Create particles, named subunit_bead (short for protein)
    vector<EDGE> subunit_edge;                                     //create edges between subunit_bead's
    vector<SUBUNIT> protein;                                  //create subunits, named protein
    vector<FACE> subunit_face;                                     //create faces of subunit_bead's on protein
    vector<PAIR> lj_pairlist;                                     //create vector to hold LJ pairings
    

    double capsomere_concentration, salt_concentration, ks, kb , fric_zeta;
    string file_name = "41part";
	double totaltime = 100;				//total time in MD units
	double ecut = 2.5;					//lennard jones cut-off distance
	double qs = 1;						//salt valency
	double delta_t = .002;					//time step in MD units
	double number_capsomeres = 8;					//number of subunits in the box
	double chain_length_real = 5;		//nose hoover chain length
	double Q = 1;						//nose hoover mass (reduced units)
	double T = 1;						//set temperature (reduced units)
   
    cout << "capsomere concentration (micromolar):" << endl;     cin >> capsomere_concentration;
    cout << "salt concentration (millimolar):" << endl;          cin >> salt_concentration;
    cout << "stretching constant (KbT):" << endl;          		 cin >> ks;
    cout << "bending constant (KbT):" << endl;             		 cin >> kb;
    cout << "friction coefficient (reduced unit):" << endl;      cin >> fric_zeta;



    generate_lattice(capsomere_concentration, number_capsomeres, file_name);     //Setting up the input file (uses user specified file to generate lattice)

    double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * 6.022e23)), 1.0 / 3.0);    //calculating box size
    VECTOR3D bxsz = VECTOR3D(box_x, box_x, box_x);

    initialize_system(subunit_bead, subunit_edge, protein, subunit_face, bxsz, lj_pairlist);


    //user-derived parameters (not edittable)
    double lb = 0.416;                                // e^2 / (4 pi Er E0 Kb T)
    double ni = salt_concentration * Avagadro * SIsigma * SIsigma * SIsigma;       //number density (1/sigma*^3)
    int count = 0;                                    //used in oligomer counting later on...


    vector<OLIGOMER> oligomers_list;                                 //setting up other objects used in analysis
    vector<int> massbins(protein.size());
    vector<vector<int> > ms_bin(totaltime / (delta_t * 1000), vector<int>(protein.size()));

    cout << endl << "Simulation will run for " << totaltime * SItime / (1e-9) << " nanoseconds with a "
         << delta_t * SItime / (1e-12) << " picosecond timestep." << endl;

    ofstream traj("energy.out", ios::out);              //setting up file outputs
    ofstream ofile("ovito.lammpstrj", ios::out);
    ofstream msdata("ms.out", ios::out);
    initialize_outputfile(traj, ofile);



/*                  ___                        __      __      ___
      /|   /|      |   \             |        /  \    /  \    |   \
     / |  / |      |    \            |       |    |  |    |   |    |
    /  | /  | ---- |     |           |       |    |  |    |   |___/
   /   |/   |      |    /            |       |    |  |    |   |
  /         |      |___/             |_____   \__/    \__/    |                       */

// Calculate initial forces
    dress_up(subunit_edge, subunit_face);
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTRA MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
     for (unsigned int i = 0; i < protein.size(); i++) 
	{
		for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
		{
			protein[i].itsB[ii]->update_stretching_force(ks, bondlength);
			protein[i].itsB[ii]->bforce = VECTOR3D(0,0,0);		//resetting bending force here
		}
		
		for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++)
		{
			if (protein[i].itsE[kk]->type != 0)					//if it is a bending edge...
				protein[i].itsE[kk]->update_bending_forces(kb);
		}
	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTER MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

    update_ES_forces(protein, lb, ni, qs);

    update_LJ_forces(subunit_bead, ecut, lj_pairlist);

    double senergy = 0;                                 //blank all the energy metrics
    double kenergy = 0;
    double benergy = 0;
    double tenergy = 0;
    double ljenergy = 0;
    double cenergy = 0;

    initialize_bead_velocities(protein, subunit_bead, T);            //assign random velocities based on initial temperature

    int mstime = -1;                                          //parameter for ms_bin filling

    cout << endl << endl << "Calculating..." << endl << endl;

    gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long int Seed = 23410981;
    gsl_rng_set(r,Seed);

    for (unsigned int a = 0; a < (totaltime / delta_t); a++)         // BEGIN MD LOOP
    {


        for (unsigned int i = 0; i < subunit_bead.size(); i++) {
            //subunit_bead[i].brownian_update_velocity(delta_t, fric_zeta);  //update velocity half step
            subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + subunit_bead[i].ljforce + subunit_bead[i].eforce;
            subunit_bead[i].vel.x +=
                    (subunit_bead[i].vel.x * (-0.5 * fric_zeta * delta_t)) + (subunit_bead[i].tforce.x * (0.5 * delta_t / subunit_bead[i].m)) +
                    sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) * (gsl_rng_uniform(r) - 0.5);
            subunit_bead[i].vel.y +=
                    (subunit_bead[i].vel.y * (-0.5 * fric_zeta * delta_t)) + (subunit_bead[i].tforce.y * (0.5 * delta_t / subunit_bead[i].m)) +
                    sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) * (gsl_rng_uniform(r) - 0.5);
            subunit_bead[i].vel.z +=
                    (subunit_bead[i].vel.z * (-0.5 * fric_zeta * delta_t)) + (subunit_bead[i].tforce.z * (0.5 * delta_t / subunit_bead[i].m)) +
                    sqrt(2 * 6 * delta_t * fric_zeta / subunit_bead[i].m) * (gsl_rng_uniform(r) - 0.5);
        }
			
		for (unsigned int i = 0; i < protein.size(); i++)
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
					{
						protein[i].itsB[ii]->update_position(delta_t);  //update position full step
					}
			}

        dress_up(subunit_edge, subunit_face);                              //update edge and face properties
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTRA MOLECULAR FORCES												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (unsigned int i = 0; i < protein.size(); i++) 
	{
		for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
		{
			protein[i].itsB[ii]->update_stretching_force(ks, bondlength);
			protein[i].itsB[ii]->bforce = VECTOR3D(0,0,0);		//resetting bending force here
		}
		
		for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++)
		{
			if (protein[i].itsE[kk]->type != 0)					//if it is a bending edge...
				protein[i].itsE[kk]->update_bending_forces(kb);
		}
	}
        
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*									INTER MOLECULAR FORCES												*/
////////////////////////////////////////////////////////////////////////////////////////////////////////// 

        update_ES_forces(protein, lb, ni, qs);

        update_LJ_forces(subunit_bead, ecut, lj_pairlist);



        for (unsigned int i = 0; i < subunit_bead.size(); i++) {
            //subunit_bead[i].brownian_update_velocity(delta_t, fric_zeta);    				//update velocity the other half step
            subunit_bead[i].tforce = subunit_bead[i].sforce + subunit_bead[i].bforce + subunit_bead[i].ljforce + subunit_bead[i].eforce;
            subunit_bead[i].vel.x += (subunit_bead[i].vel.x*(-0.5*fric_zeta*delta_t)) + (subunit_bead[i].tforce.x*(0.5*delta_t/subunit_bead[i].m)) +
                             sqrt(2*6*delta_t*fric_zeta/subunit_bead[i].m)*(gsl_rng_uniform(r)-0.5);
            subunit_bead[i].vel.y += (subunit_bead[i].vel.y*(-0.5*fric_zeta*delta_t)) + (subunit_bead[i].tforce.y*(0.5*delta_t/subunit_bead[i].m)) +
                             sqrt(2*6*delta_t*fric_zeta/subunit_bead[i].m)*(gsl_rng_uniform(r)-0.5);
            subunit_bead[i].vel.z += (subunit_bead[i].vel.z*(-0.5*fric_zeta*delta_t)) + (subunit_bead[i].tforce.z*(0.5*delta_t/subunit_bead[i].m)) +
                             sqrt(2*6*delta_t*fric_zeta/subunit_bead[i].m)*(gsl_rng_uniform(r)-0.5);
        }




/*      __                 __                        ____     ________     ____
 *     /  \    |\    |    /  \    |       \   /     /    \        |       /    \
 *    |____|   | \   |   |____|   |        \ /      \_            |       \_
 *    |    |   |  \  |   |    |   |         |         \__         |         \__
 *    |    |   |   \ |   |    |   |         |            \        |            \
 *    |    |   |    \|   |    |   |_____    |       \____/    ____|____   \____/                                */



//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								ENERGY ANALYSIS															*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (a % 100 == 0) {                                           //analysis loop (energies)
			 for (unsigned int i = 0; i < protein.size(); i++) 		//blanking out energies here 
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
				{
					protein[i].itsB[ii]->be = 0;
					protein[i].itsB[ii]->ne = 0;
					protein[i].itsB[ii]->ce = 0;
				}
			}
			
            
            for (unsigned int i = 0; i < protein.size(); i++) 		// Intramolecular Energies
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
				{
					protein[i].itsB[ii]->update_stretching_energy(ks, bondlength);
					protein[i].itsB[ii]->update_kinetic_energy();
				}
				for (unsigned int kk = 0; kk < protein[i].itsE.size(); kk++)
				{
					if (protein[i].itsE[kk]->type != 0)			//if it is a bending edge...
						protein[i].itsE[kk]->update_bending_energy(kb);
				}
			}
																	//Intermolecular Energies
            update_ES_energies(protein, lb, ni, qs);

            update_LJ_energies(subunit_bead, ecut, lj_pairlist);

           for (unsigned int i = 0; i < protein.size(); i++) 		//blanking out energies here 
			{
				for (unsigned int ii = 0; ii < protein[i].itsB.size(); ii++)
				{
					senergy += protein[i].itsB[ii]->se;                          //sum up total energies
					kenergy += protein[i].itsB[ii]->ke;
					ljenergy += protein[i].itsB[ii]->ne;
					benergy += protein[i].itsB[ii]->be;
					cenergy += protein[i].itsB[ii]->ce;
				}
			}


            tenergy = senergy + kenergy + ljenergy + benergy + cenergy ; //print info to files for data analysis

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*								PRINT ENERGY INFO TO FILE												*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

            ofile << "ITEM: TIMESTEP" << endl << a << endl << "ITEM: NUMBER OF ATOMS" << endl << subunit_bead.size() << endl
                  << "ITEM: BOX BOUNDS" << endl << -bxsz.x / 2 << setw(15) << bxsz.x / 2 << endl << -bxsz.y / 2
                  << setw(15)
                  << bxsz.y / 2 << endl << -bxsz.z / 2 << setw(15) \
				  << bxsz.z / 2 << endl << "ITEM: ATOMS index type x y z b charge" << endl;

            traj << a * delta_t << setw(15) << kenergy / subunit_bead.size() << setw(15) << senergy / subunit_bead.size() << setw(15) <<
                 benergy / subunit_bead.size() << setw(15) << ljenergy / subunit_bead.size() << setw(15) << cenergy / subunit_bead.size()
                 << setw(15) << tenergy / subunit_bead.size()
                 << setw(15) << (benergy + senergy + ljenergy + cenergy) / subunit_bead.size() << setw(15)
                 << kenergy * 2 / (3 * subunit_bead.size()) << setw(15) << endl;

            for (unsigned int b = 0; b < subunit_bead.size(); b++) {
                ofile << b + 1 << setw(15) << subunit_bead[b].type << setw(15) << subunit_bead[b].pos.x << setw(15) << subunit_bead[b].pos.y \
					  << setw(15) << subunit_bead[b].pos.z << setw(15) << subunit_bead[b].be << setw(15) << subunit_bead[b].q << endl;

                count += 1;

            }
            senergy = 0;                            //blanking out energies
            kenergy = 0;
            benergy = 0;
            tenergy = 0;
            ljenergy = 0;
            cenergy = 0;

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
                            for (int k = i + 1; k < protein.size(); k++) { //look for new growth
                                if (protein[k].itsO.size() == 0) {  //if it isn't in an oligomer yet...
                                    for (int m = 0;
                                         m < protein[g].itsB.size(); m++) { //check to see if it is in this oligomer
                                        for (int n = 0; n < protein[k].itsB.size(); n++) {
                                            if (dist(protein[g].itsB[m], protein[k].itsB[n]).GetMagnitude() < 1.5) {
                                                oligomers_list[index].itsS.push_back(&protein[k]);   //if it is attached, add it
                                                protein[k].itsO.push_back(oligomers_list[index]);   //mark subunit as bonded
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
            for (int i = 0; i < oligomers_list.size(); i++) {
                if (oligomers_list[i].itsS.size() >= 1) {
                    ms_bin[mstime][(oligomers_list[i].itsS.size() - 1)] += 1;          //fill mass bins
                }
            }

            for (int j = 0; j < ms_bin[mstime].size(); j++) {
                msdata << ms_bin[mstime][j] << setw(15);                //print mass bin data to file
            }
            msdata << endl;

            for (int i = 0; i < protein.size(); i++) {             // clear oligomer pointers from subunit
                protein[i].itsO.clear();
            }

            oligomers_list.erase(oligomers_list.begin(),oligomers_list.end());                 //erases oligomer objects
// g(r) dump data


        }//end of energy analysis loop


        double fraction_completed = ((a + 1) / (totaltime / delta_t));   //progress bar
        ProgressBar(fraction_completed);
    } //time loop end



    compute_MD_trust_factor_R(1);                   //computes R


}
