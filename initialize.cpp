//
// Created by lauren on 1/25/18.
//

#include <cstdlib>
#include "initialize.h"
#include "rand_gen.h"

using namespace std;

void initialize_system(vector<BEAD> &subunit_bead, vector<EDGE> &subunit_edge, vector<SUBUNIT> &protein, vector<FACE> &subunit_face, \
                        VECTOR3D bxsz, vector<PAIR> & lj_pairlist)
{
    /*                                            OPEN FILE                                                            */

    ifstream crds;                                      //open coordinates file
    crds.open("input.GEN.out");
    if (!crds) {                                        //check to make sure file is there
        cerr << "ERR: FILE NOT OPENED. Check directory and/or filename.";
        exit(1);
    }
    int charge, index, g1, g2, g3, norm, type;                   //count, charge, index, edge caps
    unsigned long count;
    long double x, y, z, length,radius, mass, epsilon, sigma;                                //x y z coordinates
    string dummy;                                       //dummy string not imported








/*                                                 PARTICLES                                                        */

    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " beads in the system." << endl;
    for(int i=0;i<count;i++)
    {                                                   //Read file and initialize positions, id and charge
        crds >> index >> x >> y >> z >> dummy >> charge >> type >> radius >> mass;
        subunit_bead.push_back(BEAD(VECTOR3D(x,y,z)));
        subunit_bead[index].id = index;
        subunit_bead[index].q = charge;
        subunit_bead[index].type = type;
        subunit_bead[index].sigma = radius;
        subunit_bead[index].m = mass;                           //assign mass (clone of user value)
        subunit_bead[index].bx=bxsz;
    }



/*                                                   SUBSUBUNITS                                                        */

    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There is(are) " << count << " subunit(s) in the system." << endl;
    protein.resize(count);
    for(int i=0;i<count;i++)
    {
        vector<int> pog;                                //particles of protein (POG)
        pog.resize((subunit_bead.size()/count));
        crds >> index ;
                for (int j=0; j<pog.size();j++){
                    crds >> pog[j];
                }
        crds >> charge;
        protein[index].id=index;
        for(int j=0;j<(subunit_bead.size()/count);j++)                            //assign particles to subunit and vice versa
        {
            protein[index].itsB.push_back(&subunit_bead[(pog[j])]); //first particle in the subunit, stored in a pointer vector
            subunit_bead[(pog[j])].itsS.push_back(&protein[index]);
            subunit_bead[(pog[j])].unit = protein[index].id;
        }
    }



/*                                                      EDGES                                                       */

    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " edges in the system." << endl;
    subunit_edge.resize(count);
    for(int i=0;i<count;i++)
    {                                                   //Read file and initialize positions, id and charge
        crds >> index >> g1 >> g2 >> type >> length;
        subunit_edge[index].id=index;
        subunit_edge[index].type=type;
        subunit_edge[index].len0=length;
        subunit_edge[index].itsB.push_back(&subunit_bead[g1]);     //first particle in the edge (stored in pointer vector in EDGE class)
        subunit_edge[index].itsB.push_back(&subunit_bead[g2]);
        subunit_bead[g1].itsE.push_back(&subunit_edge[index]);     //the edge on g1 (stored in pointer vector in PARTICLE class)
        subunit_bead[g2].itsE.push_back(&subunit_edge[index]);
		subunit_bead[g1].itsS[0]->itsE.push_back(&subunit_edge[index]);
    }




/*                                                      PARTICLE FACES                                              */


    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " bead faces in the system." << endl;
    subunit_face.resize(count);
    for(int i=0;i<count;i++) {
        crds >> index >> g1 >> g2 >> g3 >> type >>norm;
        subunit_face[index].id = index;
        subunit_face[index].type = type;
        subunit_face[index].itsB.push_back(&subunit_bead[g1]);
        subunit_face[index].itsB.push_back(&subunit_bead[g2]);
        subunit_face[index].itsB.push_back(&subunit_bead[g3]);
        subunit_bead[g1].itsF.push_back(&subunit_face[index]);
        subunit_bead[g2].itsF.push_back(&subunit_face[index]);
        subunit_bead[g3].itsF.push_back(&subunit_face[index]);


        for(int j=0;j<(subunit_edge.size());j++)                 //Finding edge between faces and storing the result for later use
        {
                if ((*subunit_face[index].itsB[0]).id == (*subunit_edge[j].itsB[0]).id && \
                (*subunit_face[index].itsB[1]).id == (*subunit_edge[j].itsB[1]).id) {
                    subunit_face[index].itsE.push_back(&subunit_edge[j]);
                    subunit_edge[j].itsF.push_back(&subunit_face[index]);
                } else if ((*subunit_face[index].itsB[0]).id == (*subunit_edge[j].itsB[1]).id && \
                       (*subunit_face[index].itsB[1]).id == (*subunit_edge[j].itsB[0]).id) {
                    subunit_face[index].itsE.push_back(&subunit_edge[j]);
                    subunit_edge[j].itsF.push_back(&subunit_face[index]);
                }

                if ((*subunit_face[index].itsB[0]).id == (*subunit_edge[j].itsB[0]).id && \
                (*subunit_face[index].itsB[2]).id == (*subunit_edge[j].itsB[1]).id) {
                    subunit_face[index].itsE.push_back(&subunit_edge[j]);
                    subunit_edge[j].itsF.push_back(&subunit_face[index]);
                } else if ((*subunit_face[index].itsB[0]).id == (*subunit_edge[j].itsB[1]).id && \
                       (*subunit_face[index].itsB[2]).id == (*subunit_edge[j].itsB[0]).id) {
                    subunit_face[index].itsE.push_back(&subunit_edge[j]);
                    subunit_edge[j].itsF.push_back(&subunit_face[index]);
                }

                if ((*subunit_face[index].itsB[1]).id == (*subunit_edge[j].itsB[1]).id && \
                (*subunit_face[index].itsB[2]).id == (*subunit_edge[j].itsB[0]).id) {
                    subunit_face[index].itsE.push_back(&subunit_edge[j]);
                    subunit_edge[j].itsF.push_back(&subunit_face[index]);
                } else if ((*subunit_face[index].itsB[1]).id == (*subunit_edge[j].itsB[0]).id && \
                       (*subunit_face[index].itsB[2]).id == (*subunit_edge[j].itsB[1]).id) {
                    subunit_face[index].itsE.push_back(&subunit_edge[j]);
                    subunit_edge[j].itsF.push_back(&subunit_face[index]);
                }
        }                                               //edge_btw function --assigning edges/faces to each other


    }

/*                                                  LJ PAIR                                                         */
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " types of LJ attractive pairs." << endl;
    vector<vector<int> > lj_a(5,vector<int>(count));
    for (int i=0; i<count; i++){
        crds >> index >> g1 >> g2 >> epsilon >> sigma;
        lj_a[0][i] = index;
        lj_a[1][i] = g1;
        lj_a[2][i] = g2;
        lj_a[3][i] = epsilon;
        lj_a[4][i] = sigma;
    }
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " types of LJ repulsive pairs." << endl;
    vector<vector<int> > lj_r(5,vector<int>(count));
    for (int i=0; i<count; i++){
        crds >> index >> g1 >> g2 >> epsilon >> sigma;
        lj_r[0][i] = index;
        lj_r[1][i] = g1;
        lj_r[2][i] = g2;
        lj_r[3][i] = epsilon;
        lj_r[4][i] = sigma;
    }

    crds.close();                                       //closing file


/*                                                  MAKING LJ PAIRS                                                   */

    count=0;
	bool skip_loop = false;
    for (int i=0; i<subunit_bead.size()-1; i++){
        for (int j=i+1; j<subunit_bead.size(); j++){
            if (subunit_bead[i].unit != subunit_bead[j].unit){
                for (int k=0; k<lj_a[0].size(); k++){   //make list of LJ pairs to use in simulation. Categorize attractive / repulsive pairs
                    if (subunit_bead[i].type == lj_a[1][k] && subunit_bead[j].type == lj_a[2][k]){
//                         lj_pairlist.push_back(PAIR(VECTOR3D(0,0,0)));
//                         lj_pairlist[count].type=1;
//                         lj_pairlist[count].itsB.push_back(&subunit_bead[i]);
//                         lj_pairlist[count].itsB.push_back(&subunit_bead[j]);
//                         lj_pairlist[count].epsilon=lj_a[3][k];
//                         lj_pairlist[count].sigma=lj_a[4][k];
						PAIR* new_pair = new PAIR(VECTOR3D(0,0,0));
						subunit_bead[i].itsP.push_back(new_pair);
						subunit_bead[i].itsP.back()->itsB.push_back(&subunit_bead[j]);
						subunit_bead[i].itsP.back()->type = 1;
						subunit_bead[i].itsP.back()->epsilon = lj_a[3][k];
						subunit_bead[i].itsP.back()->sigma = lj_a[4][k];
						subunit_bead[j].itsP.push_back(new_pair);
						subunit_bead[j].itsP.back()->itsB.push_back(&subunit_bead[i]);
						subunit_bead[j].itsP.back()->type = 1;
						subunit_bead[j].itsP.back()->epsilon = lj_a[3][k];
						subunit_bead[j].itsP.back()->sigma = lj_a[4][k];
                        count+=1;
						skip_loop = true;
						break;
                    }
                }
                if (skip_loop == false){
					for (int k=0; k<lj_r[0].size(); k++){
						if (subunit_bead[i].type == lj_r[1][k] && subunit_bead[j].type == lj_r[2][k]){
// 							lj_pairlist.push_back(PAIR(VECTOR3D(0,0,0)));
// 							lj_pairlist[count].type=0;
// 							lj_pairlist[count].itsB.push_back(&subunit_bead[i]);
// 							lj_pairlist[count].itsB.push_back(&subunit_bead[j]);
// 							lj_pairlist[count].epsilon=lj_r[3][k];
// 							lj_pairlist[count].sigma=lj_r[4][k];
							PAIR* new_pair = new PAIR(VECTOR3D(0,0,0));
							subunit_bead[i].itsP.push_back(new_pair);
							subunit_bead[i].itsP.back()->itsB.push_back(&subunit_bead[j]);
							subunit_bead[i].itsP.back()->type = 0;
							subunit_bead[i].itsP.back()->epsilon = lj_r[3][k];
							subunit_bead[i].itsP.back()->sigma = lj_r[4][k];
							subunit_bead[j].itsP.push_back(new_pair);
							subunit_bead[j].itsP.back()->itsB.push_back(&subunit_bead[i]);
							subunit_bead[j].itsP.back()->type = 0;
							subunit_bead[j].itsP.back()->epsilon = lj_r[3][k];
							subunit_bead[j].itsP.back()->sigma = lj_r[4][k];
							count+=1;
							skip_loop = true;
							break;
						}
					}
				}
                if (skip_loop == false) 
				{ 
					cout << "ERROR!!! Beads " << subunit_bead[i].id << " and " << subunit_bead[j].id << " are not an LJ pair!" << endl;
				} else skip_loop = false;
            }
        }
    }

    cout << "Done initializing" << endl;

}




void initialize_outputfile(ofstream &reftraj, ofstream &refofile)
{
    //create files for data analysis
    reftraj << "time" << setw(15) << "kinetic" << setw(15) << "stretching" << setw(15) << "bending" << setw(15) << "LJ" <<
               setw(15) <<  "electrostatics" << setw(15) << "total" \
            << setw(15) << "potential" << setw(15) << "temperature" << endl;



}


void generate_lattice (double capsomere_concentration,unsigned int number_capsomeres, string file_name) {

    ofstream inputfile("input.GEN.out", ios::out);

    ifstream crds;                                      //open coordinates file
    crds.open(file_name.c_str());
    if (!crds) {                                        //check to make sure file is there
        cerr << "ERR: FILE " << file_name << " NOT OPENED. Check directory and/or filename.";
        exit(1);
    }
    long double SIsigma;                                //reading in data
    string dumb;
    double dumber;
    crds >> dumb >> dumb >> dumber >> dumb >> dumb >> dumber >> dumb >> dumb >> SIsigma >> dumb >> dumb >> dumber;

    //Determine box size:
    long double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * 6.022e23)), 1.0 / 3.0);
    VECTOR3D bxsz = VECTOR3D(box_x,box_x,box_x);        //determine box size based on desired concentration

    unsigned int num_fill = int(ceil(pow((double(number_capsomeres)), 1.0 / 3.0)));
    int index = 0;

    int np;                                             //begin reading in template
    crds >> dumb >> dumb >> dumb >> dumb >> dumb >> np;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb >> dumb >> dumb >> dumb;
    long double x,y,z,charge,length,radius, mass;
    int type,name;
    long double part_template[8][np];

    for (int i=0; i < np; i++){
        crds >> name >> x >> y >> z >> type >> charge >> radius >> mass;
        part_template[0][i] = x;
        part_template[1][i] = y;
        part_template[2][i] = z;
        part_template[3][i] = type;
        part_template[4][i] = charge;
        part_template[5][i] = name;
        part_template[6][i] = radius;
        part_template[7][i] = mass;
    }
    int ne;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb >> ne;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb;
    int e1,e2;
    long double edge_template[5][ne];
    for (int i=0; i<ne; i++){
        crds >> name >> e1 >> e2 >> type >> length;
        edge_template[0][i] = name;
        edge_template[1][i] = e1;
        edge_template[2][i] = e2;
        edge_template[3][i] = type;
        edge_template[4][i] = length;
    }
    int nt;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb >> nt;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb >> dumb;
    int t1,t2,t3,normal;
    long double face_template[6][nt];
    for (int i=0; i<nt; i++){
        crds >> name >>t1 >> t2 >> t3 >> type >> normal;
        face_template[0][i] = name;
        face_template[1][i] = t1;
        face_template[2][i] = t2;
        face_template[3][i] = t3;
        face_template[4][i] = normal;
        face_template[5][i] = type;
    }
    int na=0;
    double epsilon,sigma;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb >> dumb >> na;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb;
    long double lj_a_template[5][na];
    for (int i=0; i<na; i++){
        crds >> name >> e1 >> e2 >> epsilon >> sigma;
        lj_a_template[0][i] = name;
        lj_a_template[1][i] = e1;
        lj_a_template[2][i] = e2;
        lj_a_template[3][i] = epsilon;
        lj_a_template[4][i] = sigma;
    }
    int nr=0;
    crds >> dumb >> dumb >> dumb >> dumb >> dumb >> dumb >> nr >> dumb >> dumb >> dumb >> dumb >> dumb;
    long double lj_r_template[5][nr];
    for (int i=0; i<nr; i++){
        crds >> name >> e1 >> e2 >> epsilon >> sigma;
        lj_r_template[0][i] = name;
        lj_r_template[1][i] = e1;
        lj_r_template[2][i] = e2;
        lj_r_template[3][i] = epsilon;
        lj_r_template[4][i] = sigma;
    }
                                                        // Create a master template with copies of original template for each
                                                        // new position of the subunit. This is read in by the original data-reading
                                                        // function, initialize_system
    inputfile << "# Number of Particles = " << np*number_capsomeres << endl << "Coordinates:" << endl;
    inputfile << "index x y z subunit charge type diameter mass" << endl;

    for (int i = 0; i < num_fill; i++) {
        for (int j = 0; j < num_fill; j++) {
            for (int k = 0; k < num_fill; k++) {
                index += 1;
                if (number_capsomeres > (index-1)) {
                    for (int l = 0; l < np; l++) {

                        inputfile << index*np-np+l << setw(15) << setprecision(12) <<(((double)i*bxsz.x*(1/(double)num_fill))+part_template[0][l]) << setw(15)
                                  << setprecision(12) <<(((double)j*bxsz.y*(1/(double)num_fill))+part_template[1][l]) << setw(15)
                                  << setprecision(12) <<(((double)k*bxsz.z*(1/(double)num_fill))+part_template[2][l]) << setw(15) << index << setw(15)
                                  << part_template[4][l] << setw(15) << part_template[3][l] << setw(15) << setprecision(12) <<part_template[6][l]
                                  << setw(20) << part_template[7][l] << endl;

                    }
                }
            }
        }
    }
    index=0;
    inputfile << endl << endl;
    inputfile << "# Number of Subunits = " << number_capsomeres << endl << endl << "Subunits:" << endl << endl;
    inputfile << "index " << setw(15) << "particles in unit" << setw(15) << "charge" << endl << endl;
    for (int i=0; i<number_capsomeres; i++){
        inputfile << i << setw(5);
        for (int j=0; j<np; j++){
            inputfile << (i+1)*np-np+j << setw(5);
        }
        inputfile << 0 << endl;
    }
    inputfile << endl <<endl;
    inputfile << "# Number of Edges = " << ne*number_capsomeres << endl << endl << "Edges:" << endl << endl;
    inputfile << "index" << setw(15) << "e1" << setw(15) << "e2" << setw(15) << "type" << setw(15) << "length" << endl << endl;
    for (int i=0; i<number_capsomeres; i++){
        for (int j=0; j<ne; j++){
            inputfile << index << setw(15) << edge_template[1][j]+np*(double)i << setw(15) <<
                         edge_template[2][j]+np*(double)i << setw(15) << edge_template[3][j] << setw(15) << setprecision(12) <<edge_template[4][j] << endl;
            index +=1;

        }
    }
    index =0;
    inputfile << endl << endl;
    inputfile << "# Number of Triangles = " << nt*number_capsomeres << endl << endl << "Triangles:" <<endl <<endl;
    inputfile << "index" << setw(15) <<	"t1" <<setw(15) << "t2" << setw(15) << "t3" <<setw(15) << "type" << setw(15) << "normal" <<endl <<endl;
    for (int i=0; i<number_capsomeres; i++){
        for (int j=0; j<nt; j++){
            inputfile << index << setw(5) << face_template[1][j]+np*(double)i << setw(5) << face_template[2][j]+np*(double)i
                      << setw(5) << face_template[3][j]+np*(double)i << setw(5) << face_template[5][j] << setw(5) << face_template[4][j] << endl;
            index+=1;
        }
    }
    inputfile << endl << endl << "# Number of LJ Attractions = " << na << endl << endl;
    inputfile << "index" << setw(15) << "p1" << setw(15) << "p2" << setw(15) << "epsilon" << setw(15) << "sigmaHC" << endl;
    for (int i=0; i<na; i++){
        inputfile << lj_a_template[0][i] << setw(15) << lj_a_template[1][i] << setw(15) << lj_a_template[2][i] <<
                     setw(15) << lj_a_template[3][i] << setw(15) << lj_a_template[4][i] << endl;
    }
    inputfile << endl << endl << "# Number of LJ Repulsers = " << nr << endl << endl;
    inputfile << "index" << setw(15) << "p1" << setw(15) << "p2" << setw(15) << "epsilon" << setw(15) << "sigmaHC" << endl;
    for (int i=0; i<nr; i++){
        inputfile << lj_r_template[0][i] << setw(15) << lj_r_template[1][i] << setw(15) << lj_r_template[2][i] <<
                     setw(15) << lj_r_template[3][i] << setw(15) << lj_r_template[4][i] << endl;
    }
}



void initialize_constant_bead_velocities(vector<SUBUNIT> &protein, vector<BEAD> &subunit_bead, double T){
	for (unsigned int i = 0; i < protein.size(); i+=4)
	{
		for (int j = 0; j < protein[i].itsB.size(); j++) {
            protein[i+0].itsB[j]->vel = VECTOR3D(+0.4, -0.2, -0.5);
			protein[i+1].itsB[j]->vel = VECTOR3D(-0.3, +0.3, -0.4);
			protein[i+2].itsB[j]->vel = VECTOR3D(-0.2, +0.5, +0.3);
			protein[i+3].itsB[j]->vel = VECTOR3D(+0.5, -0.4, +0.1);
        }
	}
	
	// average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
    VECTOR3D average_velocity_vector = VECTOR3D(0,0,0);
    for (unsigned int i = 0; i < protein.size(); i++) {
        for (int j=0;j<protein[i].itsB.size();j++) {
            average_velocity_vector = average_velocity_vector + protein[i].itsB[j]->vel;
        }
    }
    
    average_velocity_vector = average_velocity_vector^(1.0/subunit_bead.size());

    // subtract this computed average_velocity_vector from the velocity of each particle to ensure that the total average after this operation is 0
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
        subunit_bead[i].vel = subunit_bead[i].vel - average_velocity_vector;

    // scaling the velocity so that the initial kinetic energy corresponds to the initial temperature supplied by the user.
    double initial_ke = 0.0;
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
    {
        subunit_bead[i].update_kinetic_energy();
        initial_ke += subunit_bead[i].ke;
    }
    double intial_unscaled_T = initial_ke/(1.5*subunit_bead.size());
    double scalefactor = sqrt(T/intial_unscaled_T);

    // scale the velocities
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
    {
        subunit_bead[i].vel = (subunit_bead[i].vel^scalefactor);
        subunit_bead[i].update_kinetic_energy();
    }

    double scaled_ke = 0.0;
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
        scaled_ke = scaled_ke + subunit_bead[i].ke;

    double set_T = scaled_ke / (1.5*subunit_bead.size());
    cout << "initial velocities corrrespond to this set temperature " << set_T << endl;
}


void initialize_bead_velocities(vector<SUBUNIT> &protein, vector<BEAD> &subunit_bead, double T){ //VIKRAM MANY_PARTICLE CODE BLOCK *editted*
    double sigma = sqrt(T);// a rough estimate of how large is the spread in the velocities of the particles at a given temperature
    // Maxwell distribution width
    // assumes all lj atoms have the same mass
    RAND_GEN ugsl;

    // initialized velocities
    // choose uniformly between -0.5 and 0.5
    for (unsigned int i = 0; i < protein.size(); i++)
    {
        double rnumber;
        rnumber = gsl_rng_uniform(ugsl.r);
        double ux = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
        rnumber = gsl_rng_uniform(ugsl.r);
        double uy = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
        rnumber = gsl_rng_uniform(ugsl.r);
        double uz = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
        for (int j=0;j<protein[i].itsB.size();j++) {
            protein[i].itsB[j]->vel = VECTOR3D(ux, uy, uz);
        }
    }

    // average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
    VECTOR3D average_velocity_vector = VECTOR3D(0,0,0);
    for (unsigned int i = 0; i < protein.size(); i++) {
        for (int j=0;j<protein[i].itsB.size();j++) {
            average_velocity_vector = average_velocity_vector + protein[i].itsB[j]->vel;
        }
    }
    average_velocity_vector = average_velocity_vector^(1.0/subunit_bead.size());

    // subtract this computed average_velocity_vector from the velocity of each particle to ensure that the total average after this operation is 0
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
        subunit_bead[i].vel = subunit_bead[i].vel - average_velocity_vector;

    // scaling the velocity so that the initial kinetic energy corresponds to the initial temperature supplied by the user.
    double initial_ke = 0.0;
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
    {
        subunit_bead[i].update_kinetic_energy();
        initial_ke += subunit_bead[i].ke;
    }
    double intial_unscaled_T = initial_ke/(1.5*subunit_bead.size());
    double scalefactor = sqrt(T/intial_unscaled_T);

    // scale the velocities
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
    {
        subunit_bead[i].vel = (subunit_bead[i].vel^scalefactor);
        subunit_bead[i].update_kinetic_energy();
    }

    double scaled_ke = 0.0;
    for (unsigned int i = 0; i < subunit_bead.size(); i++)
        scaled_ke = scaled_ke + subunit_bead[i].ke;

    double set_T = scaled_ke / (1.5*subunit_bead.size());
    cout << "initial velocities corrrespond to this set temperature " << set_T << endl;

    // initial configuration
    ofstream list_initial_velocities("list_initial_velocities.out", ios::out);
    for (unsigned int i = 0; i < protein.size(); i++)
        list_initial_velocities << protein[i].itsB[0]->id << setw(15) << protein[i].itsB[0]->vel.x << "  " << protein[i].itsB[0]->vel.y << "  " <<
                                protein[i].itsB[0]->vel.z << "  " << protein[i].itsB[0]->vel.GetMagnitude() << endl;
    list_initial_velocities.close();

    return;
}
