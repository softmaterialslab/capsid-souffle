//
// Created by lauren on 1/25/18.
//

#include <cstdlib>
#include "initialize.h"
#include "rand_gen.h"

using namespace std;




void initialize_outputfile(ofstream &reftraj, ofstream &refofile)
{
    //create files for data analysis
    reftraj << "time" << setw(15) << "kinetic" << setw(15) << "stretching" << setw(15) << "bending" << setw(15) << "LJ" <<
               setw(15) <<  "electrostatics" << setw(15) << "total" \
            << setw(15) << "potential" << setw(15) << "temperature" << endl;
}





vector<vector<int> > generate_lattice (double capsomere_concentration,unsigned int number_capsomeres, string file_name, double &bondlength,  double &SIsigma,  double &SImass, double &SItime, vector<BEAD> &subunit_bead, vector<EDGE> &subunit_edge, vector<SUBUNIT> &protein, vector<FACE> &subunit_face) {

    ofstream inputfile("outfiles/input.GEN.out", ios::out);

    ifstream crds;                                      //open coordinates file
    crds.open(file_name.c_str());
    if (!crds) {                                        //check to make sure file is there
        cerr << "ERR: FILE " << file_name << " NOT OPENED. Check directory and/or filename.";
        exit(1);
    }
    //long double SIsigma, SImass, SItime, bondlength;                                //reading in data
    string dummy;
    crds >> dummy >> dummy >> bondlength >> dummy >> dummy >> SImass >> dummy >> dummy >> SIsigma >> dummy >> dummy >> SItime;

    //Determine box size:
    long double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * 6.022e23)), 1.0 / 3.0);
    VECTOR3D bxsz = VECTOR3D(box_x,box_x,box_x);        //determine box size based on desired concentration

    unsigned int num_fill = int(ceil(pow((double(number_capsomeres)), 1.0 / 3.0)));
    unsigned int index = 0;

    unsigned int np;                                             //begin reading in template
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> np;
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    long double x,y,z,charge,length,radius, mass;
    unsigned int type,name;
    long double part_template[8][np];

    for (unsigned int i=0; i < np; i++){
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
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> ne;
    crds >> dummy >> dummy >> dummy >> dummy >> dummy;
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
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> nt;
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
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
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> na;
    crds >> dummy >> dummy >> dummy >> dummy >> dummy;
    //long double lj_a_template[5][na];
	vector<vector<int> > lj_a_template(5,vector<int>(na));
    for (int i=0; i<na; i++){
        crds >> name >> e1 >> e2 >> epsilon >> sigma;
        lj_a_template[0][i] = name;
        lj_a_template[1][i] = e1;
        lj_a_template[2][i] = e2;
        lj_a_template[3][i] = epsilon;
        lj_a_template[4][i] = sigma;
    }
    int nr=0;
    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> nr >> dummy >> dummy >> dummy >> dummy >> dummy;
    int lj_r_template[5][nr];
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

	int myindex = 0;
	
	
    for (unsigned int i = 0; i < num_fill; i++) {
        for (unsigned int j = 0; j < num_fill; j++) {
            for (unsigned int k = 0; k < num_fill; k++) {
                index += 1;
                if (number_capsomeres > (index-1)) {
                    for (unsigned int l = 0; l < np; l++) {
						x = (((double)i*bxsz.x*(1/(double)num_fill))+part_template[0][l]);
						y = (((double)j*bxsz.y*(1/(double)num_fill))+part_template[1][l]);
						z = (((double)k*bxsz.z*(1/(double)num_fill))+part_template[2][l]);
						myindex = index*np-np+l;
						subunit_bead.push_back(BEAD(VECTOR3D(x,y,z)));
						subunit_bead[myindex].id = myindex;
						subunit_bead[myindex].q = part_template[4][l];
						subunit_bead[myindex].type = part_template[3][l];
						subunit_bead[myindex].sigma = part_template[6][l];
						subunit_bead[myindex].m = part_template[7][l];                           //assign mass (clone of user value)
						subunit_bead[myindex].bx=bxsz;
						subunit_bead[myindex].unit = index;
						
						inputfile << subunit_bead[myindex].id << setw(15) << setprecision(12) << subunit_bead[myindex].pos.x << setw(15) << setprecision(12) << subunit_bead[myindex].pos.y << setw(15) << setprecision(12) << subunit_bead[myindex].pos.z << setw(15) << subunit_bead[myindex].unit  << setw(15) << subunit_bead[myindex].q << setw(15) << subunit_bead[myindex].type << setw(15) << subunit_bead[myindex].sigma << setw(20)  << setprecision(12) << subunit_bead[myindex].m << endl;

                    }
                }
            }
        }
    }
    
    
    index=0;
    inputfile << endl << endl;
    inputfile << "# Number of Subunits = " << number_capsomeres << endl << endl << "Subunits:" << endl << endl;
    inputfile << "index " << setw(15) << "particles in unit" << setw(15) << "charge" << endl << endl;
	protein.resize(number_capsomeres);
    for (unsigned int i=0; i<number_capsomeres; i++){

		protein[i].id = i;
		inputfile << protein[i].id << setw(5) ;
        for (unsigned int j=0; j<np; j++){
			
			protein[i].itsB.push_back(&subunit_bead[((i+1)*np-np+j)]); //first particle in the subunit, stored in a pointer vector
            subunit_bead[((i+1)*np-np+j)].itsS.push_back(&protein[i]);
			
			inputfile << subunit_bead[(i+1)*np-np+j].id << setw(5);
			
        }
		inputfile << -11 << endl;
    }
    
    
    inputfile << endl <<endl;
    inputfile << "# Number of Edges = " << ne*number_capsomeres << endl << endl << "Edges:" << endl << endl;
    inputfile << "index" << setw(15) << "e1" << setw(15) << "e2" << setw(15) << "type" << setw(15) << "length" << endl << endl;
	subunit_edge.resize(ne*number_capsomeres);
    for (unsigned int i=0; i<number_capsomeres; i++){
        for (int j=0; j<ne; j++){
			
			int g1 = edge_template[1][j]+np*(double)i;
			int g2 = edge_template[2][j]+np*(double)i;
			
			
			subunit_edge[index].id=index;
			subunit_edge[index].type=edge_template[3][j];
			subunit_edge[index].len0=edge_template[4][j];
			subunit_edge[index].itsB.push_back(&subunit_bead[g1]);     //first particle in the edge (stored in pointer vector in EDGE class)
			subunit_edge[index].itsB.push_back(&subunit_bead[g2]);
			subunit_bead[g1].itsE.push_back(&subunit_edge[index]);     //the edge on g1 (stored in pointer vector in PARTICLE class)
			subunit_bead[g2].itsE.push_back(&subunit_edge[index]);
			subunit_bead[g1].itsS[0]->itsE.push_back(&subunit_edge[index]);
			
			inputfile << subunit_edge[index].id << setw(15) << subunit_edge[index].itsB[0]->id << setw(15) << subunit_edge[index].itsB[1]->id << setw(15) << subunit_edge[index].type << setw(15) << setprecision(12) << subunit_edge[index].len0 << endl;
			
            index +=1;

        }
    }
    
    index =0;
    inputfile << endl << endl;
    inputfile << "# Number of Triangles = " << nt*number_capsomeres << endl << endl << "Triangles:" <<endl <<endl;
    inputfile << "index" << setw(15) <<	"t1" <<setw(15) << "t2" << setw(15) << "t3" <<setw(15) << "type" << setw(15) << "normal" <<endl <<endl;
	subunit_face.resize(nt*number_capsomeres);
    for (unsigned int i=0; i<number_capsomeres; i++){
        for (int j=0; j<nt; j++){
			
			int g1 = face_template[1][j]+np*(double)i;
			int g2 = face_template[2][j]+np*(double)i;
			int g3 = face_template[3][j]+np*(double)i;
			
			subunit_face[index].id = index;
			subunit_face[index].type = face_template[5][j];
			subunit_face[index].itsB.push_back(&subunit_bead[g1]);
			subunit_face[index].itsB.push_back(&subunit_bead[g2]);
			subunit_face[index].itsB.push_back(&subunit_bead[g3]);
			subunit_bead[g1].itsF.push_back(&subunit_face[index]);
			subunit_bead[g2].itsF.push_back(&subunit_face[index]);
			subunit_bead[g3].itsF.push_back(&subunit_face[index]);
			
			inputfile << subunit_face[index].id << setw(5) << subunit_face[index].itsB[0]->id << setw(5) << subunit_face[index].itsB[1]->id << setw(5) << subunit_face[index].itsB[2]->id << setw(5) << subunit_face[index].type << setw(5) << 1 << endl;
			
			
			 for(unsigned int j=0;j<(subunit_edge.size());j++)                 //Finding edge between faces and storing the result for later use
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
        index+=1;
        }
    }
    

	return lj_a_template;
}













void initialize_constant_bead_velocities(vector<SUBUNIT> &protein, vector<BEAD> &subunit_bead, double T){
	for (unsigned int i = 0; i < protein.size(); i+=4)
	{
		for (unsigned int j = 0; j < protein[i].itsB.size(); j++) {
            protein[i+0].itsB[j]->vel = VECTOR3D(+0.4, -0.2, -0.5);
			protein[i+1].itsB[j]->vel = VECTOR3D(-0.3, +0.3, -0.4);
			protein[i+2].itsB[j]->vel = VECTOR3D(-0.2, +0.5, +0.3);
			protein[i+3].itsB[j]->vel = VECTOR3D(+0.5, -0.4, +0.1);
        }
	}
	
	// average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
    VECTOR3D average_velocity_vector = VECTOR3D(0,0,0);
    for (unsigned int i = 0; i < protein.size(); i++) {
        for (unsigned int j=0;j<protein[i].itsB.size();j++) {
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
    //double sigma = sqrt(T);// a rough estimate of how large is the spread in the velocities of the particles at a given temperature
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
        for (unsigned int j=0;j<protein[i].itsB.size();j++) {
            protein[i].itsB[j]->vel = VECTOR3D(ux, uy, uz);
        }
    }

    // average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
    VECTOR3D average_velocity_vector = VECTOR3D(0,0,0);
    for (unsigned int i = 0; i < protein.size(); i++) {
        for (unsigned int j=0;j<protein[i].itsB.size();j++) {
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
    ofstream list_initial_velocities("outfiles/list_initial_velocities.out", ios::out);
    for (unsigned int i = 0; i < protein.size(); i++)
        list_initial_velocities << protein[i].itsB[0]->id << setw(15) << protein[i].itsB[0]->vel.x << "  " << protein[i].itsB[0]->vel.y << "  " <<
                                protein[i].itsB[0]->vel.z << "  " << protein[i].itsB[0]->vel.GetMagnitude() << endl;
    list_initial_velocities.close();

    return;
}
