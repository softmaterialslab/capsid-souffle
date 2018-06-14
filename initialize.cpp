//
// Created by lauren on 1/25/18.
//

#include <cstdlib>
#include "initialize.h"
#include "rand_gen.h"

using namespace std;

void initialize_system(vector<BEAD> &rgary, vector<EDGE> &rgedge, vector<UNIT> &rgarfield, vector<FACE> &rgface, \
                        VECTOR3D bxsz, BOX & rtardis, vector<PAIR> & rgpair)
{
    /*                                            OPEN FILE                                                            */

    ifstream crds;                                      //open coordinates file
    crds.open("input.GEN.out");
    if (!crds) {                                        //check to make sure file is there
        cerr << "ERR: FILE NOT OPENED. Check directory and/or filename.";
        exit(1);
    }
    int charge, index, g1, g2, g3, norm, type;                   //count, charge, index, bond caps
    unsigned long count;
    long double x, y, z, length,radius, mass, epsilon, sigma;                                //x y z coordinates
    string dummy;                                       //dummy string not imported





/*                                                  BOX                                                             */
    cout << "Gary is in the TARDIS." << endl;
    rtardis.size=bxsz;                                  //tardis box size




/*                                                 PARTICLES                                                        */

    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " beads in the system." << endl;
    for(int i=0;i<count;i++)
    {                                                   //Read file and initialize positions, id and charge
        crds >> index >> x >> y >> z >> dummy >> charge >> type >> radius >> mass;
        rgary.push_back(BEAD(VECTOR3D(x,y,z)));
        rgary[index].id = index;
        rgary[index].q = charge;
        rgary[index].type = type;
        rgary[index].sigma = radius;
        rgary[index].m = mass;                           //assign mass (clone of user value)
        rgary[index].bx=bxsz;
        rgary[index].itsT = (&rtardis);                //assign box/boundaries
        rtardis.itsB.push_back(&rgary[index]);         // puts particles in the box
    }



/*                                                   SUBUNITS                                                        */

    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There is(are) " << count << " subunit(s) in the system." << endl;
    rgarfield.resize(count);
    for(int i=0;i<count;i++)
    {
        vector<int> pog;                                //particles of garfield (POG)
        pog.resize((rgary.size()/count));
        crds >> index ;
                for (int j=0; j<pog.size();j++){
                    crds >> pog[j];
                }
        crds >> charge;
        rgarfield[index].id=index;
        rtardis.itsU.push_back(&rgarfield[index]);
        for(int j=0;j<(rgary.size()/count);j++)                            //assign particles to subunit and vice versa
        {
            rgarfield[index].itsB.push_back(&rgary[(pog[j])]); //first particle in the subunit, stored in a pointer vector
            rgary[(pog[j])].itsU.push_back(&rgarfield[index]);
            rgary[(pog[j])].unit = rgarfield[index].id;
        }
    }



/*                                                      BONDS                                                       */

    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " edges in the system." << endl;
    rgedge.resize(count);
    for(int i=0;i<count;i++)
    {                                                   //Read file and initialize positions, id and charge
        crds >> index >> g1 >> g2 >> type >> length;
        rgedge[index].id=index;
        rgedge[index].type=type;
        rgedge[index].len0=length;
        rgedge[index].itsB.push_back(&rgary[g1]);     //first particle in the bond (stored in pointer vector in BOND class)
        rgedge[index].itsB.push_back(&rgary[g2]);
        rgary[g1].itsE.push_back(&rgedge[index]);     //the bond on g1 (stored in pointer vector in PARTICLE class)
        rgary[g2].itsE.push_back(&rgedge[index]);
    }




/*                                                      PARTICLE FACES                                              */


    crds >> dummy >> dummy >> dummy >> dummy >> dummy >> count >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    cout << "There are " << count << " bead faces in the system." << endl;
    rgface.resize(count);
    for(int i=0;i<count;i++) {
        crds >> index >> g1 >> g2 >> g3 >> type >>norm;
        rgface[index].id = index;
        rgface[index].type = type;
        rgface[index].itsB.push_back(&rgary[g1]);
        rgface[index].itsB.push_back(&rgary[g2]);
        rgface[index].itsB.push_back(&rgary[g3]);
        rgary[g1].itsF.push_back(&rgface[index]);
        rgary[g2].itsF.push_back(&rgface[index]);
        rgary[g3].itsF.push_back(&rgface[index]);


        for(int j=0;j<(rgedge.size());j++)                 //Finding edge between faces and storing the result for later use
        {
                if ((*rgface[index].itsB[0]).id == (*rgedge[j].itsB[0]).id && \
                (*rgface[index].itsB[1]).id == (*rgedge[j].itsB[1]).id) {
                    rgface[index].itsE.push_back(&rgedge[j]);
                    rgedge[j].itsF.push_back(&rgface[index]);
                } else if ((*rgface[index].itsB[0]).id == (*rgedge[j].itsB[1]).id && \
                       (*rgface[index].itsB[1]).id == (*rgedge[j].itsB[0]).id) {
                    rgface[index].itsE.push_back(&rgedge[j]);
                    rgedge[j].itsF.push_back(&rgface[index]);
                }

                if ((*rgface[index].itsB[0]).id == (*rgedge[j].itsB[0]).id && \
                (*rgface[index].itsB[2]).id == (*rgedge[j].itsB[1]).id) {
                    rgface[index].itsE.push_back(&rgedge[j]);
                    rgedge[j].itsF.push_back(&rgface[index]);
                } else if ((*rgface[index].itsB[0]).id == (*rgedge[j].itsB[1]).id && \
                       (*rgface[index].itsB[2]).id == (*rgedge[j].itsB[0]).id) {
                    rgface[index].itsE.push_back(&rgedge[j]);
                    rgedge[j].itsF.push_back(&rgface[index]);
                }

                if ((*rgface[index].itsB[1]).id == (*rgedge[j].itsB[1]).id && \
                (*rgface[index].itsB[2]).id == (*rgedge[j].itsB[0]).id) {
                    rgface[index].itsE.push_back(&rgedge[j]);
                    rgedge[j].itsF.push_back(&rgface[index]);
                } else if ((*rgface[index].itsB[1]).id == (*rgedge[j].itsB[0]).id && \
                       (*rgface[index].itsB[2]).id == (*rgedge[j].itsB[1]).id) {
                    rgface[index].itsE.push_back(&rgedge[j]);
                    rgedge[j].itsF.push_back(&rgface[index]);
                }
        }                                               //bond_btw function --assigning bonds/faces to each other


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
    for (int i=0; i<rgary.size()-1; i++){
        for (int j=i+1; j<rgary.size(); j++){
            if (rgary[i].unit != rgary[j].unit){
                for (int k=0; k<lj_a[0].size(); k++){   //make list of LJ pairs to use in simulation. Categorize attractive / repulsive pairs
                    if (rgary[i].type == lj_a[1][k] && rgary[j].type == lj_a[2][k]){
                        rgpair.push_back(PAIR(VECTOR3D(0,0,0)));
                        rgpair[count].type=1;
                        rgpair[count].itsB.push_back(&rgary[i]);
                        rgpair[count].itsB.push_back(&rgary[j]);
                        rgpair[count].epsilon=lj_a[3][k];
                        rgpair[count].sigma=lj_a[4][k];
                        count+=1;
                        goto breakpoint;
                    }
                }
                for (int k=0; k<lj_r[0].size(); k++){
                    if (rgary[i].type == lj_r[1][k] && rgary[j].type == lj_r[2][k]){
                        rgpair.push_back(PAIR(VECTOR3D(0,0,0)));
                        rgpair[count].type=0;
                        rgpair[count].itsB.push_back(&rgary[i]);
                        rgpair[count].itsB.push_back(&rgary[j]);
                        rgpair[count].epsilon=lj_r[3][k];
                        rgpair[count].sigma=lj_r[4][k];
                        count+=1;
                        goto breakpoint;
                    }
                }
                cout << "ERROR!!! Beads " << rgary[i].id << " and " << rgary[j].id << " are not an LJ pair!" << endl;
            }
            breakpoint: ;
        }
    }

}




void initialize_outputfile(ofstream &reftraj, ofstream &refofile)
{
    //create files for data analysis
    reftraj << "time" << setw(15) << "kinetic" << setw(15) << "stretching" << setw(15) << "bending" << setw(15) << "LJ" <<
               setw(15) <<  "electrostatics" << setw(15) << "total" \
            << setw(15) << "potential" << setw(15) << "temperature" << endl;



}


void allonsy (double capconc,unsigned int numden, string file_name) {

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
    long double box_x = pow((numden * 1000 / (capconc * pow(SIsigma, 3) * 6.022e23)), 1.0 / 3.0);
    VECTOR3D bxsz = VECTOR3D(box_x,box_x,box_x);        //determine box size based on desired concentration

    unsigned int num_fill = int(ceil(pow((double(numden)), 1.0 / 3.0)));
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
    inputfile << "# Number of Particles = " << np*numden << endl << "Coordinates:" << endl;
    inputfile << "index x y z subunit charge type diameter mass" << endl;

    for (int i = 0; i < num_fill; i++) {
        for (int j = 0; j < num_fill; j++) {
            for (int k = 0; k < num_fill; k++) {
                index += 1;
                if (numden > (index-1)) {
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
    inputfile << "# Number of Subunits = " << numden << endl << endl << "Subunits:" << endl << endl;
    inputfile << "index " << setw(15) << "particles in unit" << setw(15) << "charge" << endl << endl;
    for (int i=0; i<numden; i++){
        inputfile << i << setw(5);
        for (int j=0; j<np; j++){
            inputfile << (i+1)*np-np+j << setw(5);
        }
        inputfile << 0 << endl;
    }
    inputfile << endl <<endl;
    inputfile << "# Number of Edges = " << ne*numden << endl << endl << "Edges:" << endl << endl;
    inputfile << "index" << setw(15) << "e1" << setw(15) << "e2" << setw(15) << "type" << setw(15) << "length" << endl << endl;
    for (int i=0; i<numden; i++){
        for (int j=0; j<ne; j++){
            inputfile << index << setw(15) << edge_template[1][j]+np*(double)i << setw(15) <<
                         edge_template[2][j]+np*(double)i << setw(15) << edge_template[3][j] << setw(15) << setprecision(12) <<edge_template[4][j] << endl;
            index +=1;

        }
    }
    index =0;
    inputfile << endl << endl;
    inputfile << "# Number of Triangles = " << nt*numden << endl << endl << "Triangles:" <<endl <<endl;
    inputfile << "index" << setw(15) <<	"t1" <<setw(15) << "t2" << setw(15) << "t3" <<setw(15) << "type" << setw(15) << "normal" <<endl <<endl;
    for (int i=0; i<numden; i++){
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






void initialize_bead_velocities(vector<UNIT> &garfield, vector<BEAD> &gary, double T){ //VIKRAM MANY_PARTICLE CODE BLOCK *editted*
    double sigma = sqrt(T);// a rough estimate of how large is the spread in the velocities of the particles at a given temperature
    // Maxwell distribution width
    // assumes all lj atoms have the same mass
    RAND_GEN ugsl;

    // initialized velocities
    // choose uniformly between -0.5 and 0.5
    for (unsigned int i = 0; i < garfield.size(); i++)
    {
        double rnumber;
        rnumber = gsl_rng_uniform(ugsl.r);
        double ux = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
        rnumber = gsl_rng_uniform(ugsl.r);
        double uy = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
        rnumber = gsl_rng_uniform(ugsl.r);
        double uz = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
        for (int j=0;j<garfield[i].itsB.size();j++) {
            garfield[i].itsB[j]->vel = VECTOR3D(ux, uy, uz);
        }
    }

    // average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
    VECTOR3D average_velocity_vector = VECTOR3D(0,0,0);
    for (unsigned int i = 0; i < garfield.size(); i++) {
        for (int j=0;j<garfield[i].itsB.size();j++) {
            average_velocity_vector = average_velocity_vector + garfield[i].itsB[j]->vel;
        }
    }
    average_velocity_vector = average_velocity_vector^(1.0/gary.size());

    // subtract this computed average_velocity_vector from the velocity of each particle to ensure that the total average after this operation is 0
    for (unsigned int i = 0; i < gary.size(); i++)
        gary[i].vel = gary[i].vel - average_velocity_vector;

    // scaling the velocity so that the initial kinetic energy corresponds to the initial temperature supplied by the user.
    double initial_ke = 0.0;
    for (unsigned int i = 0; i < gary.size(); i++)
    {
        gary[i].update_kinetic_energy();
        initial_ke += gary[i].ke;
    }
    double intial_unscaled_T = initial_ke/(1.5*gary.size());
    double scalefactor = sqrt(T/intial_unscaled_T);

    // scale the velocities
    for (unsigned int i = 0; i < gary.size(); i++)
    {
        gary[i].vel = (gary[i].vel^scalefactor);
        gary[i].update_kinetic_energy();
    }

    double scaled_ke = 0.0;
    for (unsigned int i = 0; i < gary.size(); i++)
        scaled_ke = scaled_ke + gary[i].ke;

    double set_T = scaled_ke / (1.5*gary.size());
    cout << "initial velocities corrrespond to this set temperature " << set_T << endl;

    // initial configuration
    ofstream list_initial_velocities("list_initial_velocities.out", ios::out);
    for (unsigned int i = 0; i < garfield.size(); i++)
        list_initial_velocities << garfield[i].itsB[0]->id << setw(15) << garfield[i].itsB[0]->vel.x << "  " << garfield[i].itsB[0]->vel.y << "  " <<
                                garfield[i].itsB[0]->vel.z << "  " << garfield[i].itsB[0]->vel.GetMagnitude() << endl;
    list_initial_velocities.close();

    return;
}
