//
// Created by lauren on 1/25/18.
//

#include <cstdlib>
#include "initialize.h"
#include "rand_gen.h"

using namespace std;


void initialize_outputfile(ofstream &reftraj, ofstream &refofile) {
   if (world.rank() == 0) {
      //create files for data analysis
      reftraj << "time" << setw(15) << "kinetic" << setw(15) << "stretching" << setw(15) << "bending" << setw(15)
              << "LJ" << setw(15) << "electrostatics" << setw(15) << "total" << setw(15) << "potential" << setw(15)
              << "temperature" << setw(20) << "thermo_potential" << setw(15) << "thermo_kinetic" << setw(15) << "system_velocity_magnitude"<< endl;
   } //if
} // void fxn


vector<vector<int> > generate_lattice(double capsomere_concentration, unsigned int number_capsomeres, string file_name, 
                                      double &bondlength,double &SIsigma, double &SImass, vector<BEAD> &subunit_bead,
                                      vector<EDGE> &subunit_edge, vector<SUBUNIT> &protein, vector<FACE> &subunit_face, 
                                      bool &restartFile, int &restartStep, bool clusters, unsigned int &cluster_size) {
   ofstream inputfile("outfiles/input.GEN.out", ios::out);
   ifstream crds;                                     //open coordinates file
   crds.open(("infiles/"+file_name).c_str());
   if (!crds) {                                       //check to make sure file is there
      if (world.rank() == 0)
         cerr << "ERR: FILE " << file_name << " NOT OPENED. Check directory and/or filename.";
      exit(1);
   }
   
   string file_name2 = "infiles/xyz_T4";                      //File for T3/T4 chunks
   ifstream clust_crds;                                      //open coordinates file
   if (clusters){
      clust_crds.open(file_name2.c_str());
      if (!clust_crds) {                                        //check to make sure file is there
         if (world.rank() == 0)
            cerr << "ERR: FILE " << file_name2 << " NOT OPENED. Check directory and/or filename.";
         exit(1);
      }
   }
   
   string dummy;                                      //reading in data
   crds >> dummy >> dummy >> bondlength >> dummy >> dummy >> SImass >> dummy >> dummy >> SIsigma;
                                                      //Determine box size:
   long double box_x = pow((number_capsomeres * 1000 / (capsomere_concentration * pow(SIsigma, 3) * 6.022e23)),
                            1.0 / 3.0);               // pre factor of 1000 accounts for units
   VECTOR3D bxsz = VECTOR3D(box_x, box_x, box_x);     //determine box size based on desired concentration
   unsigned int num_fill = int(ceil(pow((double(number_capsomeres)), 1.0 / 3.0)));
   unsigned int number_clusters;
  // unsigned int cluster_size;
   if (!clusters) {
      num_fill = int(ceil(pow((double(number_capsomeres)), 1.0 / 3.0)));
   } else if (clusters) {
      clust_crds >> dummy >> dummy >> cluster_size;
      number_clusters = int(double(number_capsomeres)/double(cluster_size));
      if ( number_capsomeres % cluster_size != 0) {
         if (world.rank() == 0)
            cerr << "ERR: Number of capsomeres (" << number_capsomeres << ") is not divisible by cluster_size (" << cluster_size << ")"; 
         exit(1);
      }
      num_fill = int(ceil(pow((double(number_clusters)), 1.0 / 3.0)));
   }
   unsigned int index = 0;
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*                               READING FILE                                                          */
   ////////////////////////////////////////////////////////////////////////////////////////////////////////// 
   unsigned int np;                                   //begin reading in template
   crds >> dummy >> dummy >> dummy >> dummy >> dummy >> np;
   crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
   long double x, y, z, charge, length, radius, mass;
   unsigned int type, name;
   long double part_template[8][np];

   for (unsigned int i = 0; i < np; i++) {
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
   int e1, e2;
   long double edge_template[5][ne];
   for (int i = 0; i < ne; i++) {
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
   int t1, t2, t3, normal;
   long double face_template[6][nt];
   for (int i = 0; i < nt; i++) {
      crds >> name >> t1 >> t2 >> t3 >> type >> normal;
      face_template[0][i] = name;
      face_template[1][i] = t1;
      face_template[2][i] = t2;
      face_template[3][i] = t3;
      face_template[4][i] = normal;
      face_template[5][i] = type;
   }
   int na = 0;
   double epsilon, sigma;
   crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> na;
   crds >> dummy >> dummy >> dummy >> dummy >> dummy;
   vector<vector<int> > lj_a_template(5, vector<int>(na));
   for (int i = 0; i < na; i++) {
      crds >> name >> e1 >> e2 >> epsilon >> sigma;
      lj_a_template[0][i] = name;
      lj_a_template[1][i] = e1;
      lj_a_template[2][i] = e2;
      lj_a_template[3][i] = epsilon;
      lj_a_template[4][i] = sigma;
   }
   vector<vector<long double> > cluster_template(3,vector<long double>(1)); //[3][(np*cluster_size)];   //read in coordinates if there is a cluster file
   if (clusters){
      for (unsigned int i = 0; i < 3; i++){
         cluster_template[i].resize(np*cluster_size);              
      }
      for (unsigned int i = 0; i < (np*cluster_size); i++){
         clust_crds >> dummy >> x >> y >> z;
         cluster_template[0][i] = x;
         cluster_template[1][i] = y;
         cluster_template[2][i] = z;
      }
   }
   /* not reading repulsive LJ
   int nr = 0;
   crds >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> nr >> dummy >> dummy >> dummy >> dummy >> dummy;
   int lj_r_template[5][nr];
   for (int i = 0; i < nr; i++) {
      crds >> name >> e1 >> e2 >> epsilon >> sigma;
      lj_r_template[0][i] = name;
      lj_r_template[1][i] = e1;
      lj_r_template[2][i] = e2;
      lj_r_template[3][i] = epsilon;
      lj_r_template[4][i] = sigma;
   }
   */
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*                               MAKING INPUT.GEN FROM TEMPLATE                                         */
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
   // Create a master template with copies of original template for each
   // new position of the subunit. This is read in by the original data-reading
   // function, initialize_system
   ifstream restart;
   if (world.rank() == 0) {                                 // MAKING BEAD OBJECT
      inputfile << "# Number of Beads = " << np * number_capsomeres << endl << "Coordinates:" << endl;
      inputfile << "index x y z subunit charge type diameter mass" << endl;
   }
  
   vector<string> names = getFileNames("outfiles/"); // get all the files in the directory
   filter(names,"restart"); // filter out everything but restart files
   for (unsigned int i = 0; i < names.size(); i++) {
      if (world.rank() == 0) {
         //cout << names[i] << endl;
      }
   }
   if (names.size() == 0) {
      if (world.rank() == 0) {
         cout << "No restart files found! Starting simulation from scratch." << endl;
      }
      restartFile = false;
   }
   if (restartFile == true) {
      sort(names.begin(), names.end(), numeric_string_compare);
      if (world.rank() == 0) {
       cout << "Restarting from: " << names.back() << endl;
      }
      restart.open(("outfiles/" + names.back()).c_str());
      if (!restart) {                                       //check to make sure file is there
         if (world.rank() == 0)
            cerr << "ERR: RESTART FILE outfiles/restart.out NOT OPENED. Check directory and/or filename.";
         exit(1);
      }
      restart >> dummy >> dummy >> dummy >> dummy >> restartStep;
   }
   int myindex = 0;
   unsigned int loopIndex;
   int subIndex = 0;
   unsigned int tempIndex = 0;
   double vel_x, vel_y, vel_z;
   unsigned int check;

   for (unsigned int i = 0; i < num_fill; i++) {
      for (unsigned int j = 0; j < num_fill; j++) {
         for (unsigned int k = 0; k < num_fill; k++) {
               index += 1;
               if (!clusters) check = index - 1;
               if (clusters) check = (index -1)*cluster_size;
               if (number_capsomeres > check) {
                  if (clusters) loopIndex = np * cluster_size;
                  if (!clusters) loopIndex = np;
                  for (unsigned int l = 0; l < loopIndex; l++) {
                     if (restartFile == false && clusters == false) {
                        x = (((double) i * bxsz.x * (1 / (double) num_fill)) + part_template[0][l]);
                        y = (((double) j * bxsz.y * (1 / (double) num_fill)) + part_template[1][l]);
                        z = (((double) k * bxsz.z * (1 / (double) num_fill)) + part_template[2][l]);
                     } else if (restartFile == true) {
                        restart >> dummy >> vel_x >> vel_y >> vel_z >> x >> y >> z;
                     } else if (restartFile == false && clusters == true) {
                        x = (((double) i * bxsz.x * (1 / (double) num_fill)) + cluster_template[0][l]);
                        y = (((double) j * bxsz.x * (1 / (double) num_fill)) + cluster_template[1][l]);
                        z = (((double) k * bxsz.x * (1 / (double) num_fill)) + cluster_template[2][l]);
                     }
                     if (!clusters) {     // If there are no clusters...
                        myindex = index * np - np + l;
                        subunit_bead.push_back(BEAD(VECTOR3D(x, y, z)));
                        subunit_bead[myindex].id = myindex;
                        subunit_bead[myindex].q = part_template[4][l];
                        subunit_bead[myindex].type = part_template[3][l];
                        subunit_bead[myindex].sigma = part_template[6][l];
                        subunit_bead[myindex].m = part_template[7][l];//assign mass (clone of user value)
                        subunit_bead[myindex].bx = bxsz;
                        subunit_bead[myindex].hbx = bxsz ^ 0.5;
                        subunit_bead[myindex].unit = index;
                     }
                     if (clusters) {      // If there are clusters...
                        if (tempIndex == np) tempIndex = 0;
                        if (myindex % np == 0) subIndex ++;
                        subunit_bead.push_back(BEAD(VECTOR3D(x, y, z)));
                        subunit_bead[myindex].id = myindex;
                        subunit_bead[myindex].q = part_template[4][tempIndex];
                        subunit_bead[myindex].type = part_template[3][tempIndex];
                        subunit_bead[myindex].sigma = part_template[6][tempIndex];
                        subunit_bead[myindex].m = part_template[7][tempIndex];//assign mass (clone of user value)
                        subunit_bead[myindex].bx = bxsz;
                        subunit_bead[myindex].hbx = bxsz ^ 0.5;
                        subunit_bead[myindex].unit = subIndex;
                        tempIndex ++;
                     }
                     if (restartFile == true) {
                        subunit_bead[myindex].vel = VECTOR3D(vel_x, vel_y, vel_z);
                     }
                     if (world.rank() == 0) {
                           inputfile << subunit_bead[myindex].id << setw(15) << setprecision(12)
                                    << subunit_bead[myindex].pos.x << setw(25) << setprecision(12)
                                    << subunit_bead[myindex].pos.y << setw(25) << setprecision(12)
                                    << subunit_bead[myindex].pos.z << setw(25) << subunit_bead[myindex].unit
                                    << setw(15)
                                    << subunit_bead[myindex].q << setw(15) << subunit_bead[myindex].type << setw(15)
                                    << subunit_bead[myindex].sigma << setw(20) << setprecision(12)
                                    << subunit_bead[myindex].m << endl;
                     } //if
                     if (clusters) myindex ++;
                  } //for l
               } //if
         } // for k
      } //for j
   } //for i
      //finding center bead for later use
   VECTOR3D centroid = VECTOR3D(0,0,0);
   for (unsigned int i = 0; i < np; i++) {
      centroid += subunit_bead[i].pos; 
   }
   centroid = centroid / np; // average the positions of all the beads in a subunit in order to get centroid of cluster
   long double oldDistance = 1000;
   long double newDistance;
   int centerBeadID = 0;
   for (unsigned int i = 0; i < np; i++) {
      VECTOR3D r_vec = (subunit_bead[i].pos - centroid);
      VECTOR3D hbox = bxsz / 2;
      if (r_vec.x > hbox.x) r_vec.x -= bxsz.x;
      else if (r_vec.x < -hbox.x) r_vec.x += bxsz.x;
      if (r_vec.y > hbox.y) r_vec.y -= bxsz.y;
      else if (r_vec.y < -hbox.y) r_vec.y += bxsz.y;
      if (r_vec.z > hbox.z) r_vec.z -= bxsz.z;
      else if (r_vec.z < -hbox.z) r_vec.z += bxsz.z;
      newDistance = r_vec.GetMagnitudeSquared();
      if (newDistance < oldDistance) {    //find bead closest to the centroid
         oldDistance = newDistance;
         centerBeadID = subunit_bead[i].id; 
      }
   }
   index = 0;
   if (world.rank() == 0) {                           // MAKING SUBUNIT OBJECTS
      inputfile << endl << endl;
      inputfile << "# Number of Subunits = " << number_capsomeres << endl << endl << "Subunits:" << endl << endl;
      inputfile << "index " << setw(15) << "particles in unit" << setw(15) << "charge" << endl << endl;
   }
   protein.resize(number_capsomeres);
   for (unsigned int i = 0; i < number_capsomeres; i++) {
      protein[i].id = i;
      protein[i].centerBead = & subunit_bead[((i + 1) * np - np + centerBeadID)];
      if (world.rank() == 0)
         inputfile << protein[i].id << setw(5);
      for (unsigned int j = 0; j < np; j++) {
         protein[i].itsB.push_back(&subunit_bead[((i + 1) * np - np + j)]); //first particle in the subunit, stored in a pointer vector
         subunit_bead[((i + 1) * np - np + j)].itsS.push_back(&protein[i]);
         if (world.rank() == 0)
               inputfile << subunit_bead[(i + 1) * np - np + j].id << setw(5);
      } //for j
      if (world.rank() == 0)
         inputfile << -11 << endl;
   } //for i

   if (world.rank() == 0) {                           // MAKING EDGE OBJECTS
      inputfile << endl << endl;
      inputfile << "# Number of Edges = " << ne * number_capsomeres << endl << endl << "Edges:" << endl << endl;
      inputfile << "index" << setw(15) << "e1" << setw(15) << "e2" << setw(15) << "type" << setw(15) << "length"
                << endl << endl;
   }
   subunit_edge.resize(ne * number_capsomeres);
   for (unsigned int i = 0; i < number_capsomeres; i++) {
      for (int j = 0; j < ne; j++) {
         int g1 = edge_template[1][j] + np * (double) i;
         int g2 = edge_template[2][j] + np * (double) i;

         subunit_edge[index].id = index;
         subunit_edge[index].type = edge_template[3][j];
         subunit_edge[index].len0 = edge_template[4][j];
         subunit_edge[index].itsB.push_back(&subunit_bead[g1]);     //first particle in the edge (stored in pointer vector in EDGE class)
         subunit_edge[index].itsB.push_back(&subunit_bead[g2]);
         subunit_bead[g1].itsE.push_back(&subunit_edge[index]);     //the edge on g1 (stored in pointer vector in PARTICLE class)
         subunit_bead[g2].itsE.push_back(&subunit_edge[index]);
         subunit_bead[g1].itsS[0]->itsE.push_back(&subunit_edge[index]);

         if (world.rank() == 0) {
            inputfile << subunit_edge[index].id << setw(15) << subunit_edge[index].itsB[0]->id << setw(15)
                      << subunit_edge[index].itsB[1]->id << setw(15) << subunit_edge[index].type << setw(15)
                      << setprecision(12) << subunit_edge[index].len0 << endl;
         } //if
         index += 1;
      } //for j
   } //for i

   index = 0;
   if (world.rank() == 0) {                          // MAKING FACE OBJECTS
      inputfile << endl << endl;
      inputfile << "# Number of Faces = " << nt * number_capsomeres << endl << endl << "Faces:" << endl << endl;
      inputfile << "index" << setw(15) << "t1" << setw(15) << "t2" << setw(15) << "t3" << setw(15) << "type"
               << setw(15) << "normal" << endl << endl;
   }
   subunit_face.resize(nt * number_capsomeres);
   for (unsigned int i = 0; i < number_capsomeres; i++) {
      for (int j = 0; j < nt; j++) {
         int g1 = face_template[1][j] + np * (double) i;
         int g2 = face_template[2][j] + np * (double) i;
         int g3 = face_template[3][j] + np * (double) i;

         subunit_face[index].id = index;
         subunit_face[index].type = face_template[5][j];
         subunit_face[index].itsB.push_back(&subunit_bead[g1]);
         subunit_face[index].itsB.push_back(&subunit_bead[g2]);
         subunit_face[index].itsB.push_back(&subunit_bead[g3]);
         subunit_bead[g1].itsF.push_back(&subunit_face[index]);
         subunit_bead[g2].itsF.push_back(&subunit_face[index]);
         subunit_bead[g3].itsF.push_back(&subunit_face[index]);
         protein[i].itsF.push_back(&subunit_face[index]);

         if (world.rank() == 0) {
               inputfile << subunit_face[index].id << setw(5) << subunit_face[index].itsB[0]->id << setw(5)
                        << subunit_face[index].itsB[1]->id << setw(5) << subunit_face[index].itsB[2]->id << setw(5)
                        << subunit_face[index].type << setw(5) << 1 << endl;
         }
         for (unsigned int j = 0; j <(subunit_edge.size()); j++) //Finding edge between faces and storing the result for later use
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
         }                                             
         index += 1;
      } //for j
   } //for i
   
   /*
   for (unsigned int i = 0; i < subunit_bead.size(); i++) {
      if (subunit_bead[i].unit % 3 == 0) {
         if (subunit_bead[i].type == 1 || subunit_bead[i].type == 11) {
            subunit_bead[i].type = 10;
         } else if (subunit_bead[i].type == 2 || subunit_bead[i].type == 9) {
            subunit_bead[i].type = 8;
         }
      }
   } */

   return lj_a_template; //not a member variable, so returned from this fxn
} //generate lattice fxn



void initialize_constant_bead_velocities(vector<SUBUNIT> &protein, vector<BEAD> &subunit_bead, double T) {
   for (unsigned int i = 0; i < protein.size(); i += 4) {
      for (unsigned int j = 0; j < protein[i].itsB.size(); j++) {
         protein[i + 0].itsB[j]->vel = VECTOR3D(+0.4, -0.2, -0.5);
         protein[i + 1].itsB[j]->vel = VECTOR3D(-0.3, +0.3, -0.4);
         protein[i + 2].itsB[j]->vel = VECTOR3D(-0.2, +0.5, +0.3);
         protein[i + 3].itsB[j]->vel = VECTOR3D(+0.5, -0.4, +0.1);
      }
   }
   // average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
   VECTOR3D average_velocity_vector = VECTOR3D(0, 0, 0);
   for (unsigned int i = 0; i < protein.size(); i++) {
      for (unsigned int j = 0; j < protein[i].itsB.size(); j++) {
         average_velocity_vector = average_velocity_vector + protein[i].itsB[j]->vel;
      }
   }

   average_velocity_vector = average_velocity_vector ^ (1.0 / subunit_bead.size());

   // subtract this computed average_velocity_vector from the velocity of each particle to ensure that the total average after this operation is 0
   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      subunit_bead[i].vel = subunit_bead[i].vel - average_velocity_vector;

   // scaling the velocity so that the initial kinetic energy corresponds to the initial temperature supplied by the user.
   double initial_ke = 0.0;
   for (unsigned int i = 0; i < subunit_bead.size(); i++) {
      subunit_bead[i].update_kinetic_energy();
      initial_ke += subunit_bead[i].ke;
   }
   double intial_unscaled_T = initial_ke / (1.5 * subunit_bead.size());
   double scalefactor = sqrt(T / intial_unscaled_T);

   // scale the velocities
   for (unsigned int i = 0; i < subunit_bead.size(); i++) {
      subunit_bead[i].vel = (subunit_bead[i].vel ^ scalefactor);
      subunit_bead[i].update_kinetic_energy();
   }

   double scaled_ke = 0.0;
   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      scaled_ke = scaled_ke + subunit_bead[i].ke;

   double set_T = scaled_ke / (1.5 * subunit_bead.size());
   if (world.rank() == 0)
      cout << "initial velocities corrrespond to this set temperature " << set_T << endl;
}



void initialize_bead_velocities(vector<SUBUNIT> &protein, vector<BEAD> &subunit_bead,
                                double T, bool cluster, unsigned int cluster_size) { 
   //double sigma = sqrt(T);// a rough estimate of how large is the spread in the velocities of the particles at a given temperature
   // Maxwell distribution width
   // assumes all lj atoms have the same mass
   RAND_GEN ugsl;

   // initialized velocities
   // choose uniformly between -0.5 and 0.5
   unsigned int maxValue;
   if (!cluster) maxValue = protein.size();
   if (cluster) maxValue = int(double(protein.size())/double(cluster_size));
      
   for (unsigned int i = 0; i < maxValue; i++) {
      double rnumber;
      rnumber = gsl_rng_uniform(ugsl.r);
      double ux = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
      rnumber = gsl_rng_uniform(ugsl.r);
      double uy = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
      rnumber = gsl_rng_uniform(ugsl.r);
      double uz = 0.5 * (rnumber) + (-0.5) * (1 - rnumber); // scale to get the number between -0.5 and 0.5
      if (cluster){
         for (unsigned int j = 0; j < cluster_size; j++) {
            for (unsigned int k = 0; k < protein[( (i*cluster_size) + j)].itsB.size(); k++){
               protein[( (i*cluster_size) + j)].itsB[k]->vel = VECTOR3D(ux, uy, uz);
            }
         }
      }
      if (!cluster){
         for (unsigned int j = 0; j < protein[i].itsB.size(); j++) {
            protein[i].itsB[j]->vel = VECTOR3D(ux, uy, uz);
         }
      }
   }

   // average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
   VECTOR3D average_velocity_vector = VECTOR3D(0, 0, 0);
   if (!cluster){
      for (unsigned int i = 0; i < protein.size(); i++) {
         for (unsigned int j = 0; j < protein[i].itsB.size(); j++) {
            average_velocity_vector = average_velocity_vector + protein[i].itsB[j]->vel;
         }
      }
   }
   if (cluster){
      for (unsigned int i = 0; i < maxValue; i++) {
         for (unsigned int j = 0; j < cluster_size; j++) {
            for (unsigned int k = 0; k < protein[( (i*cluster_size) + j)].itsB.size(); k++){
               average_velocity_vector = average_velocity_vector + protein[( (i*cluster_size) + j)].itsB[k]->vel;
            }
         }
      }
   }
   
   average_velocity_vector = average_velocity_vector ^ (1.0 / subunit_bead.size());

   // subtract this computed average_velocity_vector from the velocity of each particle to ensure that the total average after this operation is 0
   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      subunit_bead[i].vel = subunit_bead[i].vel - average_velocity_vector;

   // scaling the velocity so that the initial kinetic energy corresponds to the initial temperature supplied by the user.
   double initial_ke = 0.0;
   for (unsigned int i = 0; i < subunit_bead.size(); i++) {
      subunit_bead[i].update_kinetic_energy();
      initial_ke += subunit_bead[i].ke;
   }
   double intial_unscaled_T = initial_ke / (1.5 * subunit_bead.size());
   double scalefactor = sqrt(T / intial_unscaled_T);

   // scale the velocities
   for (unsigned int i = 0; i < subunit_bead.size(); i++) {
      subunit_bead[i].vel = (subunit_bead[i].vel ^ scalefactor);
      subunit_bead[i].update_kinetic_energy();
   }

   double scaled_ke = 0.0;
   for (unsigned int i = 0; i < subunit_bead.size(); i++)
      scaled_ke = scaled_ke + subunit_bead[i].ke;

   double set_T = scaled_ke / (1.5 * subunit_bead.size());
   if (world.rank() == 0)
      cout << "initial velocities corrrespond to this set temperature " << set_T << endl;

   // initial configuration
   ofstream list_initial_velocities("outfiles/list_initial_velocities.out", ios::out);
   if (world.rank() == 0) {
      for (unsigned int i = 0; i < protein.size(); i++)
         list_initial_velocities << protein[i].itsB[0]->id << setw(15) << protein[i].itsB[0]->vel.x << "  "
                                 << protein[i].itsB[0]->vel.y << "  " <<
                                 protein[i].itsB[0]->vel.z << "  " << protein[i].itsB[0]->vel.GetMagnitude() << endl;
      list_initial_velocities.close();
   }
   return;
}

vector<BEAD> generate_big_beads(int num_big_beads, double sigma, const vector<BEAD>& existing_beads,
                               const VECTOR3D& box_size, double mass, int type, double charge,
                               unsigned int max_attempts = 10000, bool restartflag = false) {
    vector<BEAD> big_beads;
    srand(0);  // Seed random number generator

    const double box_x = box_size.x;
    const double hbx = box_x * 0.5;  // Half box size (assuming cubic box)
    const double R_squared = (sigma) * (sigma);  // Minimum distance squared between big beads

    if (restartflag == true) {
      ifstream restart;
      double vel_x, vel_y, vel_z, x, y, z;
      vector<string> names = getFileNames("outfiles/");
      filter(names, "rb"); // Filters to keep only files containing "rb"
      sort(names.begin(), names.end(), numeric_string_compare);
      restart.open(("outfiles/" + names.back()).c_str());
      for(int i = 0; i < num_big_beads; i++) {
         restart >> vel_x >> vel_y >> vel_z >> x >> y >> z;
         cout << vel_x << ", " << vel_y << ", " << vel_z << endl;
         BEAD new_bead(VECTOR3D(x, y, z));
         new_bead.sigma = sigma;
         new_bead.type = type;
         new_bead.m = mass;
         new_bead.q = charge;
         new_bead.bx = box_size;
         new_bead.hbx = VECTOR3D(hbx, hbx, hbx);
         new_bead.id = existing_beads.size() + big_beads.size();
         new_bead.eforce = VECTOR3D(0,0,0);
         new_bead.ljforce = VECTOR3D(0,0,0);
         new_bead.vel = VECTOR3D(vel_x,vel_y,vel_z);

         big_beads.push_back(new_bead);
      }
      return big_beads;
   }

    for(int i = 0; i < num_big_beads; i++) {
        bool placed = false;
        unsigned int attempts = 0;

        while(!placed && attempts < max_attempts) {
            // Generate random position within box
            double x = (static_cast<double>(rand())/RAND_MAX - 0.5) * box_size.x;
            double y = (static_cast<double>(rand())/RAND_MAX - 0.5) * box_size.y;
            double z = (static_cast<double>(rand())/RAND_MAX - 0.5) * box_size.z;

            bool collision = false;

            // Check against existing beads
           for(const BEAD& bead : existing_beads) {
              double dx = x - bead.pos.x;
              dx -= box_size.x * floor(dx/box_size.x + 0.5);  // Periodic correction

              double dy = y - bead.pos.y;
              dy -= box_size.y * floor(dy/box_size.y + 0.5);

              double dz = z - bead.pos.z;
              dz -= box_size.z * floor(dz/box_size.z + 0.5);

              double dist_sq = dx*dx + dy*dy + dz*dz;
              double min_dist = bead.sigma / 2 + sigma / 2;

              if(dist_sq < min_dist*min_dist) {
                 collision = true;
                 break;
              }
           }

            // Check against other big beads
           if(!collision) {
              for(const BEAD& big_bead : big_beads) {
                 double dx = x - big_bead.pos.x;
                 dx -= box_size.x * floor(dx/box_size.x + 0.5);

                 double dy = y - big_bead.pos.y;
                 dy -= box_size.y * floor(dy/box_size.y + 0.5);

                 double dz = z - big_bead.pos.z;
                 dz -= box_size.z * floor(dz/box_size.z + 0.5);

                 double dist_sq = dx*dx + dy*dy + dz*dz;

                 if(dist_sq < R_squared) {
                    collision = true;
                    break;
                 }
              }
           }

            if(!collision) {
                BEAD new_bead(VECTOR3D(x, y, z));
                new_bead.sigma = sigma;
                new_bead.type = type;
                new_bead.m = mass;
                new_bead.q = charge;
                new_bead.bx = box_size;
                new_bead.hbx = VECTOR3D(hbx, hbx, hbx);
                new_bead.id = existing_beads.size() + big_beads.size();
                new_bead.vel = VECTOR3D(0,0,0);

                big_beads.push_back(new_bead);
                placed = true;
            }
            attempts++;
        }

        if(attempts >= max_attempts) {
            cerr << "Error: Failed to place big bead " << i+1
                 << " after " << max_attempts << " attempts" << endl;
            exit(1);
        }
    }

    return big_beads;
}