#include <iostream>
#include "md.h"
#include "analysis.h"
#include <boost/filesystem/operations.hpp>

using namespace std;

int analyze_output(int, char **);

int main(int argc, char *argv[]) {
    //if (boost::filesystem::remove_all("outfiles") == 0) perror("Error deleting outfiles directory");
    //else cout << "Pre-existing outfiles directory successfully, replaced with an empty directory." << endl;
    //boost::filesystem::create_directory("outfiles");

	run_simulation(argc,argv);
    analyze_output(argc,argv);
}
