#include <iostream>
#include "md.h"
#include "analysis.h"

using namespace std;

int analyze_output(int, char **);

int main(int argc, char *argv[]) {
	run_simulation(argc,argv);
      analyze_output(argc,argv);
}
