#include <iostream>
#include <time.h>

#ifdef EOLC_MOSEK
#include "external\SolverWrappers\Mosek\QuadProgMosek.h"
#endif

#ifdef EOLC_GUROBI
#include "external\SolverWrappers\Gurobi\Gurobi.h"
#endif

#include "conversions.h"

#include "runner.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
	if (argc < 3) {
		cout << "Usage: " << endl;
		cout << "	" << argv[0] << "<general setting> <simulation settings>" << endl;
		cout << "where the settings arguments are json file" << endl;
		cout << endl << "Refer to the examples and json templates for further details" << endl;
		return 0;
	}

	srand(time(NULL));
    
	start_running(argv[1], argv[2]);

    return 0;
}