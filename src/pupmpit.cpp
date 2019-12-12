// pupmpit.cpp: an entry point for the project
//

#include "stdafx.h"
#include "math.h"
#include <memory>
#include "molModel.h"
#include "molDataRead.h"
#include "level.h"
#include "transition.h"
#include "modelPhysPars.h"
#include "PumpCyclesSearchMonteCarlo.h"

int main(int argc, char* argv[])
{
	string filename_lamda;
	string filename_pops;
	string filename_J;
	string phys_cond;
	int seed, close_levels;
	size_t num_of_cycles;

	{
		ifstream fin;
		fin.open("Parameters/pupmpit.txt", ios::in);
		if (!fin.good()) 	// check if we found the file or that console input is not corrupted
			throw runtime_error("can't read Parameters/pupmpit.txt");
			
		string str;

		getline(fin, str);
		getline(fin, str);
		filename_pops = trim(str);

		getline(fin, str);
		getline(fin, str);
		filename_J = trim(str);

		getline(fin, str);
		getline(fin, str);
		phys_cond = trim(str);

		getline(fin, str);
		getline(fin, str);
		filename_lamda = trim(str);

		getline(fin, str);
		seed = readline<int>(fin);

		getline(fin, str);
		num_of_cycles = readline<size_t>(fin);

		getline(fin, str);
		close_levels = readline<size_t>(fin);

		fin.close();
	}

	size_t initial_level = atoi(argv[1]);
	size_t final_level = atoi(argv[2]);

	ifstream fin;
	fin.open(phys_cond, ios::in);
	initialize_modelPhysPars(fin); 	// read physical conditions from file
	fin.close();

	molModel *mod = new molModel();
	{
		molDataRead molReader;
		molReader.read_data(filename_lamda, filename_pops, filename_J, mod); 	// read file with molecular data in LAMDA format
	}

	MonteCarloSearchCycles MC(seed, num_of_cycles, close_levels);
	mod->compute_T(MC.T, MC.isItCollisionalDominated);	// compute net population flow rates (see e.g. Sobolev & Deguchi 1994)
	delete mod;
	MC.pop_flow(initial_level, final_level); 	// search for the maser pumping cycles

	return 0;
}

