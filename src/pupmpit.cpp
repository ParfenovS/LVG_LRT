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
	vector <string> filename_lamda;
	vector <string> filename_pops;
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
		phys_cond = trim(str);

		ifstream fin_phys;
		fin_phys.open(phys_cond, ios::in);
		initialize_modelPhysPars(fin_phys); 	// read physical conditions from file
		fin_phys.close();

		getline(fin, str);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			getline(fin, str);
			filename_pops.push_back(trim(str));
		}

		getline(fin, str);
		getline(fin, str);
		filename_J = trim(str);

		getline(fin, str);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			getline(fin, str);
			filename_lamda.push_back(trim(str));
		}

		getline(fin, str);
		seed = readline<int>(fin);

		getline(fin, str);
		num_of_cycles = readline<size_t>(fin);

		getline(fin, str);
		close_levels = readline<int>(fin);

		fin.close();
	}

	if (atoi(argv[1]) < 1) {
		cerr << "The index of molecular species is 1-based and should be > 0\n";
		return 1;
	}
	size_t molSpecies = atoi(argv[1]) - 1;
	
	if (atoi(argv[2]) < 1) {
		cerr << "The index of initial level is 1-based and should be > 0\n";
		return 1;
	}
	size_t initial_level = atoi(argv[2]);

	if (atoi(argv[3]) < 1) {
		cerr << "The index of final level is 1-based and should be > 0\n";
		return 1;
	}
	size_t final_level = atoi(argv[3]);

	MonteCarloSearchCycles MC(seed, num_of_cycles, close_levels);
	for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
		molModel *mol = new molModel(ispec);
		{
			molDataRead molReader;
			molReader.read_data(filename_lamda[ispec], filename_pops[ispec], filename_J, mol); 	// read file with molecular data in LAMDA format
		}
		if (ispec == molSpecies) {
			mol->compute_T(MC.T, MC.isItCollisionalDominated);	// compute net population flow rates (see e.g. Sobolev & Deguchi 1994)
			delete mol;
			break;
		}
		delete mol;
	}
	MC.pop_flow(initial_level, final_level); 	// search for the maser pumping cycles

	return 0;
}

