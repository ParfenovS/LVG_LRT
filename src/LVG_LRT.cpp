// LVG_LRT.cpp: an entry point for the project
//

#include "stdafx.h"
#include <iostream>
#include <memory>
#include <type_traits>
#include "molDataRead.h"
#include "RT_point_iterations.h"
#include "RT_time.h"
#include "RT_time_Newt.h"
#include "modelPhysPars.h"
#include "output.h"

int main(int argc, char* argv[])
{
	// choose whether to solve solve statistical equilibrium equations (with RT_point_iterations class) or integrate kinetic equations for level populations over time (with RT_time_Newt or RT_time class)
	typedef typename std::conditional<TIME_INTEGRATION_OF_KINETIC_EQUATIONS, RT_time_Newt, RT_point_iterations>::type RTclass;
	RTclass *LRT = nullptr;

	int LRT_failure = 1; // 0 - successful program execution; 1 - failure

	try {
		if (argc == 1) {					// read model parameters from file
			ifstream fin; 
			fin.open("Parameters/PhysicalConditions.txt", ios::in);
			initialize_modelPhysPars(fin);	// read physical conditions from file
			fin.close();

			LRT = new RTclass();
		} else {							// read model parameters console input
			initialize_modelPhysPars(cin);	// read physical conditions from console input
			LRT = new RTclass(cin);
		}

		{
			molDataRead molReader;
			molReader.read_data(LRT->filename_lamda, LRT->filename_pops_in, LRT->mol); // reading file with molecular data in LAMDA format
		}

		LRT_failure = LRT->radiative_transfer();	// radiative transfer calculations
		if (LRT_failure == 0) {
			if (argc == 1) output_results(LRT->mol, LRT->filename_pops_out); 	// output results into files
			else output_results(LRT->mol, LRT->filename_pops_out, cout); 		// output results into console
		}
	}
	catch (const exception & e) {
		cerr << "#error: " << e.what() << endl;
		if (LRT != nullptr) delete LRT;
		return LRT_failure;
	}

	delete LRT;
	return LRT_failure;
}