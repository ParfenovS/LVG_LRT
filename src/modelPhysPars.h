#pragma once
#include "auxiliary.h"

struct modelPhysPars						// stores physical conditions
{
	static double Tks;						// gas kinetic temperature [K];
	static double Hdens;					// number density of H2 [1/cm3]; 
	static vector<double> fraction_H2; 		// array with fractional abundances (wrt Hdens) of collisional agents;
	static double NdV;						// specific column density [cm-3 s]
	static double abundance; 				// molecular abundance (wrt H2)
};

double modelPhysPars::Tks;
vector<double> modelPhysPars::fraction_H2;
double modelPhysPars::Hdens;
double modelPhysPars::NdV;
double modelPhysPars::abundance;

template <typename T>
void initialize_modelPhysPars(T & fin)
{
	string str;

	if (!fin.good()) 								// check if we found the file or that the console input is not corrupted
		throw runtime_error("can't read physical conditions; check Parameters/PhysicalConditions.txt or console input");

	getline(fin, str);
	modelPhysPars::Tks = readline<double>(fin); 	// readline function is taken from auxiliary.h

	getline(fin, str);
	modelPhysPars::Hdens = readline<double>(fin);

	getline(fin, str);
	modelPhysPars::NdV = readline<double>(fin);

	getline(fin, str);
	modelPhysPars::abundance = readline<double>(fin);

	getline(fin,str);
	while (getline(fin, str)) {
		str = trim(str); 							// trim is taken from auxiliary.h
		if (str.size() != 0) modelPhysPars::fraction_H2.push_back( stod(str) );
		else break;
	}
}