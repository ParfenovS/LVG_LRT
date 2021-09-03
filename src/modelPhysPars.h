#pragma once
#include "auxiliary.h"

struct modelPhysPars						// stores physical conditions
{
	static size_t nSpecies;					// number of molecular species
	static size_t collPartCounter;			// index used to assign fractional abundances of collisional agents to a given molecular species
	static double Tks;						// gas kinetic temperature [K];
	static double Hdens;					// number density of H2 [1/cm3]; 
	static vector<double> fraction_H2;		// array with fractional abundances (wrt Hdens) of collisional agents;
	static vector<double> NdV;				// specific column density [cm-3 s]
	static vector<double> abundance;		// molecular abundance (wrt H2)
};

size_t modelPhysPars::nSpecies;
size_t modelPhysPars::collPartCounter;
double modelPhysPars::Tks;
double modelPhysPars::Hdens;
vector<double> modelPhysPars::fraction_H2;
vector<double> modelPhysPars::NdV;
vector<double> modelPhysPars::abundance;

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
	modelPhysPars::nSpecies = readline<size_t>(fin);

	getline(fin, str);
	for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) modelPhysPars::NdV.push_back(readline<double>(fin));

	getline(fin, str);
	for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) modelPhysPars::abundance.push_back(readline<double>(fin));

	getline(fin,str);
	while (getline(fin, str)) {
		str = trim(str); 							// trim is taken from auxiliary.h
		if (str.size() != 0) modelPhysPars::fraction_H2.push_back( stod(str) );
		else break;
	}

	modelPhysPars::collPartCounter = 0;
}
