#pragma once
#include "molModel.h"
#include <stdio.h>

void output_results(molModel *mod, const string & filename) 	// output results of radiative transfer calculations into files
{
	FILE *fout;
	errno_t err;

	if (filename.length() > 1) {
		err = fopen_s(&fout, filename.c_str(), "w");
		if (err !=0 ) throw runtime_error("can't open file to save populations");
		for (size_t i = 0; i < mod->levels.size(); i++)
			fprintf(fout, "%ld %.17e\n", mod->levels[i].id, mod->levels[i].pop);
		fclose(fout);
	}

	err = fopen_s(&fout, "Results/J_emission.txt", "w");
	if (err != 0) throw runtime_error("can't open file to save mean intensity");
	for (size_t i = 0; i < mod->rad_trans.size(); i++)
		fprintf(fout, "%ld %ld\t%.17e\n", mod->rad_trans[i].up_level+1, mod->rad_trans[i].low_level+1, mod->rad_trans[i].J);
	fclose(fout);

	err = fopen_s(&fout, "Results/tau.txt", "w");
	if (err != 0) throw runtime_error("can't open file to save optical depths");
	fprintf(fout, "#id \t up->low \t freq., GHz \t tau \t Excitation temperature\n");
	for (size_t i = 0; i < mod->rad_trans.size(); i++) {
		fprintf(fout, "%ld \t %ld -> %ld \t ", mod->rad_trans[i].trans_id, mod->rad_trans[i].up_level+1, mod->rad_trans[i].low_level+1);
		fprintf(fout, "%le \t %le \t %le\n", mod->rad_trans[i].nu * 1.e-9, mod->rad_trans[i].tau, mod->rad_trans[i].Tex);
	}
	fclose(fout);

	err = fopen_s(&fout, "Results/Tb.txt", "w");
	if (err != 0) throw runtime_error("can't open file to save brightness temperature");
	fprintf(fout, "#id \t up->low \t freq., GHz \t Brightness temperature\n");
	for (size_t i = 0; i < mod->rad_trans.size(); i++) {
		fprintf(fout, "%ld \t %ld -> %ld \t ", mod->rad_trans[i].trans_id, mod->rad_trans[i].up_level+1, mod->rad_trans[i].low_level+1);
		fprintf(fout, "%le \t %le\n", mod->rad_trans[i].nu * 1.e-9, mod->rad_trans[i].Tbr);
	}
	fclose(fout);

	err = fopen_s(&fout, "Results/masers.txt", "w");
	if (err != 0) throw runtime_error("can't open file to save maser transitions");
	fprintf(fout, "#id \t up->low \t freq., GHz \t tau \t Brightness temperature \t Excitation temperature\n");
	for (size_t i = 0; i < mod->rad_trans.size(); i++) {
		if (mod->rad_trans[i].tau < -0.01) {
			fprintf(fout, "%ld \t %ld -> %ld \t ", mod->rad_trans[i].trans_id, mod->rad_trans[i].up_level + 1, mod->rad_trans[i].low_level + 1);
			fprintf(fout, "%le \t %le \t ", mod->rad_trans[i].nu * 1.e-9, mod->rad_trans[i].tau);
			fprintf(fout, "%le \t %le\n", mod->rad_trans[i].Tbr, mod->rad_trans[i].Tex);
		}
	}
	fclose(fout);
}

void output_results(molModel *mod, const string & filename, ostream & cout) 	// output results of radiative transfer calculations into console output
{
	cout << "#id \tup->low \tfreq., GHz \ttau \tBrightness temperature \tExcitation temperature\n";
	for (size_t i = 0; i < mod->rad_trans.size(); i++) {
		cout << mod->rad_trans[i].trans_id << "\t" << mod->rad_trans[i].up_level + 1 << "->" << mod->rad_trans[i].low_level + 1 << "\t\t";
		cout << mod->rad_trans[i].nu * 1.e-9 << "\t\t" << mod->rad_trans[i].tau << "\t\t";
		cout << mod->rad_trans[i].Tbr << "\t\t" << mod->rad_trans[i].Tex << endl;
	}

	if (filename.length() > 1) {
		ofstream fout;
		fout.open(filename.c_str(), ios::out);
		for (size_t i = 0; i < mod->levels.size(); i++)
			fout << mod->levels[i].id << "\t" << mod->levels[i].pop << endl;
		fout.close();
	}
}
