#pragma once
#include "molModel.h"
#include <stdio.h>

void output_results(molModel *mod, const string & filename) 	// output results of radiative transfer calculations into files
{
	FILE *fout;

	if (filename.length() > 1) {
		fout = fopen(filename.c_str(), "w");
		for (size_t i = 0; i < mod->levels.size(); i++)
			fprintf(fout, "%ld %.17e\n", mod->levels[i].id, mod->levels[i].pop);
		fclose(fout);
	}

	fout = fopen("Results/J_emission.txt", "w");
	for (size_t i = 0; i < mod->rad_trans.size(); i++)
		fprintf(fout, "%ld %ld\t%.17e\n", mod->rad_trans[i].up_level+1, mod->rad_trans[i].low_level+1, mod->rad_trans[i].J);
	fclose(fout);

	fout = fopen("Results/tau.txt", "w");
	fprintf(fout, "#id \tup->low \tfreq., GHz \ttau \tExcitation temperature\n");
	for (size_t i = 0; i < mod->rad_trans.size(); i++) {
		fprintf(fout, "%ld %ld -> %ld\t\t", mod->rad_trans[i].trans_id, mod->rad_trans[i].up_level+1, mod->rad_trans[i].low_level+1);
		fprintf(fout, "%le \t %le \t %le\n", mod->rad_trans[i].nu * 1.e-9, mod->rad_trans[i].tau, mod->rad_trans[i].Tex);
	}
	fclose(fout);

	fout = fopen("Results/Tb.txt", "w");
	fprintf(fout, "#id \tup->low \tfreq., GHz \tBrightness temperature\n");
	for (size_t i = 0; i < mod->rad_trans.size(); i++) {
		fprintf(fout, "%ld %ld -> %ld\t\t", mod->rad_trans[i].trans_id, mod->rad_trans[i].up_level+1, mod->rad_trans[i].low_level+1);
		fprintf(fout, "%le \t %le\n", mod->rad_trans[i].nu * 1.e-9, mod->rad_trans[i].Tbr);
	}
	fclose(fout);

	fout = fopen("Results/masers.txt", "w");
	fprintf(fout, "#id \tup->low \tfreq., GHz \ttau \tBrightness temperature \tExcitation temperature\n");
	for (size_t i = 0; i < mod->rad_trans.size(); i++) {
		if (mod->rad_trans[i].tau < -0.01) {
			fprintf(fout, "%ld\t%ld -> %ld\t", mod->rad_trans[i].trans_id, mod->rad_trans[i].up_level + 1, mod->rad_trans[i].low_level + 1);
			fprintf(fout, "%le \t %le \t", mod->rad_trans[i].nu * 1.e-9, mod->rad_trans[i].tau);
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