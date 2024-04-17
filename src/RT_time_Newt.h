#pragma once
#include "RT.h"

class RT_time_Newt : public RT		// integrates kinetic equations for level populations over time; the non-linear system of equations at each time step is solved with Newton method
{
private:

	void prepare_results_for_output(beta_LVG & LVG_beta, const double & time)
	{
		// preparing quantities for output
		double dummy_S = 0.0;
		double dummy_beta = 0.0;
		double dummy_betaS = 0.0;
		double dummy_kabs = 0.0;
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				compute_tau(i, &mols[ispec]); 		// computing final optical depths that can be used for output
				compute_J_S_beta(&mols[ispec], i, LVG_beta, dummy_S, dummy_beta, dummy_betaS, dummy_kabs); 	// computing final mean intensities that can be used for output
				compute_Tex(i, &mols[ispec]); 		// computing excitation temperature that can be used for output
				compute_brightness_temperature(i, time, &mols[ispec]); 	// computing brightness temperature and intensity of the emission
			}
		}
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				// the output optical depth is an optical depth along the line of sight, thus one need to take into account beaming
				mols[ispec].rad_trans[i].tau *= beamH;
			}
		}
	}

	void write_rad_trans_data_into_binary_file(beta_LVG & LVG_beta, const double & time, ofstream & binTbrfile, ofstream & binTexfile, ofstream & bintaufile)
	{
		binTbrfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
		binTexfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
		bintaufile.write(reinterpret_cast<const char*>(&time), sizeof(double));
		prepare_results_for_output(LVG_beta, time);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			double temp_var = 0.0;
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				temp_var = mols[ispec].rad_trans[i].Tbr;
				binTbrfile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
				temp_var = mols[ispec].rad_trans[i].Tex;
				binTexfile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
				temp_var = mols[ispec].rad_trans[i].tau;
				bintaufile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
			}
		}
	}

	void read_output_binary_file(const string & filename) // just an example of how one can read binary files with time-dependence of level populations and/or radiation properties produced with the code
	{
		vector <double> atime;
		vector <vector <double> > data;
		vector <size_t> ntrans;
		size_t nspecies = 1;

		ifstream rfile;
		rfile.open(filename, ios::binary | ios::in);

		if (rfile.is_open()) {
			rfile.seekg(0, ios::end);
			size_t fsize = (size_t) rfile.tellg();
			rfile.seekg(0, ios::beg);

			rfile.read(reinterpret_cast<char*> (&nspecies), sizeof(size_t));
			ntrans.reserve(nspecies);
			rfile.read(reinterpret_cast<char*> (&ntrans[0]), nspecies*sizeof(ntrans[0]));

			size_t sum_of_trans = 0;
			for (size_t i = 0; i < nspecies ; i++ ) sum_of_trans += ntrans[i];
			vector <double> temp(sum_of_trans);
			double temp_time = 0;
			size_t itime = 0;
			//while (!rfile.eof()) {
			while (rfile.tellg() < fsize) {
				rfile.read(reinterpret_cast<char*> (&temp_time), sizeof(temp_time));
				atime.push_back(temp_time);
				rfile.read(reinterpret_cast<char*> (&temp[0]), sum_of_trans*sizeof(temp[0]));
				data.push_back(vector <double>());
				for (size_t itrans = 0; itrans < sum_of_trans; itrans++) data[itime].push_back(temp[itrans]);
				itime += 1;
			}
			rfile.close();
			cout << itime << endl;
		}
		else throw runtime_error("can't open binary file"); 
		//for (size_t i = 0; i < atime.size(); i++) cout << atime[i] << '\t' << data[i][1378] << endl;
	}

	void clear_mem_close_files(vector <double*> A, vector <double*> pop, vector <double*> dpop_dt, vector <vector <vector <double> > > & oldpops_Ng, vector <vector <vector <double> > > & oldpops_time, ofstream & binpopfile, ofstream & binTbrfile, ofstream & binTexfile, ofstream & bintaufile)
	{
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			delete[] A[ispec]; delete[] pop[ispec]; delete[] dpop_dt[ispec];
		}
		oldpops_Ng.clear(); oldpops_time.clear();
		if (cerr_output_iter_progress) {
			binpopfile.close(); binTbrfile.close(); binTexfile.close(); bintaufile.close();
		}
	}

	void clear_mem_close_files(vector <double*> A, vector <double*> pop, vector <double*> Jac, vector <double*> dpop_dt, vector <vector <vector <double> > > & oldpops_Ng, vector <vector <vector <double> > > & oldpops_time, ofstream & binpopfile, ofstream & binTbrfile, ofstream & binTexfile, ofstream & bintaufile)
	{
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			delete[] A[ispec]; delete[] Jac[ispec]; delete[] pop[ispec]; delete[] dpop_dt[ispec];
		}
		oldpops_Ng.clear(); oldpops_time.clear();
		if (cerr_output_iter_progress) {
			binpopfile.close(); binTbrfile.close(); binTexfile.close(); bintaufile.close();
		}
	}

	void compute_brightness_temperature(const size_t & i, const double & time, molModel *mol)		// computes brightness temperature for radiative transition i
	{ 	// it is somewhat similar to the last Equation for Tbr in Appendix A of Sobolev et al. 1997
		const double & nu = mol->rad_trans[i].nu;

		const double Tcloud_cont = exp(-mol->rad_trans[i].taud_in*beamH) * dust_HII_CMB_Jext_emission->continuum_behind_maser_region(nu, time) +
									oneMinusExp(mol->rad_trans[i].taud_in*beamH) * dust_HII_CMB_Jext_emission->inner_dust_source_function(nu, time); // contribution of the maser cloud into continuum emission
		const double Tdust_infront_cont = exp(-dust_HII_CMB_Jext_emission->tau_dust_LOS(nu, time)) * Tcloud_cont + dust_HII_CMB_Jext_emission->external_dust_layer_emission(nu, time); // absorption of continuum and emission by the external dust in front of the maser region at the line-of-sight
		const double THii_infront_cont = exp(-dust_HII_CMB_Jext_emission->tau_HII_infront(nu)) * Tdust_infront_cont; // absorption of continuum by the HII region in front of the maser region

		const double Tcloud = exp(-mol->rad_trans[i].tau*beamH) * dust_HII_CMB_Jext_emission->continuum_behind_maser_region(nu, time) +
									oneMinusExp(mol->rad_trans[i].tau*beamH) * compute_source_function(i, mol); // contribution of the maser cloud into total emission
		const double Tdust_infront = exp(-dust_HII_CMB_Jext_emission->tau_dust_LOS(nu, time)) * Tcloud + dust_HII_CMB_Jext_emission->external_dust_layer_emission(nu, time); // absorption and emission by the external dust in front of the maser region at the line-of-sight
		const double THii_infront = exp(-dust_HII_CMB_Jext_emission->tau_HII_infront(nu)) * Tdust_infront; // absorption of continuum by the HII region in front of the maser region

		mol->rad_trans[i].Tbr = (THii_infront - THii_infront_cont) * (pow(SPEED_OF_LIGHT/nu, 2.0) / (2. * BOLTZMANN_CONSTANT));
	}

	void update_external_emission(const double & time, molModel *mol)
	{
		// set the external emission mean intensity and dust emission coefficient
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			mol->rad_trans[i].JExtHII = dust_HII_CMB_Jext_emission->compute_JextHII(mol->rad_trans[i].nu, time); //external emission from HII region, should be separated from other types of emission because of maser beaming
			mol->rad_trans[i].JExt = dust_HII_CMB_Jext_emission->compute_Jext_dust_CMB_file(mol->rad_trans[i].nu, time); //external emission from dust, cosmic microwave background or file
			mol->rad_trans[i].emiss_dust = mol->rad_trans[i].kabs_dust * dust_HII_CMB_Jext_emission->inner_dust_source_function(mol->rad_trans[i].nu, time); //emission coefficient of the dust inside the maser region
		}
	}

	double populate_matrix_vector(double A[], double B[], double Jac[], beta_LVG & LVG_beta, const vector <vector <double> > & oldpops, const double BDF_coeffs[], const double & h, double dpop_dt[], molModel *mol)	// fill the matrix A and vector B from the non-linear system of equations A*X=B which is obtained at each step of integration over time and Jacobian used to solve the system; compute vector of derivatives of populations over time and its norm
	{
		const size_t & n = mol->levels.size();
		
		// A[i * j*n] - rate of transition from j-th to i-th level
		// Collisional transitions
		// Note that diagonal elements of C were computed in compute_C function in molModel.h, Cii = sum{k=1,Nlevel}(Cik)
		for (size_t i = 0; i < n; i++) {
			A[i + i*n] = - mol->coll_trans[i][i];
			Jac[i + i*n] = A[i + i*n];
			for (size_t j = i+1; j < n; j++) {
				A[i + j*n] = mol->coll_trans[j][i];
				A[j + i*n] = mol->coll_trans[i][j];
				Jac[i + j*n] = A[i + j*n];
				Jac[j + i*n] = A[j + i*n];
			}
		}
		
        // Radiative transitions
		double temp_var, temp_var_Jac, Sf, beta, Jin, kabs;
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			compute_tau(i, mol);
			compute_J_S_beta(mol, i, LVG_beta, Sf, beta, Jin, kabs);
            
			const size_t & up = mol->rad_trans[i].up_level;
			const size_t & low = mol->rad_trans[i].low_level;
	
			temp_var = mol->rad_trans[i].Blu * mol->rad_trans[i].J;
			A[up + low*n]  += temp_var;                                   	// A[up][low] += Blu * J
			A[low + low*n] -= temp_var;    									// A[low][low] -= Blu * J

			double blend_demiss_dnlow = 0.0;								// derivative of emission coefficient on population of the lower level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_demiss_dnlow += (int)(low == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].A * mol->rad_trans[i].blends[j].fac;
			}
			double blend_dkabs_dnlow = mol->rad_trans[i].Blu;				// derivative of absorption coefficient on population of the lower level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_dkabs_dnlow += (int)(low == mol->rad_trans[mol->rad_trans[i].blends[j].id].low_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Blu * mol->rad_trans[i].blends[j].fac;
				blend_dkabs_dnlow -= (int)(low == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level  && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Bul * mol->rad_trans[i].blends[j].fac;
			}
			double dS_dnlow_beta = HC4PI * invlineWidth * modelPhysPars::n_mol[mol->idspec] * (blend_demiss_dnlow - Sf * blend_dkabs_dnlow);	// derivative of source function on population of the lower level of radiative transition * (1 - beta)
			if (fabs(kabs) > 0.0) dS_dnlow_beta *= (1. - beta) / kabs;
			else dS_dnlow_beta *= 0.5 * modelPhysPars::NdV[mol->idspec] * lineWidth / (modelPhysPars::n_mol[mol->idspec]); // (1-b)/tau -> 0.5 for tau -> 0.0
			double dtau_dnlow = HC4PI * modelPhysPars::NdV[mol->idspec] * blend_dkabs_dnlow;													// derivative of optical depth on population of the lower level of radiative transition
			temp_var_Jac = mol->levels[low].pop * mol->rad_trans[i].Blu * (
					dS_dnlow_beta + dtau_dnlow * (
					(mol->rad_trans[i].JExt - Sf) * LVG_beta.DbetaDtau(mol->rad_trans[i].tau) +
					dust_HII_CMB_Jext_emission->HII_region_at_LOS * LVG_beta.DbetaHIIDtau_LOS(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII +
					(1 - dust_HII_CMB_Jext_emission->HII_region_at_LOS) * LVG_beta.DbetaHIIDtau_pump(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII
				)
			);
			Jac[up + low*n] += temp_var + temp_var_Jac;
			Jac[low + low*n] -= (temp_var + temp_var_Jac);

			temp_var = (mol->rad_trans[i].A + mol->rad_trans[i].Bul * mol->rad_trans[i].J);
			A[low + up*n]  += temp_var;                                   	// A[low][up] += Aul + Bul * J
			A[up + up*n] -= temp_var;       								// A[up][up] -= Aul + Bul * J , where Aul, Bul, Blu - Einstein coefficients, J - mean intensity

			double blend_demiss_dnup = mol->rad_trans[i].A;								// derivative of emission coefficient on population of the upper level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_demiss_dnup += (int)(up == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].A * mol->rad_trans[i].blends[j].fac;
			}
			double blend_dkabs_dnup = - mol->rad_trans[i].Bul;				// derivative of absorption coefficient on population of the upper level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_dkabs_dnup += (int)(up == mol->rad_trans[mol->rad_trans[i].blends[j].id].low_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Blu * mol->rad_trans[i].blends[j].fac;
				blend_dkabs_dnup -= (int)(up == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level  && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Bul * mol->rad_trans[i].blends[j].fac;
			}
			double dS_dnup_beta = HC4PI * invlineWidth * modelPhysPars::n_mol[mol->idspec] * (blend_demiss_dnup - Sf * blend_dkabs_dnup);	// derivative of source function on population of the upper level of radiative transition * (1 - beta)
			if (fabs(kabs) > 0.0) dS_dnup_beta *= (1. - beta) / kabs;
			else dS_dnup_beta *= 0.5 * modelPhysPars::NdV[mol->idspec] * lineWidth / (modelPhysPars::n_mol[mol->idspec]); // (1-b)/tau -> 0.5 for tau -> 0.0
			double dtau_dnup = HC4PI * modelPhysPars::NdV[mol->idspec] * blend_dkabs_dnup;													// derivative of optical depth on population of the upper level of radiative transition
			temp_var_Jac = mol->levels[up].pop * mol->rad_trans[i].Bul * (
					dS_dnup_beta + dtau_dnup * (
					(mol->rad_trans[i].JExt - Sf) * LVG_beta.DbetaDtau(mol->rad_trans[i].tau) +
					dust_HII_CMB_Jext_emission->HII_region_at_LOS * LVG_beta.DbetaHIIDtau_LOS(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII +
					(1 - dust_HII_CMB_Jext_emission->HII_region_at_LOS) * LVG_beta.DbetaHIIDtau_pump(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII
				)
			);
			Jac[low + up*n] += temp_var + temp_var_Jac;
			Jac[up + up*n] -= (temp_var + temp_var_Jac);
		}
		
		double norm = 0.0;
		//double pops_sum = 0.0;
		for (size_t i = n; i-- > 0; ) {
			dpop_dt[i] = 0.0;
			for (size_t j = n; j-- > 0; ) {
				dpop_dt[i] += A[i + j*n] * mol->levels[j].pop;
				A[i + j*n] = Jac[i + j*n];
			}
			norm += dpop_dt[i] * dpop_dt[i];
			B[i] = (BDF_coeffs[0] * mol->levels[i].pop + BDF_coeffs[1] * oldpops[i][0] + BDF_coeffs[2] * oldpops[i][1]) / h - dpop_dt[i];		// BDF method formula
			A[i + i*n] -= BDF_coeffs[0] / h;
			//A[0 + i*n] = 1.0;	// A[0][i] = 1.0 - the equation for the first level is replaced by the particle number conservation law, i.e. the sum of populations should be = partition functions ratio or = 1
			//pops_sum += mol->levels[i].pop;
		}
		//B[0] = this->partition_function_ratio[mol->idspec] - pops_sum; // the sum of populations should be = partition functions ratio or = 1
		return sqrt(norm);
	}

public:
	
	int radiative_transfer() override		// main function which should be called to perform radiative transfer calculations
	{
		initial_solution();							// get the initial values of populations

		double time = 0.0;

		// prepare to output the dependence of level populations on time into the binary file
		ofstream binpopfile, binTbrfile, binTexfile, bintaufile;
		if (cerr_output_iter_progress) {
			const size_t nspec = modelPhysPars::nSpecies;
			binpopfile.open("pops_vs_time.bin", ios::out | ios::binary);
			binpopfile.write(reinterpret_cast<const char*>(&nspec), sizeof(size_t));
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				const size_t nlevs = mols[ispec].levels.size();
				binpopfile.write(reinterpret_cast<const char*>(&nlevs), sizeof(size_t));
			}
			binTbrfile.open("Tbr_vs_time.bin", ios::out | ios::binary);
			binTexfile.open("Tex_vs_time.bin", ios::out | ios::binary);
			bintaufile.open("tau_vs_time.bin", ios::out | ios::binary);
			binTbrfile.write(reinterpret_cast<const char*>(&nspec), sizeof(size_t));
			binTexfile.write(reinterpret_cast<const char*>(&nspec), sizeof(size_t));
			bintaufile.write(reinterpret_cast<const char*>(&nspec), sizeof(size_t));
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				const size_t ntrans = mols[ispec].rad_trans.size();
				binTbrfile.write(reinterpret_cast<const char*>(&ntrans), sizeof(size_t));
				binTexfile.write(reinterpret_cast<const char*>(&ntrans), sizeof(size_t));
				bintaufile.write(reinterpret_cast<const char*>(&ntrans), sizeof(size_t));
			}
		}

		// set the external emission mean intensity, optical depth, absorption and emission coefficients
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				mols[ispec].rad_trans[i].taud_in = dust_HII_CMB_Jext_emission->tau_dust_in(mols[ispec].rad_trans[i].nu, lineWidth); //optical depth of the dust inside the maser region
				mols[ispec].rad_trans[i].kabs_dust = mols[ispec].rad_trans[i].taud_in * modelPhysPars::H2dens * invlineWidth / modelPhysPars::max_NH2dV; //absorption coefficient of the dust inside the maser region
				// mols[ispec].rad_trans[i].JExtHII will be zero if external emission will be taken from file
			}
		}
		
		if (lineWidth > DBL_EPSILON) find_blends(); // find overlapping lines

		beta_LVG LVG_beta (beamH); 					// LVG escape probability, see beta_LVG.h

		vector <vector <vector <double> > > oldpops_Ng;		// stores populations used for Ng acceleration
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			oldpops_Ng.push_back(vector <vector <double> >());
			for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
				oldpops_Ng[ispec].push_back(vector <double>());
				for (size_t j = 0; j < (Ng_order + 2); j++) {
					oldpops_Ng[ispec][i].push_back(mols[ispec].levels[i].pop);
				}
			}
		}

		vector <vector <vector <double> > > oldpops_time;		// stores populations at previous time steps
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			oldpops_time.push_back(vector <vector <double> >());
			for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
				oldpops_time[ispec].push_back(vector <double>());
				for (size_t j = 0; j < 2; j++) {
					oldpops_time[ispec][i].push_back(mols[ispec].levels[i].pop);
				}
			}
		}

		if (cerr_output_iter_progress) {
			binpopfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
					double temp_var = oldpops_time[ispec][i][0];
					binpopfile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
				}
			}
			write_rad_trans_data_into_binary_file(LVG_beta, time, binTbrfile, binTexfile, bintaufile);
		}

		vector <double*> pop;
		pop.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			pop.push_back(new double[mols[ispec].levels.size()]);
		}

		vector <double*> A; 	// reserve space for matrix A from the statistical equilibrium equations system A*X=B
		A.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			A.push_back(new double[mols[ispec].levels.size()*mols[ispec].levels.size()]);
		}

		vector <double*> Jac; 	// reserve space for Jacobian matrix
		Jac.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			Jac.push_back(new double[mols[ispec].levels.size()*mols[ispec].levels.size()]);
		}

		vector <double*> dpop_dt;
		dpop_dt.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			dpop_dt.push_back(new double[mols[ispec].levels.size()]);
		}

		double h = INITIAL_TIME_STEP, h_old = h;	// INITIAL_TIME_STEP is from hiddenParameters.h
		unsigned long ntimesteps = 0;
		double BDF_coeffs[3] = {1.0, -1.0, 0.0};	// the coefficients from BDF method formula; the first time step corresponds to implicit Euler method
		
		double F_norm = 0.0;
		double MaxRPopDiff = 0.0;
		size_t levelWithMaxRPopDiff = 0;
		size_t speciesWithMaxRPopDiff = 0;
		// integrate the kinetic equations system with the second order BDF method using variable time step (see Jannelli & Fazio 2006, Journal of Computational and Applied Mathematics, 191, 246)
		do {
			unsigned int iter_in = 0;
			double rn = 1.0;
			bool there_were_bad_levels = false;
			bool solution_failed = false;
			vector <double> pop_norm(modelPhysPars::nSpecies);
			size_t solveStatEqSuccess = 0;
			double MaxRPopDiff_at_prev_iter = 0.0;
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) update_external_emission(time + h, &mols[ispec]);
			if (USE_PREDICTOR && ntimesteps > 0) {
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					populate_matrix_vector(A[ispec], pop[ispec], Jac[ispec], LVG_beta, oldpops_time[ispec], BDF_coeffs, h, dpop_dt[ispec], &mols[ispec]);
				}
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					for (size_t i = 0; i < mols[ispec].levels.size(); i++) mols[ispec].levels[i].pop = fma(dpop_dt[ispec][i], h, mols[ispec].levels[i].pop);
				}
			}
			do {							// solve non-linear system of equations at a given time step with Newton method
				solveStatEqSuccess = 0;
				F_norm = 0.0;
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					double ispecF_norm = populate_matrix_vector(A[ispec], pop[ispec], Jac[ispec], LVG_beta, oldpops_time[ispec], BDF_coeffs, h, dpop_dt[ispec], &mols[ispec]);
					if (ispecF_norm > F_norm) F_norm = ispecF_norm;
					solveStatEqSuccess += solve_eq_sys(A[ispec], pop[ispec], &mols[ispec]);	// solve linear system of equations with LU decomposition, solution is stored in pop
				}
				MaxRPopDiff = 0.0;
				levelWithMaxRPopDiff = 0;
				speciesWithMaxRPopDiff = 0;
				there_were_bad_levels = false;
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					double pops_sum = 0.0;
					for (size_t i = mols[ispec].levels.size(); i-- > 1; ) {
						pop[ispec][i] = mols[ispec].levels[i].pop + pop[ispec][i];
						pops_sum += pop[ispec][i];
					}
					pop[ispec][0] = this->partition_function_ratio[ispec] - pops_sum;
					if (solveStatEqSuccess == 0) {
						there_were_bad_levels = ( update_check_pops(&mols[ispec], pop[ispec], speciesWithMaxRPopDiff, levelWithMaxRPopDiff, MaxRPopDiff, iter_in, oldpops_Ng[ispec], pop_norm[ispec]) || there_were_bad_levels);
					}
				}
				if (solveStatEqSuccess == 0 && cerr_output_iter_progress) {
					cerr << iter_in << " max.dev.= " << MaxRPopDiff << " mol/level with max.dev.= " << speciesWithMaxRPopDiff << " / " << levelWithMaxRPopDiff << endl;
				}
				if ( there_were_bad_levels || solveStatEqSuccess != 0 || fabs(MaxRPopDiff_at_prev_iter - MaxRPopDiff) / MaxRPopDiff < MAX_LOCAL_ACCURACY) {
					solution_failed = true;
					break;
				}
				iter_in += 1;
				MaxRPopDiff_at_prev_iter = MaxRPopDiff;
			} while (MaxRPopDiff > MAX_LOCAL_ACCURACY && iter_in < MAX_NUM_INNER_STEPS);
			
			double max_monitor_function = -1.e30;
			if (!solution_failed && iter_in < MAX_NUM_INNER_STEPS) {
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					double monitor_function = 0.e0;
					for (size_t i = mols[ispec].levels.size(); i-- > 0; ) {
						monitor_function += (oldpops_time[ispec][i][0] - mols[ispec].levels[i].pop) * (oldpops_time[ispec][i][0] - mols[ispec].levels[i].pop);
					}
					monitor_function = sqrt(monitor_function) / pop_norm[ispec];
					if (monitor_function > max_monitor_function) max_monitor_function = monitor_function;
				}
			}
			if (max_monitor_function > MON_FUN_UPPER_LIMIT || max_monitor_function < 0) {	// one need to roll-back populations and repeat integration with a smaller time step
				h *= TIME_STEP_DECREASE_FAC;
				if (h >= MIN_TIME_STEP) {
					if (solveStatEqSuccess == 0) {
						for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
							for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
								mols[ispec].levels[i].pop = oldpops_time[ispec][i][0];
							}
						}
					}
					if (ntimesteps > 0) {
						rn = h / h_old;
						BDF_coeffs[0] = (1+2*rn)/(1+rn); BDF_coeffs[1] = -(1+rn); BDF_coeffs[2] = rn*rn/(1+rn);
					}
				} else {
					clear_mem_close_files(A, pop, Jac, dpop_dt, oldpops_Ng, oldpops_time, binpopfile, binTbrfile, binTexfile, bintaufile);
					cerr << "#error: cant decrease timestep more, info = " << there_were_bad_levels << "," << solveStatEqSuccess << endl;
					return 1;
				}
			} else {				// proceed the integration over time
				ntimesteps += 1;
				time += h;
				h_old = h;
				if (cerr_output_iter_progress) {
					write_rad_trans_data_into_binary_file(LVG_beta, time, binTbrfile, binTexfile, bintaufile);
					cerr << "tint= " << ntimesteps << " time= " << time << " h= " << h << " n= " << F_norm << " max.dev.= " << MaxRPopDiff << " mol/level with max.dev.= " << speciesWithMaxRPopDiff << " / " << levelWithMaxRPopDiff << endl;
					binpopfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
				}
				if (max_monitor_function < MON_FUN_LOWER_LIMIT) {
					h = min(MAX_TIME_STEP, h * TIME_STEP_INCREASE_FAC);
				}
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
						oldpops_time[ispec][i][1] = oldpops_time[ispec][i][0];
						oldpops_time[ispec][i][0] = mols[ispec].levels[i].pop;
						if (cerr_output_iter_progress) {
							double temp_var = oldpops_time[ispec][i][0];
							binpopfile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
						}
					}
				}
				rn = h / h_old;
				BDF_coeffs[0] = (1+2*rn)/(1+rn); BDF_coeffs[1] = -(1+rn); BDF_coeffs[2] = rn*rn/(1+rn);
			}
		} while (F_norm > MAX_DpopsDt_EPS && ntimesteps < maxNumberOfIterations);

		clear_mem_close_files(A, pop, Jac, dpop_dt, oldpops_Ng, oldpops_time, binpopfile, binTbrfile, binTexfile, bintaufile);

		if (ntimesteps == maxNumberOfIterations) cerr << "#warning: maximum number of time steps has been exceeded " << "n= " << F_norm << " max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;
		
		prepare_results_for_output(LVG_beta, time);
		
		return 0;
	}

	RT_time_Newt() : RT()
	{}

	RT_time_Newt(istream & cin) : RT(cin)
	{}
};
