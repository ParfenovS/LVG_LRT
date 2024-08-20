#pragma once
#include "RT.h"

class RT_time_Newt : public RT		// integrates kinetic equations for level populations over time; the non-linear system of equations at each time step is solved with simple iterations
{
private:

	void prepare_results_for_output(beta_LVG & LVG_beta, const double & time)
	{
		// preparing quantities for output
		double dummy_S = 0.0;
		double dummy_beta = 0.0;
		double dummy_kabs = 0.0;
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				compute_tau(i, &mols[ispec]); 		// computing final optical depths that can be used for output
				compute_J_S_beta(&mols[ispec], i, LVG_beta, dummy_S, dummy_beta, dummy_kabs); 	// computing final mean intensities that can be used for output
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

	void close_files(ofstream & binpopfile, ofstream & binTbrfile, ofstream & binTexfile, ofstream & bintaufile)
	{
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

	pair <double, double> populate_matrix_vector(vector <double> & A, vector <double> & B, beta_LVG & LVG_beta, const vector <vector <double> > & oldpops, const vector <double> & BDF_coeffs, const double & h, vector <double> & dpop_dt, vector <double> & dpop_dt_old, molModel *mol, const size_t & stage, const size_t & maxPopid)	// fill the matrix A and vector B from the non-linear system of equations A*X=B which is obtained at each step of integration over time and Jacobian used to solve the system; compute vector of derivatives of populations over time and its norm
	{
		get_A(A, LVG_beta, mol);

		const size_t & n = mol->levels.size();
		
		double norm = 0.0;
		double norm_local = 0.0;
		double pops_sum = 0.0;
		dpop_dt[maxPopid] = 0.0;
		for (size_t j = n; j-- > 0; ) {
			dpop_dt[maxPopid] += A[maxPopid + j*n] * mol->levels[j].pop;
		}
		if (stage == 0) dpop_dt_old[maxPopid] = dpop_dt[maxPopid];
		for (size_t i = n; i-- > 0; ) {
			if (i != maxPopid) {
				dpop_dt[i] = 0.0;
				for (size_t j = n; j-- > 0; ) {
					dpop_dt[i] += A[i + j*n] * mol->levels[j].pop;
				}
				if (stage == 0) dpop_dt_old[i] = dpop_dt[i];
				B[i] = (BDF_coeffs[1] * oldpops[i][0] + BDF_coeffs[2] * oldpops[i][1]) / h + BDF_coeffs[3] * dpop_dt_old[i];		// TRBDF2 method formula	
				A[i + i*n] -= BDF_coeffs[0] / h;
				norm_local += (dpop_dt[i] - B[i]) * (dpop_dt[i] - B[i]);
			}
			norm += dpop_dt[i] * dpop_dt[i];
			A[maxPopid + i*n] = 1.0;	// A[maxPopid][i] = 1.0 - the equation for the level with maximum initial population is replaced by the particle number conservation law, i.e. the sum of populations should be = partition functions ratio or = 1
			pops_sum += mol->levels[i].pop;
		}
		B[maxPopid] = this->partition_function_ratio[mol->idspec]; // the sum of populations should be = partition functions ratio or = 1
		norm_local += (B[maxPopid] - pops_sum) * (B[maxPopid] - pops_sum);
		return make_pair(sqrt(norm), sqrt(norm_local));
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

		vector <vector <vector <double> > > backup_pops;		// stores populations at previous time steps
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			backup_pops.push_back(vector <vector <double> >());
			for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
				backup_pops[ispec].push_back(vector <double>());
				for (size_t j = 0; j < 2; j++) {
					backup_pops[ispec][i].push_back(mols[ispec].levels[i].pop);
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

		vector <vector <double> > pop(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			pop[ispec].resize(mols[ispec].levels.size());
		}

		vector <vector <double> > A(modelPhysPars::nSpecies); // reserve space for matrix A from the kinetic equations system A*X=B
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			A[ispec].resize(mols[ispec].levels.size() * mols[ispec].levels.size());
		}

		vector <vector <double> > dpop_dt(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			dpop_dt[ispec].resize(mols[ispec].levels.size());
		}

		vector <vector <double> > dpop_dt_old(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			dpop_dt_old[ispec].resize(mols[ispec].levels.size());
		}

		vector <size_t> max_pop_id(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {;
			max_pop_id[ispec] = 0;
			double max_pop = mols[ispec].levels[0].pop;
			for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
				if (mols[ispec].levels[i].pop > max_pop) {
					max_pop = mols[ispec].levels[i].pop;
					max_pop_id[ispec] = i;
				}
			}
		}

		double h = INITIAL_TIME_STEP;	// INITIAL_TIME_STEP is from hiddenParameters.h
		unsigned long ntimesteps = 0;
		vector <double> BDF_coeffs = {2.0 / TRBDF2_GAMMA, -2.0 / TRBDF2_GAMMA, 0.0, -1.0};	// the coefficients from TRBDF2 method formula
		
		double F_norm = 0.0;
		double F_norm_local = 0.0;
		double MaxRPopDiff = 0.0;
		size_t levelWithMaxRPopDiff = 0;
		size_t speciesWithMaxRPopDiff = 0;
		// integrate the kinetic equations system with the second order BDF method using variable time step (see Jannelli & Fazio 2006, Journal of Computational and Applied Mathematics, 191, 246)
		do {
			unsigned int iter_in = 1;
			bool there_were_bad_levels = false;
			bool solution_failed = false;
			vector <double> pop_norm(modelPhysPars::nSpecies);
			size_t solveStatEqSuccess = 0;
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) update_external_emission(time + h, &mols[ispec]);
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				populate_matrix_vector(A[ispec], pop[ispec], LVG_beta, oldpops_time[ispec], BDF_coeffs, h, dpop_dt[ispec], dpop_dt_old[ispec], &mols[ispec], 0, max_pop_id[ispec]);
			}
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				for (size_t i = 0; i < mols[ispec].levels.size(); i++) mols[ispec].levels[i].pop = fma(dpop_dt[ispec][i], h, mols[ispec].levels[i].pop);
			}

			for (size_t trbdf_stage = 1; trbdf_stage <= 2; trbdf_stage++) {
				if (solution_failed) break;
				if (trbdf_stage == 1) BDF_coeffs = {2.0 / TRBDF2_GAMMA, -2.0 / TRBDF2_GAMMA, 0.0, -1.0};
				else BDF_coeffs = {(2.0 - TRBDF2_GAMMA) / (1.0 - TRBDF2_GAMMA), -1.0 / (TRBDF2_GAMMA * (1.0 - TRBDF2_GAMMA)), (1.0 - TRBDF2_GAMMA) / TRBDF2_GAMMA, 0.0};
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					max_pop_id[ispec] = 0;
					double max_pop = mols[ispec].levels[0].pop;
					for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
						if (mols[ispec].levels[i].pop > max_pop) {
							max_pop = mols[ispec].levels[i].pop;
							max_pop_id[ispec] = i;
						}
					}
				}
				do {							// solve non-linear system of equations at a given time step with simple iterations
					solveStatEqSuccess = 0;
					F_norm = 0.0;
					F_norm_local = 0.0;
					for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
						auto [ispecF_norm, ispecF_norm_local] = populate_matrix_vector(A[ispec], pop[ispec], LVG_beta, oldpops_time[ispec], BDF_coeffs, h, dpop_dt[ispec], dpop_dt_old[ispec], &mols[ispec], trbdf_stage, max_pop_id[ispec]);
						if (ispecF_norm > F_norm) F_norm = ispecF_norm;
						if (ispecF_norm_local > F_norm_local) F_norm_local = ispecF_norm_local;
						solveStatEqSuccess += solve_eq_sys(A[ispec].data(), pop[ispec].data(), &mols[ispec]);	// solve linear system of equations with LU decomposition, solution is stored in pop
					}
					MaxRPopDiff = 0.0;
					levelWithMaxRPopDiff = 0;
					speciesWithMaxRPopDiff = 0;
					there_were_bad_levels = false;
					if (solveStatEqSuccess == 0) {
						for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
							there_were_bad_levels = (update_check_pops(&mols[ispec], pop[ispec], speciesWithMaxRPopDiff, levelWithMaxRPopDiff, MaxRPopDiff, iter_in, oldpops_Ng[ispec], pop_norm[ispec]) || there_were_bad_levels);
						}
						if (cerr_output_iter_progress) {
							cerr << iter_in << " Fnorm= " << F_norm_local << " max.dev.= " << MaxRPopDiff << " mol/level with max.dev.= " << speciesWithMaxRPopDiff << " / " <<levelWithMaxRPopDiff << endl;
						}
					}
					if (there_were_bad_levels || solveStatEqSuccess != 0 || MaxRPopDiff == 0.0) {
						solution_failed = true;
						break;
					}
					iter_in += 1;
				} while (MaxRPopDiff > MAX_LOCAL_ACCURACY && iter_in <= MAX_NUM_INNER_STEPS);
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
						oldpops_time[ispec][i][1] = oldpops_time[ispec][i][0];
						oldpops_time[ispec][i][0] = mols[ispec].levels[i].pop;
					}
				}
			}
			
			double max_monitor_function = -1.e30;
			if (!solution_failed && iter_in <= MAX_NUM_INNER_STEPS) {
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
								mols[ispec].levels[i].pop = backup_pops[ispec][i][0];
								oldpops_time[ispec][i][0] = backup_pops[ispec][i][0];
								oldpops_time[ispec][i][1] = backup_pops[ispec][i][1];
							}
						}
					}
				} else {
					close_files(binpopfile, binTbrfile, binTexfile, bintaufile);
					cerr << "#error: cant decrease timestep more, info = " << there_were_bad_levels << "," << solveStatEqSuccess << endl;
					return 1;
				}
			} else {				// proceed the integration over time
				ntimesteps += 1;
				time += h;
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
						backup_pops[ispec][i][1] = backup_pops[ispec][i][0];
						backup_pops[ispec][i][0] = mols[ispec].levels[i].pop;
						if (cerr_output_iter_progress) {
							double temp_var = oldpops_time[ispec][i][0];
							binpopfile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
						}
					}
				}
			}
		} while (F_norm > MAX_DpopsDt_EPS && ntimesteps < maxNumberOfIterations);

		close_files(binpopfile, binTbrfile, binTexfile, bintaufile);

		if (ntimesteps == maxNumberOfIterations) cerr << "#warning: maximum number of time steps has been exceeded " << "n= " << F_norm << " max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;
		
		prepare_results_for_output(LVG_beta, time);
		
		return 0;
	}

	RT_time_Newt() : RT()
	{}

	RT_time_Newt(istream & cin) : RT(cin)
	{}
};
