#pragma once
#include "RT.h"

class RT_time_Newt : public RT		// integrates kinetic equations for level populations over time; the non-linear system of equations at each time step is solved with Newton method
{
private:

	double populate_matrix_vector(double A[], double B[], double Jac[], beta_LVG & LVG_beta, const vector <vector <double> > & oldpops, const double BDF_coeffs[], const double & h, double dpop_dt[])	// fill the matrix A and vector B from the non-linear system of equations A*X=B which is obtained at each step of integration over time and Jacobian used to solve the system; compute vector of derivatives of populations over time and its norm
	{
		const size_t & n = mol->levels.size();
		
		// Collisional transitions
		// Note that diagonal elements of C were computed in compute_C function in molModel.h, Cii = - sum{k=1,Nlevel}(Cik)
		for (size_t i = 0; i < n; i++) {
			A[i + i*n] = mol->coll_trans[i][i];
			Jac[i + i*n] = A[i + i*n];
			for (size_t j = i+1; j < n; j++) {
				A[i + j*n] = mol->coll_trans[j][i];
				A[j + i*n] = mol->coll_trans[i][j];
				Jac[i + j*n] = A[i + j*n];
				Jac[j + i*n] = A[j + i*n];
			}
		}
		
        // Radiative transitions
		double temp_var, Sf, beta, Jin;
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			compute_tau(i);
			compute_J_S_beta(i, LVG_beta, Sf, beta, Jin);
            
			const size_t & up = mol->rad_trans[i].up_level;
			const size_t & low = mol->rad_trans[i].low_level;
	
			temp_var = mol->rad_trans[i].Blu * mol->rad_trans[i].J;
			A[up + low*n]  += temp_var;                                   	// A[up][low] += Blu * J
			A[low + low*n] -= temp_var;    									// A[low][low] -= Blu * J

			temp_var = (mol->rad_trans[i].A + mol->rad_trans[i].Bul * mol->rad_trans[i].J);
			A[low + up*n]  += temp_var;                                   	// A[low][up] += Aul + Bul * J
			A[up + up*n] -= temp_var;       								// A[up][up] -= Aul + Bul * J , where Aul, Bul, Blu - Einstein coefficients, J - mean intensity

			temp_var = LVG_beta.tauDbetaDtau(mol->rad_trans[i].tau) * (Sf - mol->rad_trans[i].JExt - mol->rad_trans[i].JExtHII);

			Jac[up + low*n]  = A[up + low*n]  - mol->rad_trans[i].Blu * (Jin + temp_var);		// Jac[up][low] -= Blu * J + Blu*tau*DbetaDtau * (S-Jext)
			Jac[low + up*n]  = A[low + up*n]  - (mol->rad_trans[i].A * (1.- beta) + mol->rad_trans[i].Bul * (Jin + temp_var)); // Jac[low][up] -= (1-beta)*Aul + Bul * J + Bul*tau*DbetaDtau * (S-Jext)
			Jac[up + up*n]   = A[up + up*n]   + (mol->rad_trans[i].A * (1.- beta) + mol->rad_trans[i].Bul * (Jin + temp_var)); // Jac[up][up] += (1-beta)*Aul + Bul * J + Bul*tau*DbetaDtau * (S-Jext)
			Jac[low + low*n] = A[low + low*n] + mol->rad_trans[i].Blu * (Jin + temp_var); // Jac[low][low] += Blu * J + Blu*tau*DbetaDtau * (S-Jext)
		}

		double norm = 0.0;
		for (size_t i = 0; i < n; i++) {
			dpop_dt[i] = 0.0;
			for (size_t j = 0; j < n; j++) {
				dpop_dt[i] += A[i + j*n] * mol->levels[j].pop;
				A[i + j*n] = h * Jac[i + j*n];
			}
			norm += dpop_dt[i] * dpop_dt[i];
			B[i] = BDF_coeffs[0] * mol->levels[i].pop + BDF_coeffs[1] * oldpops[i][0] + BDF_coeffs[2] * oldpops[i][1] - h * dpop_dt[i];		// BDF method formula
			A[i + i*n] -= BDF_coeffs[0];
			//A[0 + i*n] = 1.0;	// A[0][i] = 1.0 - the equation for the first level is replaced by the particle number conservation law, i.e. the sum of populations should be = partition functions ratio or = 1
		}
		//B[0] = this->partition_function_ratio; // the sum of populations should be = partition functions ratio or = 1
		return sqrt(norm);
	}

public:
	
	int radiative_transfer() override		// main function which should be called to perform radiative transfer calculations
	{
		const size_t nlevs = mol->levels.size();

		// prepare to output the dependence of level populations on time into the binary file
		ofstream binpopfile;
		if (cerr_output_iter_progress) {
			binpopfile.open("pops_vs_time_NR.bin", ios::out | ios::binary);
			binpopfile.write(reinterpret_cast<const char*>(&nlevs), sizeof(size_t));
		}

		// set the external emission mean intensity, optical depth, absorption and emission coefficients
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			mol->rad_trans[i].JExt = dust_HII_CMB_Jext_emission->compute_Jext_dust_CMB_file(mol->rad_trans[i].nu); //external emission from dust, cosmic microwave background or file
			mol->rad_trans[i].JExtHII = dust_HII_CMB_Jext_emission->compute_JextHII(mol->rad_trans[i].nu); //external emission from HII region, should be separated from other types of emission because of maser beaming
			mol->rad_trans[i].taud_in = dust_HII_CMB_Jext_emission->tau_dust_in(mol->rad_trans[i].nu, lineWidth); //optical depth of the dust inside the maser region
			mol->rad_trans[i].kabs_dust = mol->rad_trans[i].taud_in * (4.*PI/SPEED_OF_LIGHT/PLANK_CONSTANT) / modelPhysPars::NdV; //absorption coefficient of the dust inside the maser region
			mol->rad_trans[i].emiss_dust = mol->rad_trans[i].kabs_dust * dust_HII_CMB_Jext_emission->inner_dust_source_function(mol->rad_trans[i].nu); //absorption coefficient of the dust inside the maser region
			// mol->rad_trans[i].JExtHII will be zero if external emission will be taken from file
		}
		
		initial_solution();							// get the initial values of populations

		if (lineWidth > DBL_EPSILON) find_blends(); // find overlapping lines

		beta_LVG LVG_beta (beamH); 					// LVG escape probability, see beta_LVG.h
		
		vector <vector <double> > oldpops_Ng;		// stores populations used for Ng acceleration
		vector <vector <double> > oldpops_time;
		double time = 0.0;
		if (cerr_output_iter_progress) binpopfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
		for (size_t i = 0; i < nlevs; i++) {
			oldpops_Ng.push_back(vector <double>());
			oldpops_time.push_back(vector <double>());
			for (size_t j = 0; j < (Ng_order + 2); j++) {
				oldpops_Ng[i].push_back(mol->levels[i].pop);
			}
			for (size_t j = 0; j < 2; j++) {
				oldpops_time[i].push_back(mol->levels[i].pop);
			}
			if (cerr_output_iter_progress) {
				double temp_pop = oldpops_time[i][0];
				binpopfile.write(reinterpret_cast<const char*>(&temp_pop), sizeof(double));
			}
		}

		double MaxRPopDiff = 0.0;
		double F_norm;
		size_t levelWithMaxRPopDiff = 0;

		double *A = new double[nlevs*nlevs]; 	// reserve space for matrix from the non-linear system of equations A*X=B obtained at each step of integration over time
		double *Jac = new double[nlevs*nlevs]; 	// reserve space for Jacobian matrix
		double *pop = new double[nlevs];
		double *dpop_dt = new double[nlevs];

		double h = INITIAL_TIME_STEP, h_old = h;	// INITIAL_TIME_STEP is from hiddenParameters.h
		unsigned long ntimesteps = 0;
		double BDF_coeffs[3] = {1.0, -1.0, 0.0};	// the coefficients from BDF method formula; the first time step corresponds to implicit Euler method
		
		// integrate the kinetic equations system with the second order BDF method using variable time step (see Jannelli & Fazio 2006, Journal of Computational and Applied Mathematics, 191, 246)
		do {
			unsigned int iter_in = 0;
			double rn = 1.0;
			bool there_were_bad_levels = false;
			bool solution_failed = false;
			double pop_norm = 1.e-30;
			size_t solveStatEqSuccess = 0;
			do {							// solve non-linear system of equations at a given time step with Newton method
				F_norm = populate_matrix_vector(A, pop, Jac, LVG_beta, oldpops_time, BDF_coeffs, h, dpop_dt);
				solveStatEqSuccess = solve_eq_sys(A, pop);	// solve linear system of equations with LU decomposition, solution is stored in pop
				double pops_sum = 0.0;
				for (size_t i = 1; i < nlevs; i++) {
					pop[i] = mol->levels[i].pop + pop[i];
					pops_sum += pop[i];
				}
				pop[0] = this->partition_function_ratio - pops_sum;
				if (solveStatEqSuccess == 0) {
					there_were_bad_levels = update_check_pops(pop, levelWithMaxRPopDiff, MaxRPopDiff, iter_in, oldpops_Ng, pop_norm, 1.0);
					if (cerr_output_iter_progress) {
						cerr << iter_in << " max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;
					}
				}
				if ( there_were_bad_levels || solveStatEqSuccess != 0 ) {
					solution_failed = true;
					break;
				}
				iter_in += 1;
			} while (MaxRPopDiff > MAX_LOCAL_ACCURACY && iter_in < MAX_NUM_INNER_STEPS);
			double monitor_function = 1.e30;
			if (!solution_failed && iter_in < MAX_NUM_INNER_STEPS) {
				monitor_function = 0.e0;
				for (size_t i = 0; i < nlevs; i++) {
					monitor_function += (oldpops_time[i][0] - mol->levels[i].pop) * (oldpops_time[i][0] - mol->levels[i].pop);
				}
				monitor_function = sqrt(monitor_function) / sqrt(pop_norm);
			}
			if (monitor_function > MON_FUN_UPPER_LIMIT) {	// one need to roll-back populations and repeat integration with a smaller time step
				if (h > MIN_TIME_STEP) {
					h *= TIME_STEP_DECREASE_FAC;
					if (solveStatEqSuccess == 0) {
						for (size_t i = 0; i < nlevs; i++) mol->levels[i].pop = oldpops_time[i][0];
					}
					if (ntimesteps > 0) {
						rn = h / h_old;
						BDF_coeffs[0] = (1+2*rn)/(1+rn); BDF_coeffs[1] = -(1+rn); BDF_coeffs[2] = rn*rn/(1+rn);
					}
				} else {
					delete[] A; delete[] pop; delete[] Jac; delete[] dpop_dt;
					oldpops_Ng.clear(); oldpops_time.clear();
					if (cerr_output_iter_progress) binpopfile.close();
					cerr << "#error: cant decrease timestep more, info = " << there_were_bad_levels << "," << solveStatEqSuccess << endl;
					return 1;
				}
			} else {				// proceed the integration over time
				ntimesteps += 1;
				time += h;
				h_old = h;
				if (cerr_output_iter_progress) {
					cerr << "tint= " << ntimesteps << " time= " << time << " h= " << h << " n= " << F_norm << " max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;
					binpopfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
				}
				if (monitor_function < MON_FUN_LOWER_LIMIT) {
					h = min(MAX_TIME_STEP, h * TIME_STEP_INCREASE_FAC);
				}
				for (size_t i = 0; i < nlevs; i++) {
					oldpops_time[i][1] = oldpops_time[i][0];
					oldpops_time[i][0] = mol->levels[i].pop;
					if (cerr_output_iter_progress) {
						double temp_pop = oldpops_time[i][0];
						binpopfile.write(reinterpret_cast<const char*>(&temp_pop), sizeof(double));
					}
					mol->levels[i].pop += dpop_dt[i] * h;
				}
				rn = h / h_old;
				BDF_coeffs[0] = (1+2*rn)/(1+rn); BDF_coeffs[1] = -(1+rn); BDF_coeffs[2] = rn*rn/(1+rn);
			}
		} while (F_norm > MAX_DpopsDt_EPS && ntimesteps < maxNumberOfIterations);

		delete[] A; delete[] pop; delete[] Jac; delete[] dpop_dt;
		oldpops_Ng.clear(); oldpops_time.clear();
		if (cerr_output_iter_progress) binpopfile.close();

		if (ntimesteps == maxNumberOfIterations) cerr << "#warning: maximum number of time steps has been exceeded " << "n= " << F_norm << " max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;

		prepare_results_for_output(LVG_beta);
		
		return 0;
	}

	RT_time_Newt() noexcept : RT()
	{}

	RT_time_Newt(istream & cin) : RT(cin)
	{}

	RT_time_Newt(const unsigned short & initialSolutionSource, const double & MAX_POPS_EPS, const unsigned long & maxNumberOfIterations, const double & beamH, const double & lineWidth) : RT(initialSolutionSource, MAX_POPS_EPS, maxNumberOfIterations, beamH, lineWidth)
	{}
};
