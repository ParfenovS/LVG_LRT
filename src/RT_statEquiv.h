#pragma once
#include "RT.h"
#include <cstring>

class RT_statEquiv : public RT		// solves statistical equilibrium equations for level populations using Newton method or fixed point iterations
{
private:
	
	double get_condition_number(double A0[], molModel *mol)
	{
		const size_t & n = mol->levels.size();
		double *s = new double[n];
		double *e = new double[n];
		double *work = new double[n];
		double *A = new double[n*n];
		std::memcpy(A, A0, n*n*sizeof(double));
		dsvdc ( A, n, n, s, e, work);
		double cond_number = s[0] / s[n-1];
		delete[] s;
		delete[] e;
		delete[] work;
		delete[] A;
		return cond_number;
	}

	double populate_matrix_vector(vector <double> & A, vector <double> & Jac, vector <double> & F, vector <double> & B, const vector <double> & fscale2, const size_t & maxPopid, beta_LVG& LVG_beta, molModel* mol)	// fill the matrix A, vector B from the statistical equilibrium equations system A*pop=B, and Jacobian Jac from the non-linear system of equations Jac*dpop=-F
	{
		get_A_Jac(A, Jac, LVG_beta, mol);

		const size_t& n = mol->levels.size();

		double pops_sum = 0.0;
		double Fnorm = 0.0;
		for (size_t i = n; i-- > 0; ) {
			if (i != maxPopid) {
				F[i] = 0.0;
				for (size_t j = n; j-- > 0; ) {
					F[i] += A[i + j * n] * mol->levels[j].pop;
				}
				F[i] = - F[i];
				B[i] = 0.0;
				Fnorm += F[i] * F[i] * fscale2[i];
			}
			pops_sum += mol->levels[i].pop;
			Jac[maxPopid + i * n] = 1.0;
			A[maxPopid + i * n] = 1.0;
		}
		F[maxPopid] = this->partition_function_ratio[mol->idspec] - pops_sum;
		Fnorm += F[maxPopid] * F[maxPopid] * fscale2[maxPopid];
		B[maxPopid] = this->partition_function_ratio[mol->idspec]; // the sum of populations should be = partition functions ratio or = 1 multiplied by A[0][0] for numerical stability
		return sqrt(Fnorm);
	}

	double getF(vector <double> & A, vector <double> & pop, vector <double> & temp_pop, const vector <double> & dpop, const double & rate, const vector <double> & fscale2, const size_t & maxPopid, beta_LVG& LVG_beta, molModel* mol)	// fill the matrix A, vector B from the statistical equilibrium equations system A*pop=B, and Jacobian Jac from the non-linear system of equations Jac*dpop=-F
	{
		const size_t& n = mol->levels.size();

		for (size_t i = 0; i < n; i++) {
			pop[i] = fma(dpop[i], rate, mol->levels[i].pop);
			temp_pop[i] = mol->levels[i].pop;
			mol->levels[i].pop = pop[i];
		}

		get_A(A, LVG_beta, mol);

		double Fnorm = 0.0;
		double F;
		double pops_sum = 0.0;
		for (size_t i = n; i-- > 0; ) {
			if (i != maxPopid) {
				F = 0.0;
				for (size_t j = n; j-- > 0; ) F += A[i + j * n] * pop[j];
				Fnorm += F * F * fscale2[i];
			}
			pops_sum += pop[i];
			mol->levels[i].pop = temp_pop[i];
		}
		F = this->partition_function_ratio[mol->idspec] - pops_sum;
		Fnorm += F * F * fscale2[maxPopid];
		return sqrt(Fnorm);
	}

	void prepare_results_for_output(beta_LVG & LVG_beta)
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
				compute_brightness_temperature(i, &mols[ispec]); 	// computing brightness temperature and intensity of the emission
			}
		}
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				// the output optical depth is an optical depth along the line of sight, thus one need to take into account beaming
				mols[ispec].rad_trans[i].tau *= beamH;
			}
		}
	}

public:

	int radiative_transfer() override		// main function which should be called to perform radiative transfer calculations
	{
		initial_solution();							// get the initial values of populations
		
		// set the external emission mean intensity, optical depth, absorption and emission coefficients
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				mols[ispec].rad_trans[i].JExt = dust_HII_CMB_Jext_emission->compute_Jext_dust_CMB_file(mols[ispec].rad_trans[i].nu); //external emission from dust, cosmic microwave background or file
				mols[ispec].rad_trans[i].JExtHII = dust_HII_CMB_Jext_emission->compute_JextHII(mols[ispec].rad_trans[i].nu); //external emission from HII region, should be separated from other types of emission because of maser beaming
				mols[ispec].rad_trans[i].taud_in = dust_HII_CMB_Jext_emission->tau_dust_in(mols[ispec].rad_trans[i].nu, lineWidth); //optical depth of the dust inside the maser region
				mols[ispec].rad_trans[i].kabs_dust = mols[ispec].rad_trans[i].taud_in * modelPhysPars::H2dens * invlineWidth / modelPhysPars::max_NH2dV; //absorption coefficient of the dust inside the maser region
				mols[ispec].rad_trans[i].emiss_dust = mols[ispec].rad_trans[i].kabs_dust * dust_HII_CMB_Jext_emission->inner_dust_source_function(mols[ispec].rad_trans[i].nu); //emission coefficient of the dust inside the maser region
				// mols[ispec].rad_trans[i].JExtHII will be zero if external emission will be taken from file
			}
		}

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

		vector <vector <double> > pop(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			pop[ispec].resize(mols[ispec].levels.size());
		}

		vector <vector <double> > dpop(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			dpop[ispec].resize(mols[ispec].levels.size());
		}

		vector <vector <double> > A(modelPhysPars::nSpecies); // reserve space for matrix A from the statistical equilibrium equations system A*X=B
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			A[ispec].resize(mols[ispec].levels.size() * mols[ispec].levels.size());
		}

		vector <vector <double> > Jac(modelPhysPars::nSpecies); // reserve space for Jacobian
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			Jac[ispec].resize(mols[ispec].levels.size() * mols[ispec].levels.size());
		}

		vector <double> Fnorm_old(modelPhysPars::nSpecies);

		vector <size_t> max_pop_id(modelPhysPars::nSpecies);
		vector <vector <double> > fscale2(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			Fnorm_old[ispec] = 1.e60;
			max_pop_id[ispec] = 0;
			fscale2[ispec].resize(mols[ispec].levels.size());
			for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
				double max_pop = mols[ispec].levels[0].pop;
				if (mols[ispec].levels[i].pop > max_pop) {
					max_pop = mols[ispec].levels[i].pop;
					max_pop_id[ispec] = i;
				}
				fscale2[ispec][i] = 1.0;
			}
			populate_matrix_vector(A[ispec], Jac[ispec], dpop[ispec], pop[ispec], fscale2[ispec], max_pop_id[ispec], LVG_beta, &mols[ispec]);
			if (solve_eq_sys(Jac[ispec].data(), dpop[ispec].data(), &mols[ispec]) == 0) {
				for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
					double temp_scale = dpop[ispec][i] * dpop[ispec][i];
					if (temp_scale > 0.0) fscale2[ispec][i] = 1.0 / temp_scale;
				}
			}
		}

		double MaxRPopDiff;
		size_t levelWithMaxRPopDiff;
		size_t speciesWithMaxRPopDiff;

		unsigned int iter = 1;
		double Fnorm = 0.0;
		bool there_are_blended_lines = false;
		do {
			Fnorm = 0.0;
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				double tempFnorm = populate_matrix_vector(A[ispec], Jac[ispec], dpop[ispec], pop[ispec], fscale2[ispec], max_pop_id[ispec], LVG_beta, &mols[ispec]);
				if (tempFnorm > Fnorm) Fnorm = tempFnorm;
				Fnorm_old[ispec] = tempFnorm;
				//double cond_number = get_condition_number(A[ispec].data(), &mols[ispec]);
				size_t solveStatEqSuccess = solve_eq_sys(Jac[ispec].data(), dpop[ispec].data(), &mols[ispec]);		// solve the equation system Jac*dpop = -F
				if (solveStatEqSuccess != 0) {						// the solution can't be found
					cerr << "#error: Newton solve_eq_sys failed, info = " << solveStatEqSuccess << endl;
					return 1;
				}
				double rate = 1.0;
				// see kinsol package for the step length calculations according to the bounds on the solution vector, https://github.com/LLNL/sundials/tree/main/src/kinsol
				for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
					double cur_pop = dpop[ispec][i] + mols[ispec].levels[i].pop;
					if ((cur_pop <= 0.0 || cur_pop > 1.0) && fabs(dpop[ispec][i]) > 0.0) {
						rate = min(rate, mols[ispec].levels[i].pop / fabs(dpop[ispec][i]));
					}
				}
				if (rate < 1.0) rate *= 0.9;
				bool switch_to_simple_iteration = true;
				if (rate > MIN_NEWT_SCALE) {
					double rate1 = MIN_NEWT_SCALE;
					double rate2 = rate;
					double rate_step = (rate2 - rate1) / NUM_OF_NEWT_SCALE_GRID_POINTS;
					double temp_rate = rate1;
					double min_Fnorm = 1.e100;

					for (size_t refine_step = 0; refine_step < NUM_OF_NEWT_SCALE_GRID_REFINEMENTS; ++refine_step) {
						do {
							double temp_Fnorm = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], temp_rate, fscale2[ispec], max_pop_id[ispec], LVG_beta, &mols[ispec]);
							if (temp_Fnorm < min_Fnorm) {
								min_Fnorm = temp_Fnorm;
								rate = temp_rate;
							}
							temp_rate += rate_step;
						} while (temp_rate <= rate2);
						rate1 = max(rate1, rate - rate_step);
						rate2 = min(rate2, rate + rate_step);
						rate_step = (rate2 - rate1) / NUM_OF_NEWT_SCALE_GRID_POINTS;
						temp_rate = rate1;
						min_Fnorm = 1.e100;
					}

					double temp_Fnorm = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], rate, fscale2[ispec], max_pop_id[ispec], LVG_beta, &mols[ispec]);
					if (Fnorm_old[ispec] > (MIN_FNORM_RELATIVE_DECREASE + 1) * temp_Fnorm) { 
						for (size_t i = 0; i < mols[ispec].levels.size(); i++) pop[ispec][i] = fma(dpop[ispec][i], rate, mols[ispec].levels[i].pop);
						switch_to_simple_iteration = false;
					}
				}
				if (switch_to_simple_iteration) { // use simple iteration if Newton step is too small or relative change in norm of Dpop/Dt is too small
					populate_matrix_vector(A[ispec], Jac[ispec], dpop[ispec], pop[ispec], fscale2[ispec], max_pop_id[ispec], LVG_beta, &mols[ispec]);
					solveStatEqSuccess = solve_eq_sys(A[ispec].data(), pop[ispec].data(), &mols[ispec]);				// solve statistical equilibrium equation with LU decomposition
					if (solveStatEqSuccess != 0) {						// the solution can't be found
						cerr << "#error: simple iteration solve_eq_sys failed, info = " << solveStatEqSuccess << endl;
						return 1;
					}
				}
			}
			double dummy_pop_norm = 1.e-30;
			levelWithMaxRPopDiff = 0;
			speciesWithMaxRPopDiff = 0;
			MaxRPopDiff = -1.0;
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				update_check_pops(&mols[ispec], pop[ispec], speciesWithMaxRPopDiff, levelWithMaxRPopDiff, MaxRPopDiff, iter, oldpops_Ng[ispec], dummy_pop_norm);
			}
			if (cerr_output_iter_progress) {
				cerr << iter << " Fnorm= " << Fnorm << " max.dev.= " << MaxRPopDiff << " mol/level with max.dev.= " << speciesWithMaxRPopDiff << " / " <<levelWithMaxRPopDiff << endl;
			}
			iter += 1;
			if (Fnorm <= MAX_DpopsDt_EPS && !there_are_blended_lines) {
				if (lineWidth > DBL_EPSILON) there_are_blended_lines = find_blends(); // find overlapping lines
				if (there_are_blended_lines) Fnorm = 1.1 * MAX_DpopsDt_EPS;
			} 
		} while (Fnorm > MAX_DpopsDt_EPS && iter <= maxNumberOfIterations);

		if (iter > maxNumberOfIterations) cerr << "#warning: maximum number of iterations has exceeded Fnorm= " << Fnorm << " max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;
		
		prepare_results_for_output(LVG_beta);
		
		return 0;
	}

	RT_statEquiv() : RT()
	{}

	RT_statEquiv(istream & cin) : RT(cin)
	{}
};
