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

	double populate_matrix_vector(double A[], double Jac[], double F[], double B[], beta_LVG& LVG_beta, molModel* mol)	// fill the matrix A, vector B from the statistical equilibrium equations system A*pop=B, and Jacobian Jac from the non-linear system of equations Jac*dpop=-F
	{
		const size_t& n = mol->levels.size();

		// A[i + j*n] - rate of transition from j-th to i-th level
		// Collisional transitions
		// Note that diagonal elements of C were computed in compute_C function in molModel.h, Cii = sum{k=1,Nlevel}(Cik)
		for (size_t i = 0; i < n; i++) {
			A[i + i * n] = -mol->coll_trans[i][i];
			Jac[i + i * n] = A[i + i * n];
			for (size_t j = i + 1; j < n; j++) {
				A[i + j * n] = mol->coll_trans[j][i];
				A[j + i * n] = mol->coll_trans[i][j];
				Jac[i + j * n] = A[i + j * n];
				Jac[j + i * n] = A[j + i * n];
			}
		}

		// Radiative transitions
		double temp_var, temp_var_Jac, Sf, beta, kabs;
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			compute_tau(i, mol);
			compute_J_S_beta(mol, i, LVG_beta, Sf, beta, kabs);
			double common_multiplier = (mol->rad_trans[i].JExt - Sf) * LVG_beta.DbetaDtau(mol->rad_trans[i].tau) +
					dust_HII_CMB_Jext_emission->HII_region_at_LOS * LVG_beta.DbetaHIIDtau_LOS(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII +
					(1 - dust_HII_CMB_Jext_emission->HII_region_at_LOS) * LVG_beta.DbetaHIIDtau_pump(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII;

			const size_t& up = mol->rad_trans[i].up_level;
			const size_t& low = mol->rad_trans[i].low_level;

			temp_var = mol->rad_trans[i].Blu * mol->rad_trans[i].J;
			A[up + low * n] += temp_var;                                   	// A[up][low] += Blu * J
			A[low + low * n] -= temp_var;    									// A[low][low] -= Blu * J

			double blend_demiss_dnlow = 0.0;								// derivative of emission coefficient on population of the lower level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_demiss_dnlow += (int)(low == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].A * mol->rad_trans[i].blends[j].fac;
			}
			double blend_dkabs_dnlow = mol->rad_trans[i].Blu;				// derivative of absorption coefficient on population of the lower level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_dkabs_dnlow += (int)(low == mol->rad_trans[mol->rad_trans[i].blends[j].id].low_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Blu * mol->rad_trans[i].blends[j].fac;
				blend_dkabs_dnlow -= (int)(low == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Bul * mol->rad_trans[i].blends[j].fac;
			}
			double dS_dnlow_beta = HC4PI * invlineWidth * modelPhysPars::n_mol[mol->idspec] * (blend_demiss_dnlow - Sf * blend_dkabs_dnlow);	// derivative of source function on population of the lower level of radiative transition * (1 - beta)
			if (fabs(kabs) > 0.0) dS_dnlow_beta *= (1. - beta) / kabs;
			else dS_dnlow_beta *= 0.5 * modelPhysPars::NdV[mol->idspec] * lineWidth / (modelPhysPars::n_mol[mol->idspec]); // (1-b)/tau -> 0.5 for tau -> 0.0
			double dtau_dnlow = HC4PI * modelPhysPars::NdV[mol->idspec] * blend_dkabs_dnlow;													// derivative of optical depth on population of the lower level of radiative transition
			temp_var_Jac = (mol->levels[low].pop * mol->rad_trans[i].Blu - mol->levels[up].pop * mol->rad_trans[i].Bul) * (dS_dnlow_beta + dtau_dnlow * common_multiplier);
			Jac[up + low * n] += temp_var + temp_var_Jac;
			Jac[low + low * n] -= (temp_var + temp_var_Jac);

			temp_var = (mol->rad_trans[i].A + mol->rad_trans[i].Bul * mol->rad_trans[i].J);
			A[low + up * n] += temp_var;                                   	// A[low][up] += Aul + Bul * J
			A[up + up * n] -= temp_var;       								// A[up][up] -= Aul + Bul * J , where Aul, Bul, Blu - Einstein coefficients, J - mean intensity

			double blend_demiss_dnup = mol->rad_trans[i].A;								// derivative of emission coefficient on population of the upper level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_demiss_dnup += (int)(up == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].A * mol->rad_trans[i].blends[j].fac;
			}
			double blend_dkabs_dnup = -mol->rad_trans[i].Bul;				// derivative of absorption coefficient on population of the upper level of radiative transition taking into account line blending
			for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) {
				blend_dkabs_dnup += (int)(up == mol->rad_trans[mol->rad_trans[i].blends[j].id].low_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Blu * mol->rad_trans[i].blends[j].fac;
				blend_dkabs_dnup -= (int)(up == mol->rad_trans[mol->rad_trans[i].blends[j].id].up_level && mol->idspec == mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[mol->rad_trans[i].blends[j].id].Bul * mol->rad_trans[i].blends[j].fac;
			}
			double dS_dnup_beta = HC4PI * invlineWidth * modelPhysPars::n_mol[mol->idspec] * (blend_demiss_dnup - Sf * blend_dkabs_dnup);	// derivative of source function on population of the upper level of radiative transition * (1 - beta)
			if (fabs(kabs) > 0.0) dS_dnup_beta *= (1. - beta) / kabs;
			else dS_dnup_beta *= 0.5 * modelPhysPars::NdV[mol->idspec] * lineWidth / (modelPhysPars::n_mol[mol->idspec]); // (1-b)/tau -> 0.5 for tau -> 0.0
			double dtau_dnup = HC4PI * modelPhysPars::NdV[mol->idspec] * blend_dkabs_dnup;													// derivative of optical depth on population of the upper level of radiative transition
			temp_var_Jac = (mol->levels[up].pop * mol->rad_trans[i].Bul - mol->levels[low].pop * mol->rad_trans[i].Blu) * (dS_dnup_beta + dtau_dnup * common_multiplier);
			Jac[low + up * n] += temp_var + temp_var_Jac;
			Jac[up + up * n] -= (temp_var + temp_var_Jac);
		}

		double pops_sum = 0.0;
		double Fnorm = 0.0;
		for (size_t i = n; i-- > 1; ) {
			F[i] = 0.0;
			for (size_t j = n; j-- > 0; ) F[i] += A[i + j * n] * mol->levels[j].pop;
			F[i] = - F[i];
			Jac[0 + i * n] = 1.0;
			pops_sum += mol->levels[i].pop;
			A[0 + i * n] = A[0];
			B[i] = 0.0;
			Fnorm += F[i] * F[i];
		}
		Jac[0] = 1.0;
		F[0] = this->partition_function_ratio[mol->idspec] - (pops_sum + mol->levels[0].pop);
		Fnorm += F[0] * F[0];
		B[0] = A[0] * this->partition_function_ratio[mol->idspec]; // the sum of populations should be = partition functions ratio or = 1 multiplied by A[0][0] for numerical stability
		return sqrt(Fnorm);
	}

	double getF(double A[], double pop[], double temp_pop[], double dpop[], const double & rate, beta_LVG& LVG_beta, molModel* mol)	// fill the matrix A, vector B from the statistical equilibrium equations system A*pop=B, and Jacobian Jac from the non-linear system of equations Jac*dpop=-F
	{
		const size_t& n = mol->levels.size();
		for (size_t i = 0; i < n; i++) {
			pop[i] = mol->levels[i].pop + rate * dpop[i];
			temp_pop[i] = mol->levels[i].pop;
			mol->levels[i].pop = pop[i];
		}

		// A[i + j*n] - rate of transition from j-th to i-th level
		// Collisional transitions
		// Note that diagonal elements of C were computed in compute_C function in molModel.h, Cii = sum{k=1,Nlevel}(Cik)
		for (size_t i = 0; i < n; i++) {
			A[i + i * n] = -mol->coll_trans[i][i];
			for (size_t j = i + 1; j < n; j++) {
				A[i + j * n] = mol->coll_trans[j][i];
				A[j + i * n] = mol->coll_trans[i][j];
			}
		}

		// Radiative transitions
		double temp_var, Sf, beta, kabs;
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			compute_tau(i, mol);
			compute_J_S_beta(mol, i, LVG_beta, Sf, beta, kabs);
			
			const size_t& up = mol->rad_trans[i].up_level;
			const size_t& low = mol->rad_trans[i].low_level;

			temp_var = mol->rad_trans[i].Blu * mol->rad_trans[i].J;
			A[up + low * n] += temp_var;                                   	// A[up][low] += Blu * J
			A[low + low * n] -= temp_var;    									// A[low][low] -= Blu * J

			temp_var = (mol->rad_trans[i].A + mol->rad_trans[i].Bul * mol->rad_trans[i].J);
			A[low + up * n] += temp_var;                                   	// A[low][up] += Aul + Bul * J
			A[up + up * n] -= temp_var;       								// A[up][up] -= Aul + Bul * J , where Aul, Bul, Blu - Einstein coefficients, J - mean intensity
		}

		double pops_sum = 0.0;
		double Fnorm = 0.0;
		double F;
		for (size_t i = n; i-- > 1; ) {
			F = 0.0;
			for (size_t j = n; j-- > 0; ) F += A[i + j * n] * mol->levels[j].pop;
			pops_sum += mol->levels[i].pop;
			Fnorm += F * F;
		}
		double F0 = this->partition_function_ratio[mol->idspec] - (pops_sum + mol->levels[0].pop);
		Fnorm += F0 * F0;

		for (size_t i = 0; i < n; i++) mol->levels[i].pop = temp_pop[i];

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

		vector <double*> pop;
		pop.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			pop.push_back(new double[mols[ispec].levels.size()]);
		}

		vector <double*> dpop;
		dpop.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			dpop.push_back(new double[mols[ispec].levels.size()]);
		}

		vector <double*> A; 	// reserve space for matrix A from the statistical equilibrium equations system A*X=B
		A.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			A.push_back(new double[mols[ispec].levels.size()*mols[ispec].levels.size()]);
		}

		vector <double*> Jac; 	// reserve space for Jacobian
		Jac.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			Jac.push_back(new double[mols[ispec].levels.size()*mols[ispec].levels.size()]);
		}

		double MaxRPopDiff;
		size_t levelWithMaxRPopDiff;
		size_t speciesWithMaxRPopDiff;

		unsigned int iter = 1;
		double Fnorm = 0.0;
		do {
			Fnorm = 0.0;
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				double tempFnorm = populate_matrix_vector(A[ispec], Jac[ispec], dpop[ispec], pop[ispec], LVG_beta, &mols[ispec]);
				if (tempFnorm > Fnorm) Fnorm = tempFnorm;
				//double cond_number = get_condition_number(A[ispec], &mols[ispec]);
				size_t solveStatEqSuccess = solve_eq_sys(Jac[ispec], dpop[ispec], &mols[ispec]);		// solve the equation system Jac*dpop = -F
				if (solveStatEqSuccess != 0) {						// the solution can't be found
					for (size_t ispec1 = 0; ispec1 < modelPhysPars::nSpecies; ispec1++) {
						delete[] pop[ispec1];
						delete[] dpop[ispec1];
						delete[] A[ispec1];
						delete[] Jac[ispec1];
					}
					oldpops_Ng.clear();
					cerr << "#error: Newton solve_eq_sys failed, info = " << solveStatEqSuccess << endl;
					return 1;
				}
				double rate = 0.99; // see kinsol package for the step length calculations according to the bounds on the solution vector, https://github.com/LLNL/sundials/tree/main/src/kinsol
				for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
					double cur_pop = dpop[ispec][i] + mols[ispec].levels[i].pop;
					if ((cur_pop <= 0.0 || cur_pop > 1.0) && fabs(dpop[ispec][i]) > 0.0) {
						rate = min(rate, mols[ispec].levels[i].pop / fabs(dpop[ispec][i]));
					}
				}
				if (rate < 0.99) rate *= 0.9;
				if (rate >= MIN_NEWT_SCALE) {
					double rate1 = MIN_NEWT_SCALE;
					double rate2 = rate;
					double rate_step = (rate2 - rate1) * MAX_NEWT_SCALE_STEP;
					do {
						double Fnorm1 = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], rate1, LVG_beta, &mols[ispec]);
						double Fnorm2 = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], (rate1 + 2 * rate_step), LVG_beta, &mols[ispec]);
						double DFnormDrate1 = 0.5 * (Fnorm2 - Fnorm1) / rate_step;
						Fnorm1 = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], (rate2 - 2 * rate_step), LVG_beta, &mols[ispec]);
						Fnorm2 = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], rate2, LVG_beta, &mols[ispec]);
						double DFnormDrate2 = 0.5 * (Fnorm2 - Fnorm1) / rate_step;
						if (DFnormDrate1 * DFnormDrate2 > 0.0 || DFnormDrate1 > 0.0 || DFnormDrate2 < 0.0) break;
						rate = (rate1 + rate2) * 0.5;
						Fnorm1 = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], (rate - rate_step), LVG_beta, &mols[ispec]);
						Fnorm2 = getF(A[ispec], Jac[ispec], pop[ispec], dpop[ispec], (rate + rate_step), LVG_beta, &mols[ispec]);
						double DFnormDrate = 0.5 * (Fnorm2 - Fnorm1) / rate_step;
						if (DFnormDrate < 0.0) rate1 = rate;
						else rate2 = rate;
					} while (2 * fabs(rate2 - rate1) > NEWT_SCALE_ACCURACY * (rate2 + rate1));
					for (size_t i = 0; i < mols[ispec].levels.size(); i++) pop[ispec][i] = rate * dpop[ispec][i] + mols[ispec].levels[i].pop;
				} else { // use simple iteration if Newton step is too small
					solveStatEqSuccess = solve_eq_sys(A[ispec], pop[ispec], &mols[ispec]);				// solve statistical equilibrium equation with LU decomposition
					if (solveStatEqSuccess != 0) {						// the solution can't be found
						for (size_t ispec1 = 0; ispec1 < modelPhysPars::nSpecies; ispec1++) {
							delete[] pop[ispec1];
							delete[] dpop[ispec1];
							delete[] A[ispec1];
							delete[] Jac[ispec1];
						}
						oldpops_Ng.clear();
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
		} while (Fnorm > FNORM_STOPPING && MaxRPopDiff > MAX_DpopsDt_EPS && iter <= maxNumberOfIterations);

		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			delete[] pop[ispec];
			delete[] dpop[ispec];
			delete[] A[ispec];
			delete[] Jac[ispec];
		}
		oldpops_Ng.clear();

		if (iter > maxNumberOfIterations) cerr << "#warning: maximum number of iterations has exceeded Fnorm= " << Fnorm << " max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;
		
		prepare_results_for_output(LVG_beta);
		
		return 0;
	}

	RT_statEquiv() : RT()
	{}

	RT_statEquiv(istream & cin) : RT(cin)
	{}
};
