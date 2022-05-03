#pragma once
#include "RT.h"
#include <cstring>

class RT_point_iterations : public RT		// solves statistical equilibrium equations for level populations using fixed point iterations
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

	void populate_matrix_vector(double A[], double B[], beta_LVG & LVG_beta, molModel *mol)	// fill the matrix A and vector B from the statistical equilibrium equations system A*X=B
	{
		const size_t & n = mol->levels.size();
		
		// Collisional transitions
		// Note that diagonal elements of C were computed in compute_C function in molModel.h, Cii = sum{k=1,Nlevel}(Cik)
		for (size_t i = 0; i < n; i++) {
			A[i + i*n] = - mol->coll_trans[i][i];
			for (size_t j = i+1; j < n; j++) {
				A[i + j*n] = mol->coll_trans[j][i];
				A[j + i*n] = mol->coll_trans[i][j];
			}
		}
		
        // Radiative transitions
		double temp_var, dummy_S, dummy_beta, dummy_betaS;
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			compute_tau(i, mol);
			compute_J_S_beta(mol, i, LVG_beta, dummy_S, dummy_beta, dummy_betaS);
            
			const size_t & up = mol->rad_trans[i].up_level;
			const size_t & low = mol->rad_trans[i].low_level;
			
			temp_var = mol->rad_trans[i].Blu * mol->rad_trans[i].J;
			A[up + low*n]  += temp_var;                                   	// A[up][low] += Blu * J
			A[low + low*n] -= temp_var;    									// A[low][low] -= Blu * J

			temp_var = (mol->rad_trans[i].A + mol->rad_trans[i].Bul * mol->rad_trans[i].J);
			A[low + up*n]  += temp_var;                                   	// A[low][up] += Aul + Bul * J
			A[up + up*n] -= temp_var;       								// A[up][up] -= Aul + Bul * J , where Aul, Bul, Blu - Einstein coefficients, J - mean intensity
		}
		
		for (size_t i = 0; i < n; i++) {
			A[0 + i*n] = 1.0;	// A[0][i] = 1.0 - the equation for the first level is replaced by the particle number conservation law, i.e. the sum of populations should be = 1 or = partition functions ratio
			B[i] = 0.0;
		}
		B[0] = this->partition_function_ratio[mol->idspec]; // the sum of populations should be = 1 or = partition functions ratio
	}

	void prepare_results_for_output(beta_LVG & LVG_beta)
	{
		// preparing quantities for output
		double dummy_S = 0.0;
		double dummy_beta = 0.0;
		double dummy_betaS = 0.0;
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			for (size_t i = 0; i < mols[ispec].rad_trans.size(); i++) {
				compute_tau(i, &mols[ispec]); 		// computing final optical depths that can be used for output
				compute_J_S_beta(&mols[ispec], i, LVG_beta, dummy_S, dummy_beta, dummy_betaS); 	// computing final mean intensities that can be used for output
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
				mols[ispec].rad_trans[i].kabs_dust = mols[ispec].rad_trans[i].taud_in * modelPhysPars::Hdens * invlineWidth / modelPhysPars::max_NH2dV; //absorption coefficient of the dust inside the maser region
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

		vector <double*> A; 	// reserve space for matrix A from the statistical equilibrium equations system A*X=B
		A.reserve(modelPhysPars::nSpecies);
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			A.push_back(new double[mols[ispec].levels.size()*mols[ispec].levels.size()]);
		}

		double MaxRPopDiff;
		size_t levelWithMaxRPopDiff;
		size_t speciesWithMaxRPopDiff;

		unsigned int iter = 0;
		do {
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				populate_matrix_vector(A[ispec], pop[ispec], LVG_beta, &mols[ispec]);
				//double cond_number = get_condition_number(A[ispec], &mols[ispec]);
				size_t solveStatEqSuccess = solve_eq_sys(A[ispec], pop[ispec], &mols[ispec]);	// solve statistical equilibrium equation with LU decomposition, solution is stored in mol->levels[i].pop
				if (solveStatEqSuccess != 0) {						// the solution can't be found
					for (size_t ispec1 = 0; ispec1 < modelPhysPars::nSpecies; ispec1++) {
						delete[] pop[ispec1];
						delete[] A[ispec1];
					}
					oldpops_Ng.clear();
					cerr << "#error: solve_stat_equilibrium failed, info = " << solveStatEqSuccess << endl;
					return 1;
				}
			}
			iter += 1;
			double dummy_pop_norm = 1.e-30;
			levelWithMaxRPopDiff = 0;
			speciesWithMaxRPopDiff = 0;
			MaxRPopDiff = -1.0;
			for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
				update_check_pops(ispec, mols, pop[ispec], speciesWithMaxRPopDiff, levelWithMaxRPopDiff, MaxRPopDiff, iter, oldpops_Ng[ispec], dummy_pop_norm);
			}
			if (cerr_output_iter_progress) {
				cerr << iter << " max.dev.= " << MaxRPopDiff << " mol/level with max.dev.= " << speciesWithMaxRPopDiff << " / " <<levelWithMaxRPopDiff << endl;
			}
		} while (MaxRPopDiff > MAX_DpopsDt_EPS && iter < maxNumberOfIterations);

		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			delete[] pop[ispec];
			delete[] A[ispec];
		}
		oldpops_Ng.clear();

		if (iter == maxNumberOfIterations) cerr << "#warning: maximum number of iterations has exceeded max.dev.= " << MaxRPopDiff << " level with max.dev.= " << levelWithMaxRPopDiff << endl;
		
		prepare_results_for_output(LVG_beta);
		
		return 0;
	}

	RT_point_iterations() : RT()
	{}

	RT_point_iterations(istream & cin) : RT(cin)
	{}
};

/*
void compute_J_S_beta(const size_t & i, beta_LVG & LVG_beta, double & S, double & beta)		//computes mean intensity, source function and escape probability for radiative transition i
	{ //see also equation for Jav in Appendix A of Sobolev et al. 1997
		auto kabsf = [&](const size_t & i) { //this lambda function returns absorption coefficient for i-th transition without taking into account line overlapping
			const size_t & up = mol->rad_trans[i].up_level;
			const size_t & low = mol->rad_trans[i].low_level;
			return mol->levels[low].pop * mol->rad_trans[i].Blu - mol->levels[up].pop * mol->rad_trans[i].Bul;
		};
		auto emissf = [&](const size_t & i) { //this lambda function returns emission coefficient for i-th transition without taking into account line overlapping
			return mol->rad_trans[i].A * mol->levels[mol->rad_trans[i].up_level].pop;
		};

		double kabs = kabsf(i); //part of the line absorption coefficient 
		double emiss = emissf(i); //part of the line emission coefficient
		for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) { //take into account line overlapping
			kabs += kabsf(mol->rad_trans[i].blends[j].id) * mol->rad_trans[i].blends[j].fac;
			emiss += emissf(mol->rad_trans[i].blends[j].id) * mol->rad_trans[i].blends[j].fac;
		}

		S = 0.0; //source function
		if (fabs(kabs) > 0.0) S = emiss / kabs;

		beta = LVG_beta.beta(mol->rad_trans[i].tau); //escape probability = beta(tau)
		mol->rad_trans[i].Jin = (1. - beta) * S; //internal radiation
		mol->rad_trans[i].J = mol->rad_trans[i].Jin + beta * mol->rad_trans[i].JExt + LVG_beta.betaHII(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII; //sum of internal and external radiation
	}

void get_A_and_Jacobian(double A[], double Jac[], beta_LVG & LVG_beta)	//fill the matrix A from the statistical equilibrium equation A*X=B
{
    const size_t & lda = mol->levels.size();
    double * err_cs = new double[lda];
    
    //Collisional transitions
    //Note that diagonal elements of C were computed in compute_C function in molModel.h, Cii = - sum{k=1,Nlevel}(Cik)
    for (size_t i = 0; i < mol->levels.size(); i++) {
        for (size_t j = 0; j < mol->levels.size(); j++) {
            A[i + j*lda] = mol->coll_trans[j][i].C;
            Jac[i + j*lda] = 0.0;
        }
        err_cs[i] = 0.0;
    }
    
    //Radiative transitions
    size_t up, low;
    double beta; // = escape probability
    double Sf; //source function
    double aux_var;
    for (size_t i = 0; i < mol->rad_trans.size(); i++) {
        compute_tau(i);
        compute_J_S_beta(i, LVG_beta, Sf, beta);

        up = mol->rad_trans[i].up_level;
        low = mol->rad_trans[i].low_level;

        add_rad_trans_to_A(A, err_cs, lda, i, up, low);

        aux_var = LVG_beta.tauDbetaDtau(mol->rad_trans[i].tau) * (Sf - mol->rad_trans[i].JExt - mol->rad_trans[i].JExtHII);

        Jac[up + low*lda]  = A[up + low*lda]  - mol->rad_trans[i].Blu * (mol->rad_trans[i].Jin + aux_var);		// Jac[up][low] -= Blu * J + Blu*tau*DbetaDtau * (S-Jext)
        Jac[low + up*lda]  = A[low + up*lda]  - (mol->rad_trans[i].A * (1.- beta) + mol->rad_trans[i].Bul * (mol->rad_trans[i].Jin + aux_var)); // Jac[low][up] -= (1-beta)*Aul + Bul * J + Bul*tau*DbetaDtau * (S-Jext)
        Jac[up + up*lda]   = A[up + up*lda]   + (mol->rad_trans[i].A * (1.- beta) + mol->rad_trans[i].Bul * (mol->rad_trans[i].Jin + aux_var)); // Jac[up][up] += (1-beta)*Aul + Bul * J + Bul*tau*DbetaDtau * (S-Jext)
        Jac[low + low*lda] = A[low + low*lda] + mol->rad_trans[i].Blu * (mol->rad_trans[i].Jin + aux_var); // Jac[low][low] += Blu * J + Blu*tau*DbetaDtau * (S-Jext)
    }
    
    for (size_t i = 0; i < mol->levels.size(); i++) {
        A[i + i*lda] += err_cs[i];
        Jac[i + i*lda] += err_cs[i];
        A[0 + i*lda] = 1.0e0;  // A[0][i] = 1.0 - the equation for the first level is replaced by particle number conservation law, i.e. the sum of populations should be = 1
        Jac[0 + i*lda] = 1.0e0;
    }
    delete[] err_cs;
}

double *Jac = new double[mol->levels.size()*mol->levels.size()]; //reserve space for Jacobian
double F_norm0;
bool there_are_bad_levels = false;
iter = 0;

do { //finds solution of statistical equilibrium equations A*X=B with Newton-Rapshon method
    get_A_and_Jacobian(A, Jac, LVG_beta);		//fill the matrix A from the statistical equilibrium equation A*X=B and Jacobian
    F_norm0 = get_F(A, pop);
    solveStatEqSuccess = solve_eq_sys(Jac, pop);
    if (solveStatEqSuccess != 0) {	//the solution can't be found
        delete[] A; delete[] Jac; delete[] pop;
        cerr << "#error: solve_stat_equilibrium failed, info = " << solveStatEqSuccess << endl;
        return 1;
    }

    double rate  = 10.0;
    do {
        rate *= 0.1;
        there_are_bad_levels = false;
        for (size_t i = 0; i < mol->levels.size(); i++) {
            mol->levels[i].pop = mol->levels[i].oldpop + rate * pop[i];
            if (mol->levels[i].pop >= (1.0-1000*DBL_EPSILON) || mol->levels[i].pop <= 0.0) {
                mol->levels[i].pop = mol->levels[i].oldpop;
                pop[i] = 0.0;
                there_are_bad_levels = true;
            }
        }
        populate_matrix(A, LVG_beta, false);
        F_norm = get_F_norm(A);
        for (size_t i = 0; i < mol->levels.size(); i++) mol->levels[i].pop = mol->levels[i].oldpop;
        if (rate < MIN_NEWTON_LENGTH) break;
    } while (F_norm > F_norm0 * (1. - 1.e-4*rate) || there_are_bad_levels); //Armijo rule
    

    iter += 1;
    convAccel.doConvCheck_newton(pop, mol, rate, LTE_pops);
    if (cerr_output_iter_progress) {
        cerr << "newt rate= " << F_norm0 << " iter= " << iter << " max.dev.= " << convAccel.MaxRPopDiff << " level with max.dev.= " << convAccel.levelWithMaxRPopDiff << endl;
    }
} while (iter < MAX_NEWTON_ITER_NUMBER);
delete[] Jac;
*/
