#pragma once
#include "Dust_HII_CMB_Jext_Radiation.h"
#include "beta_LVG.h"
#include "linpack/linpack_d.hpp"

class RT		// base class containing functions for radiative transfer calculations;
{
protected:

	unsigned short initialSolutionSource;						// =0 - initial solution is taken from file ; =1 - initial solution is computed for optically thin case ; =2 - LTE initial solution
	double MAX_DpopsDt_EPS;										// the computations will stop if the maximum length of the vector with components of Dn/Dt (derivative of populations with time) or the relative populations difference between point iterations <= MAX_DpopsDt_EPS
	unsigned long maxNumberOfIterations;						// maximum number of iterations
	double beamH; 												// beaming factor; this is eps^-1 = D(ln r)/D(ln V) quantity from e.g. Sobolev et al. 1997, Cragg et al. 2005
	double lineWidth; 											// spectral width of the lines; this will be used to determine blended lines
	double invlineWidth;										// = 1 / lineWidth
	dust_HII_CMB_Jext_radiation *dust_HII_CMB_Jext_emission; 	// this object will be used for calculations of external emission
	bool cerr_output_iter_progress; 							// = true - the current iteration number, maximum relative pop difference and corresponding level number will be printed on the standard cerr pipe
	string line_profile_shape; 									// ="g" - Gaussian line profile; ="r" - rectangular line profile
	const vector<double> input_full_partition_function = PARTITION_FUNCTION;	// full partition rotational function (see hiddenParameters.h)

	vector <double> full_partition_function;
	vector <double> partition_function_ratio;
	
	template <typename T>
	void read_parameters(T & fin) 								// read parameters of radiative transfer calculations from Parameters/RadiativeTransfer.txt file or console input
	{
		istringstream sfin;
		string str;

		if (!fin.good()) 	// check if we found the file or that console input is not corrupted
			throw runtime_error("can't read radiative transfer parameters; check Parameters/RadiativeTransfer.txt file or console input");

		getline(fin, str);
		getline(fin, str);
		initialSolutionSource = readline<unsigned short>(fin);	// readline function is taken from auxiliary.h

		getline(fin, str);
		for (size_t ispec = 0; ispec < mols.size(); ispec++) {
			getline(fin, str);
			filename_pops_in.push_back(trim(str));	// trim function is taken from auxiliary.h
		}

		getline(fin, str);
		for (size_t ispec = 0; ispec < mols.size(); ispec++) {
			getline(fin, str);
			filename_pops_out.push_back(trim(str));
		}

		getline(fin, str);
		for (size_t ispec = 0; ispec < mols.size(); ispec++) {
			getline(fin, str);
			filename_lamda.push_back(trim(str));
		}

		for (short i = 0; i < 4; i++) getline(fin, str);
		sfin.str(trim(str));
		sfin >> MAX_DpopsDt_EPS >> maxNumberOfIterations;
		sfin.clear();
		if (MAX_DpopsDt_EPS < DBL_EPSILON)
			throw runtime_error("maximum length of Dn/Dt vector in Parameters/RadiativeTransfer.txt file or console input should be larger than machine epsilon");

		for (short i = 0; i < 3; i++) getline(fin, str);
		beamH = readline<double>(fin);

		for (short i = 0; i < 3; i++) getline(fin, str);
		sfin.str(trim(str));
		sfin >> lineWidth >> line_profile_shape;
		lineWidth *= 1.e5; // [km/s] -> [cm/s]
		invlineWidth = 1. / lineWidth;
	}

	/*
** reading file with levels population
*/
	void read_pops(string filename, molModel *mod)
	{
		string str;
		ifstream fin;
		istringstream sfin;

		fin.open(filename.c_str(), ios::in);
		sfin.clear();

		if ( !fin.good() ) { 	// check if we found the file
			fin.close();
			return;
		}

		size_t id;
		double popn;
		for (size_t i = 0; i < mod->levels.size(); i++) {
			getline(fin, str);
			str = trim(str);
			if (str.size() == 0) break;
			sfin.str(str);
			sfin >> id >> popn;
			sfin.clear();
			mod->set_level(id, "pop", popn);
		}

		fin.close();
	}
    
	void initial_solution()									// return LTE or optically thin initial solution if populations haven't been read from file, i.e. if the population of the first and last levels < 0
	{
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			double temperature = 1.e-30;
			double pops_sum = 0.0;

			if (initialSolutionSource == 0) read_pops(filename_pops_in[ispec], &mols[ispec]);
			
			partition_function_ratio[ispec] = 0.0; 	// partition function calculated for a given set of levels (the full, actual partition function should include all levels)
			for (size_t i = mols[ispec].levels.size(); i-- > 0; ) partition_function_ratio[ispec] += mols[ispec].levels[i].g * exp( -SPEED_OF_LIGHT*PLANK_CONSTANT*mols[ispec].levels[i].E / (BOLTZMANN_CONSTANT*modelPhysPars::Tks) );
			if (USE_PARTITION_FUNCTIONS_RATIO) partition_function_ratio[ispec] = min(1.0, partition_function_ratio[ispec] / full_partition_function[ispec]);
			else partition_function_ratio[ispec] = 1.0;

			switch (initialSolutionSource) {
				case 0:
					if (mols[ispec].levels[0].pop < 0.e0 || mols[ispec].levels[mols[ispec].levels.size() - 1].pop < 0.e0) throw runtime_error("file with initial solution haven't been read");
					for (size_t i = 0; i < mols[ispec].levels.size(); i++)	{
						mols[ispec].levels[i].pop = max(MIN_POP, mols[ispec].levels[i].pop);
						pops_sum += mols[ispec].levels[i].pop;
					}
					pops_sum = partition_function_ratio[ispec] / pops_sum;
					for (size_t i = 0; i < mols[ispec].levels.size(); i++) mols[ispec].levels[i].pop *= pops_sum;
					continue;
					break;
				case 1: 	// optically thin case
					temperature = 2.728; 	// cosmic microwave background temperature
					break;
				case 2: 	// LTE
					temperature = modelPhysPars::Tks;
					break;
				default:
					throw runtime_error("wrong value of initial solution parameter, check Parameters/RadiativeTransfer.txt file or console input");
					break;
			}

			double z = 0.0; 	// partition function calculated for a given set of levels (the full, actual partition function should include all levels)
			for (size_t i = mols[ispec].levels.size(); i-- > 0; ) z += mols[ispec].levels[i].g * exp( -(SPEED_OF_LIGHT*PLANK_CONSTANT/BOLTZMANN_CONSTANT) * mols[ispec].levels[i].E / temperature );
			z = partition_function_ratio[ispec] / z;

			for (size_t i = 0; i < mols[ispec].levels.size(); i++) {
				mols[ispec].levels[i].pop = max(MIN_POP, mols[ispec].levels[i].g * exp( -(SPEED_OF_LIGHT*PLANK_CONSTANT/BOLTZMANN_CONSTANT) * mols[ispec].levels[i].E / temperature ) * z);
			}
		}
	}

	double abs_coeff(const size_t & i, molModel *mol)	// absorption coefficient for i-th transition
	{
		auto kabsi = [&](const size_t & i, const size_t ispec) { 	//absorption coefficient for i-th transition of ispec molecule without taking into account line overlapping
			const size_t & up = mols[ispec].rad_trans[i].up_level;
			const size_t & low = mols[ispec].rad_trans[i].low_level;
			return modelPhysPars::Hdens * modelPhysPars::abundance[ispec] * (mols[ispec].levels[low].pop * mols[ispec].rad_trans[i].Blu - mols[ispec].levels[up].pop * mols[ispec].rad_trans[i].Bul);
		};

		double kabs = kabsi(i, mol->idspec);
		for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) { 	// take into account line overlapping
			kabs += kabsi(mol->rad_trans[i].blends[j].id, mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[i].blends[j].fac;
		}
		return (HC4PI * invlineWidth) * kabs + mol->rad_trans[i].kabs_dust;
	}

	double emiss_coeff(const size_t & i, molModel *mol) // emission coefficient for i-th transition
	{
		auto emissi = [&](const size_t & i, const size_t ispec) { 	// emission coefficient for i-th transition of ispec molecule without taking into account line overlapping
			return modelPhysPars::Hdens * modelPhysPars::abundance[ispec] * mols[ispec].rad_trans[i].A * mols[ispec].levels[mols[ispec].rad_trans[i].up_level].pop;
		};

		double emiss = emissi(i, mol->idspec);
		for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) { 	// take into account line overlapping
			emiss += emissi(mol->rad_trans[i].blends[j].id, mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[i].blends[j].fac;
		}
		return (HC4PI * invlineWidth) * emiss + mol->rad_trans[i].emiss_dust;
	}

	void compute_tau(const size_t & i, molModel *mol)		// optical depth of the maser region for radiative transition i
	{ // see e.g. Appendix A in Sobolev et al. 1997
		auto tauif = [&](const size_t & i, const size_t ispec) { 	//optical depth for i-th transition of ispec molecule without taking into account line overlapping
			const size_t & up = mols[ispec].rad_trans[i].up_level;
			const size_t & low = mols[ispec].rad_trans[i].low_level;
			return modelPhysPars::NdV[ispec] * (mols[ispec].levels[low].pop * mols[ispec].rad_trans[i].Blu - mols[ispec].levels[up].pop * mols[ispec].rad_trans[i].Bul);
		};

		mol->rad_trans[i].tau = tauif(i, mol->idspec);
		for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) { 	// take into account line overlapping
			mol->rad_trans[i].tau += tauif(mol->rad_trans[i].blends[j].id, mol->rad_trans[i].blends[j].ispec) * mol->rad_trans[i].blends[j].fac;
		}
		mol->rad_trans[i].tau *= HC4PI;
		mol->rad_trans[i].tau += mol->rad_trans[i].taud_in;
		if (mol->rad_trans[i].tau < MIN_TAU) mol->rad_trans[i].tau = MIN_TAU; 	// MIN_TAU is defined in hiddenParameters.h
	}
	
	void compute_J_S_beta(molModel *mol, const size_t & i, beta_LVG & LVG_beta, double & S, double & beta, double & beta_S)		//computes mean intensity, source function, escape probability, and their product for radiative transition i
	{ // see also equation for Jav in Appendix A of Sobolev et al. 1997
		const double emiss = emiss_coeff(i, mol);
		const double kabs = abs_coeff(i, mol);

		beta = 1.0e00; 		// escape probability = beta(tau) -> 1 for tau -> 0.0
		beta_S = 0.0e00; 		// (1 - beta) * source function;
		if (fabs(kabs) > 0.0) {
			beta = LVG_beta.beta(mol->rad_trans[i].tau);
			S = emiss / kabs;
			beta_S = S * (1.0e00 - beta);
		} else {
			S = 0.0;
			beta_S = 0.5 * emiss * modelPhysPars::NdV[mol->idspec] * lineWidth / (modelPhysPars::Hdens * modelPhysPars::abundance[mol->idspec]); 	// note, that (1-b)/tau -> 0.5 for tau -> 0.0
		}

		mol->rad_trans[i].J = beta_S + beta * mol->rad_trans[i].JExt + 
							  dust_HII_CMB_Jext_emission->HII_region_at_LOS * LVG_beta.betaHII_LOS(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII + 
							  (1 - dust_HII_CMB_Jext_emission->HII_region_at_LOS) * LVG_beta.betaHII_pump(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII; 	// sum of internal and external radiation mean intensities
	}

	void compute_Tex(const size_t & i, molModel *mol)		// computes excitation temperature for radiative transition i
	{
		const size_t & up = mol->rad_trans[i].up_level;
		const size_t & low = mol->rad_trans[i].low_level;
		const double level_pop_ratio = mol->levels[low].g*mol->levels[up].pop / (mol->levels[up].g * mol->levels[low].pop); 	// = nu*gl / (nl*gu), where nu,nl - level populations, gl,gu - statistical weights
		if (fabs(level_pop_ratio-1.) > DBL_EPSILON) mol->rad_trans[i].Tex = - (PLANK_CONSTANT/BOLTZMANN_CONSTANT)*mol->rad_trans[i].nu / log( level_pop_ratio );
		else mol->rad_trans[i].Tex = modelPhysPars::Tks;  	// if populations of two levels are equal then it is assumed that the transition between them is in LTE
	}

	double compute_source_function(const size_t & i, molModel *mol)		// computes source function
	{
		const double emiss = emiss_coeff(i, mol);
		const double kabs = abs_coeff(i, mol);

		double S = 0.0; 	// source function;
		if (fabs(kabs) > 0.0) S = emiss / kabs;
		return S;
	}

	size_t solve_eq_sys(double A[], double B[], molModel *mol)			// solves linear system of equations A*X = B with LU decomposition
	{
		const size_t & n = mol->levels.size();
		size_t *ipvt = new size_t[n];

		const size_t info = dgefa( A, n, n, ipvt ); 	// factorization of the matrix A, degfa is the LINPACK function, see linpack/linpack_d.hpp
		if (info > 0) {
			delete[] ipvt;
			return info;
		}

		dgesl( A, n, n, ipvt, B ); 						// obtaining the solution X using factorized matrix A, degsl is the LINPACK function, see linpack/linpack_d.hpp

		delete[] ipvt;
		return info;
	}

	void compute_brightness_temperature(const size_t & i, molModel *mol)		// computes brightness temperature for radiative transition i
	{ 	// it is somewhat similar to the last Equation for Tbr in Appendix A of Sobolev et al. 1997
		const double & nu = mol->rad_trans[i].nu;

		const double Tcloud_cont = exp(-mol->rad_trans[i].taud_in*beamH) * dust_HII_CMB_Jext_emission->continuum_behind_maser_region(nu) +
									oneMinusExp(mol->rad_trans[i].taud_in*beamH) * dust_HII_CMB_Jext_emission->inner_dust_source_function(nu); // contribution of the maser cloud into continuum emission
		const double Tdust_infront_cont = exp(-dust_HII_CMB_Jext_emission->tau_dust_LOS(nu)) * Tcloud_cont; // absorption of continuum by the external dust in front of the maser region at the line-of-sight

		const double Tcloud = exp(-mol->rad_trans[i].tau*beamH) * dust_HII_CMB_Jext_emission->continuum_behind_maser_region(nu) +
									oneMinusExp(mol->rad_trans[i].tau*beamH) * compute_source_function(i, mol); // contribution of the maser cloud into total emission
		const double Tdust_infront = exp(-dust_HII_CMB_Jext_emission->tau_dust_LOS(nu)) * Tcloud; // absorption by the external dust in front of the maser region at the line-of-sight

		mol->rad_trans[i].Tbr = (Tdust_infront - Tdust_infront_cont) * (pow(SPEED_OF_LIGHT/nu, 2.0) / (2. * BOLTZMANN_CONSTANT));
	}

	void find_blends() 	// searching for overlapped lines, only local overlapping is taken into account
	{
		std::function<double(const double &)> profile_shape = [this](const double & velf) -> double { return 0.0; };
		if (line_profile_shape == "r") { 	// choosing line profile shape
			profile_shape = [this](const double & velf) -> double {
				return 1. - velf / lineWidth;
			};
		} else if (line_profile_shape == "g") {
			profile_shape = [this](const double & velf) -> double {
				return 1. - exp(-velf*velf / (lineWidth*lineWidth) * 0.5);
			};
		} else {
			throw runtime_error("unknown parameter for line profile shape in find_blends function in RT.h");
		}
		double vel_fac = 0.0;
		for (size_t ispec0 = 0; ispec0 < modelPhysPars::nSpecies; ispec0++) {
			for (size_t i = 0; i < mols[ispec0].rad_trans.size(); i++) {
				for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
					for (size_t j = 0; j < mols[ispec].rad_trans.size(); j++) {
						vel_fac = fabs(mols[ispec0].rad_trans[i].nu - mols[ispec].rad_trans[j].nu) / mols[ispec0].rad_trans[i].nu * SPEED_OF_LIGHT;
						if (vel_fac < lineWidth && ((i != j && ispec0 == ispec) || (ispec0 != ispec)) )
							mols[ispec0].rad_trans[i].add_overlapped_line(ispec, j, profile_shape(vel_fac));
					}
				}
			}
		}
	}

	void Ng_acceleration(double pop[], vector <vector <double> > & oldpops, molModel *mol)
	{
		constexpr const size_t m = Ng_order + 1;

		double dnm, dnmj, dnmk; 	// population differences between different iterations

		double Cjk[Ng_order * Ng_order]; // the matrix of coefficients that will be used for Ng acceleration
		double bk[Ng_order]; 	// the vector of coefficients that will be used for Ng acceleration
		size_t ipvt[Ng_order]; 	// auxiliary vector that will be used for Ng acceleration

		for (size_t j = 0; j < Ng_order; j++) {
			for (size_t k = j; k < Ng_order; k++) {
				Cjk[j + Ng_order*k] = 0.0;
				for (size_t i = mol->levels.size(); i-- > 0; ) {
					dnm = pop[i] - oldpops[i][m];
					dnmj = oldpops[i][m-j-1] - oldpops[i][m-j-2];
					dnmk = oldpops[i][m-k-1] - oldpops[i][m-k-2];
					Cjk[j + Ng_order*k] += (dnm - dnmj) * (dnm - dnmk) / mol->levels[i].pop; 	// equation 8.113 from Gray's maser book
				}
				Cjk[k + Ng_order*j] = Cjk[j + Ng_order*k]; 	// matrix Cjk is symmetric
			}
			bk[j] = 0.0;
			for (size_t i = mol->levels.size(); i-- > 0; ) {
				dnm = pop[i] - oldpops[i][m];
				dnmj = oldpops[i][m-j-1] - oldpops[i][m-j-2];
				bk[j] += dnm * (dnm - dnmj) / mol->levels[i].pop; 	// equation 8.114 from Gray's maser book
			}
		}

		// solving equation 8.116 from Gray's maser book
		const size_t info = dgefa( Cjk, Ng_order, Ng_order, ipvt ); 	// factorization of the matrix Cjk, degfa is the LINPACK function, see linpack/linpack_d.hpp
		if (info > 0) {
			//cerr << "#warning: searching for Ng acceleration coefficients in RT.h have failed ";
			return;
		}
		dgesl( Cjk, Ng_order, Ng_order, ipvt, bk ); 	// obtaining the solution alpha using factorized matrix Cjk, degsl is the LINPACK function, see linpack/linpack_d.hpp

		double sum_b = 0.0;
		for (size_t k = 0; k < Ng_order; k++) sum_b += bk[k];
		sum_b = 1.0 - sum_b;

		for (size_t i = 0; i < mol->levels.size(); i++ ) {
			double old_pop = pop[i];
			pop[i] = 0.0;
			for (size_t k = 0; k < Ng_order; k++) {
				pop[i] += bk[k] * oldpops[i][m-k];
			}
			pop[i] = fma(sum_b, old_pop, pop[i]); 	// equation 8.107 from Gray's maser book
		}
	}

	bool update_check_pops(molModel *mol, double pop[], size_t & speciesWithMaxRPopDiff, size_t & levelWithMaxRPopDiff, double & MaxRPopDiff, const unsigned int & iter, vector <vector <double> > & oldpops_Ng, double & pop_norm)	// updates old populations and checks the populations, finds maximum relative difference of pops between iterations; computes norm of vector populations
	{
		double popRDiff;
		bool there_were_bad_levels = false;
		pop_norm = 1.e-30;
		for (size_t i = 0; i < mol->levels.size(); i++) {
			// find maximum relative difference of populations between succesive iterations
			popRDiff = fabs((pop[i] - mol->levels[i].pop) / mol->levels[i].pop);
			if (popRDiff > MaxRPopDiff && fabs(pop[i]) > MIN_POP_FOR_DIFF_CALC && i != 0 && i != mol->levels.size() - 1) { // note: the first and last levels are not taken into account because pop[0]=1-pop[1]-pop[2]-...-pop[n] and the last level population is usually low and may oscillate strongly
				MaxRPopDiff = popRDiff;
				levelWithMaxRPopDiff = i + 1; // note conversion from 0-based level indexing to 1-based
				speciesWithMaxRPopDiff = mol->idspec + 1; // note conversion from 0-based level indexing to 1-based
			}
			if (pop[i] >= MAX_POP || pop[i] <= 0.0) there_were_bad_levels = true;
		}
		if (MaxRPopDiff < 0) there_were_bad_levels = true;
		if (!there_were_bad_levels) {
			if (DoNg && iter > Ng_start) {
				vector<double> temp_pop;
				for (size_t i = 0; i < mol->levels.size(); i++) temp_pop.push_back(max(pop[i], MIN_POP));
				Ng_acceleration(pop, oldpops_Ng, mol);
				for (size_t i = 0; i < mol->levels.size(); i++) oldpops_Ng[i][Ng_order + 1] = temp_pop[i];
				temp_pop.clear();
			}
			vector<bool> underrelax_level(mol->levels.size(), false);
			double minimum_tau = 0.0;
			for (size_t i = 0; i < mol->rad_trans.size(); i++) {
				if (mol->rad_trans[i].tau < MAX_TAU_FOR_TRANSITIONS_TO_UNDERELAX) {
					underrelax_level[mol->rad_trans[i].up_level] = true;
					underrelax_level[mol->rad_trans[i].low_level] = true;
				}
				if (mol->rad_trans[i].tau < minimum_tau) minimum_tau = mol->rad_trans[i].tau;
			}
			/* more levels to underralax
			for (size_t i = 0; i < mol->rad_trans.size(); i++) {
				if (underrelax_level[mol->rad_trans[i].up_level]) underrelax_level[mol->rad_trans[i].low_level] = true;
				if (underrelax_level[mol->rad_trans[i].low_level]) underrelax_level[mol->rad_trans[i].up_level] = true;
			}
			*/
			const double var_under_relax_fac = 1.0 / (ceil(fabs(minimum_tau)));
			double pops_sum = 0;
			for (size_t i = mol->levels.size(); i-- > 0; ) {
				if (underrelax_level[i]) mol->levels[i].pop = max(pop[i], MIN_POP) * var_under_relax_fac + (1. - var_under_relax_fac) * mol->levels[i].pop;
				else mol->levels[i].pop = max(pop[i], MIN_POP);
				pops_sum += mol->levels[i].pop;
			}
			pops_sum = partition_function_ratio[mol->idspec] / pops_sum;
			for (size_t i = mol->levels.size(); i-- > 0; ) {
				mol->levels[i].pop *= pops_sum;
				pop_norm += mol->levels[i].pop * mol->levels[i].pop;
				for (size_t olp_i = 0; olp_i < (Ng_order + 1); olp_i++) oldpops_Ng[i][olp_i] = oldpops_Ng[i][olp_i + 1];
				if (!(DoNg && iter > Ng_start)) {
					oldpops_Ng[i][Ng_order + 1] = max(pop[i], MIN_POP);
				}
			}
			pop_norm = sqrt(pop_norm);
		}
		return there_were_bad_levels;
	}

	void init_some_parameters()
	{
		if (this->input_full_partition_function.size() > 1) {
			if (this->input_full_partition_function.size() != modelPhysPars::nSpecies) runtime_error("size of input_full_partition_function != 1 or to a number of molecular species in RT.h and hiddenParameters.h");
		}
		for (size_t ispec = 0; ispec < modelPhysPars::nSpecies; ispec++) {
			this->partition_function_ratio.push_back(1.0);
			this->mols.push_back(molModel(ispec));
			if (this->input_full_partition_function.size() > 1) this->full_partition_function.push_back(this->input_full_partition_function[ispec]);
			else this->full_partition_function.push_back(this->input_full_partition_function[0]);
		}
		this->lineWidth = 0.0;
		this->invlineWidth = 1.0;
		this->maxNumberOfIterations = 1;
		this->initialSolutionSource = 2;
		this->MAX_DpopsDt_EPS = 1.e-6;
		this->beamH = 1.0;
	}

public:

	vector <molModel> mols;				// data on molecular levels and transitions, see molModel.h
	vector <string> filename_lamda; 	// name of the file with molecular spectroscopic data in LAMDA format
	vector <string> filename_pops_in; 	// name of the file with initial level populations
	vector <string> filename_pops_out; 	// name of the file to store final populations
	
	virtual int radiative_transfer()	// main function which should be called to perform radiative transfer calculations
	{
		return 0;
	}

	RT()
	{
		this->init_some_parameters();
		ifstream fin;
		fin.open("Parameters/RadiativeTransfer.txt", ios::in);
		this->read_parameters(fin);
		fin.close();
		fin.open("Parameters/Dust_HII_CMB_Jext_Radiation.txt", ios::in);
		this->dust_HII_CMB_Jext_emission = new dust_HII_CMB_Jext_radiation(fin);
		fin.close();
		this->cerr_output_iter_progress = true;
	}
	
	RT(istream & cin)
	{
		this->init_some_parameters();
		this->read_parameters(cin);
		this->dust_HII_CMB_Jext_emission = new dust_HII_CMB_Jext_radiation(cin);
		this->cerr_output_iter_progress = false;
	}

	virtual ~RT()
	{
		delete this->dust_HII_CMB_Jext_emission;
	}
};
