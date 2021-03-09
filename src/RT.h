#pragma once
#include "Dust_HII_CMB_Jext_Radiation.h"
#include "beta_LVG.h"
#include "linpack/linpack_d.hpp"

class RT		// base class containing functions for radiative transfer calculations;
{
protected:

	unsigned short initialSolutionSource;						// =0 - initial solution is taken from file ; =1 - initial solution is computed for optically thin case ; =2 - LTE initial solution
	double MAX_DpopsDt_EPS;										// the computations will stop if the maximum length of the vector with components of Dn/Dt (derivative of populations with time)
	unsigned long maxNumberOfIterations;						// maximum number of iterations
	double beamH; 												// beaming factor; this is eps^-1 = D(ln r)/D(ln V) quantity from e.g. Sobolev et al. 1997, Cragg et al. 2005
	double lineWidth; 											// spectral width of the lines; this will be used to determine blended lines
	dust_HII_CMB_Jext_radiation *dust_HII_CMB_Jext_emission; 	// this object will be used for calculations of external emission
	bool cerr_output_iter_progress; 							// = true - the current iteration number, maximum relative pop difference and corresponding level number will be printed on the standard cerr pipe
	string line_profile_shape; 									// ="g" - Gaussian line profile; ="r" - rectangular line profile
	const double full_partition_function = PARTITION_FUNCTION;	// full partition rotational function (see hiddenParameters.h)

	double partition_function_ratio;
	
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
		getline(fin, str);
		filename_pops_in = trim(str);	// trim function is taken from auxiliary.h

		getline(fin, str);
		getline(fin, str);
		filename_pops_out = trim(str);

		getline(fin, str);
		getline(fin, str);
		filename_lamda = trim(str);

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
	}
    
	void initial_solution()									// return LTE or optically thin initial solution if populations haven't been read from file, i.e. if the population of the first and last levels < 0
	{
		double temperature = 1.e-30;
		double pops_sum = 0.0;

		partition_function_ratio = 0.0; 	// partition function calculated for a given set of levels (the full, actual partition function should include all levels)
		for (size_t i = mol->levels.size(); i-- > 0; ) partition_function_ratio += mol->levels[i].g * exp( -SPEED_OF_LIGHT*PLANK_CONSTANT*mol->levels[i].E / (BOLTZMANN_CONSTANT*modelPhysPars::Tks) );
		if (USE_PARTITION_FUNCTIONS_RATIO) partition_function_ratio = min(1.0, partition_function_ratio / full_partition_function);
		else partition_function_ratio = 1.0;

		switch (initialSolutionSource) {
			case 0:
				if (mol->levels[0].pop < 0.e0 || mol->levels[mol->levels.size() - 1].pop < 0.e0) throw runtime_error("file with initial solution haven't been read");
				for (size_t i = 0; i < mol->levels.size(); i++)	{
					mol->levels[i].pop = max(MIN_POP, mol->levels[i].pop);
					pops_sum += mol->levels[i].pop;
				}
				pops_sum = partition_function_ratio / pops_sum;
				for (size_t i = 0; i < mol->levels.size(); i++)	mol->levels[i].pop *= pops_sum;
				return;
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
		for (size_t i = mol->levels.size(); i-- > 0; ) z += mol->levels[i].g * exp( -(SPEED_OF_LIGHT*PLANK_CONSTANT/BOLTZMANN_CONSTANT) * mol->levels[i].E / temperature );
		z = partition_function_ratio / z;

		for (size_t i = 0; i < mol->levels.size(); i++) {
			mol->levels[i].pop = max(MIN_POP, mol->levels[i].g * exp( -(SPEED_OF_LIGHT*PLANK_CONSTANT/BOLTZMANN_CONSTANT) * mol->levels[i].E / temperature ) * z);
		}
	}

	double abs_coeff(const size_t & i)	// absorption coefficient for i-th transition
	{
		auto kabsi = [&](const size_t & i) { 	//absorption coefficient for i-th transition without taking into account line overlapping
			const size_t & up = mol->rad_trans[i].up_level;
			const size_t & low = mol->rad_trans[i].low_level;
			//return mol->levels[low].pop * mol->rad_trans[i].Blu - mol->levels[up].pop * mol->rad_trans[i].Bul;
			const __float128 a = (__float128)(mol->levels[low].pop * mol->rad_trans[i].Blu);
			const __float128 b = (__float128)(mol->levels[up].pop * mol->rad_trans[i].Bul);
			const __float128 c = a - b;
			return (double)c;
		};

		double kabs = kabsi(i);
		for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) { 	// take into account line overlapping
			kabs += kabsi(mol->rad_trans[i].blends[j].id) * mol->rad_trans[i].blends[j].fac;
		}
		return kabs + mol->rad_trans[i].kabs_dust;
	}

	double emiss_coeff(const size_t & i) // emission coefficient for i-th transition
	{
		auto emissi = [&](const size_t & i) { 	// emission coefficient for i-th transition without taking into account line overlapping
			return mol->rad_trans[i].A * mol->levels[mol->rad_trans[i].up_level].pop;
		};

		double emiss = emissi(i);
		for (size_t j = 0; j < mol->rad_trans[i].blends.size(); j++) { 	// take into account line overlapping
			emiss += emissi(mol->rad_trans[i].blends[j].id) * mol->rad_trans[i].blends[j].fac;
		}
		return emiss + mol->rad_trans[i].emiss_dust;
	}

	void compute_tau(const size_t & i)		// optical depth of the maser region for radiative transition i
	{ // see e.g. Appendix A in Sobolev et al. 1997
		mol->rad_trans[i].tau = (PLANK_CONSTANT * SPEED_OF_LIGHT / 4. / PI) * modelPhysPars::NdV * abs_coeff(i) ;
		if (mol->rad_trans[i].tau < MIN_TAU) mol->rad_trans[i].tau = MIN_TAU; 	// MIN_TAU is defined in hiddenParameters.h
	}
	
	void compute_J_S_beta(const size_t & i, beta_LVG & LVG_beta, double & S, double & beta, double & beta_S)		//computes mean intensity, source function, escape probability, and their product for radiative transition i
	{ // see also equation for Jav in Appendix A of Sobolev et al. 1997
		const double emiss = emiss_coeff(i);
		const double kabs = abs_coeff(i);

		beta = 1.0e00; 		// escape probability = beta(tau) -> 1 for tau -> 0.0
		beta_S = 0.0e00; 		// (1 - beta) * source function;
		if (fabs(kabs) > 0.0) {
			beta = LVG_beta.beta(mol->rad_trans[i].tau);
			S = emiss / kabs;
			beta_S = S * (1.0e00 - beta);
		} else {
			S = 0.0;
			beta_S = 0.5 * emiss * modelPhysPars::NdV / (4.*PI/SPEED_OF_LIGHT/PLANK_CONSTANT); 	// note, that (1-b)/tau -> 0.5 for tau -> 0.0
		}

		mol->rad_trans[i].J = beta_S + beta * mol->rad_trans[i].JExt + 
							  dust_HII_CMB_Jext_emission->HII_region_at_LOS * LVG_beta.betaHII_LOS(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII + 
							  (1 - dust_HII_CMB_Jext_emission->HII_region_at_LOS) * LVG_beta.betaHII_pump(mol->rad_trans[i].tau, beamH) * mol->rad_trans[i].JExtHII; 	// sum of internal and external radiation mean intensities
	}

	void compute_Tex(const size_t & i)		// computes excitation temperature for radiative transition i
	{
		const size_t & up = mol->rad_trans[i].up_level;
		const size_t & low = mol->rad_trans[i].low_level;
		const double level_pop_ratio = mol->levels[low].g*mol->levels[up].pop / (mol->levels[up].g * mol->levels[low].pop); 	// = nu*gl / (nl*gu), where nu,nl - level populations, gl,gu - statistical weights
		if (fabs(level_pop_ratio-1.) > DBL_EPSILON) mol->rad_trans[i].Tex = - (PLANK_CONSTANT/BOLTZMANN_CONSTANT)*mol->rad_trans[i].nu / log( level_pop_ratio );
		else mol->rad_trans[i].Tex = modelPhysPars::Tks;  	// if populations of two levels are equal then it is assumed that the transition between them is in LTE
	}

	double compute_source_function(const size_t & i)		// computes source function
	{
		const double emiss = emiss_coeff(i);
		const double kabs = abs_coeff(i);

		double S = 0.0; 	// source function;
		if (fabs(kabs) > 0.0) S = emiss / kabs;
		return S;
	}

	size_t solve_eq_sys(double A[], double B[])			// solves linear system of equations A*X = B with LU decomposition
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

	void compute_brightness_temperature(const size_t & i)		// computes brightness temperature and intensity for radiative transition i
	{ 	// it is similar to the last Equation for Tbr in Appendix A of Sobolev et al. 1997
		const double & nu = mol->rad_trans[i].nu;
		mol->rad_trans[i].Tbr = exp(-dust_HII_CMB_Jext_emission->tau_dust_LOS(nu)) * (
			oneMinusExp(mol->rad_trans[i].tau*beamH) * compute_source_function(i) +
			(exp(-mol->rad_trans[i].tau*beamH) - exp(-mol->rad_trans[i].taud_in*beamH)) * dust_HII_CMB_Jext_emission->continuum(nu) -
			oneMinusExp(mol->rad_trans[i].taud_in*beamH) * dust_HII_CMB_Jext_emission->inner_dust_source_function(nu)
		) * (pow(SPEED_OF_LIGHT/nu, 2.0) / (2. * BOLTZMANN_CONSTANT));
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
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			for (size_t j = 0; j < mol->rad_trans.size(); j++) {
				if (i != j) {
					vel_fac = fabs(mol->rad_trans[i].nu - mol->rad_trans[j].nu) / mol->rad_trans[i].nu * SPEED_OF_LIGHT;
					if (vel_fac < lineWidth)
						mol->rad_trans[i].add_overlapped_line(j, profile_shape(vel_fac));
				}
			}
		}
	}

	void Ng_acceleration(const double pop[], vector <vector <double> > & oldpops, double & pop_norm)
	{
		constexpr const size_t m = Ng_order + 1;

		//double dnm, dnmj, dnmk; 	// population differences between different iterations
		double err_c = 0.0;

		double Cjk[Ng_order * Ng_order]; // the matrix of coefficients that will be used for Ng acceleration
		double bk[Ng_order]; 	// the vector of coefficients that will be used for Ng acceleration
		size_t ipvt[Ng_order]; 	// auxiliary vector that will be used for Ng acceleration

		for (size_t j = 0; j < Ng_order; j++) {
			for (size_t k = j; k < Ng_order; k++) {
				Cjk[j + Ng_order*k] = 0.0;
				err_c = 0.0;
				for (size_t i = mol->levels.size(); i-- > 0; ) {
					//dnm = pop[i] - oldpops[i][m];
					//dnmj = oldpops[i][m-j-1] - oldpops[i][m-j-2];
					//dnmk = oldpops[i][m-k-1] - oldpops[i][m-k-2];
					//cascade_summation(err_c, Cjk[j + Ng_order*k], (dnm - dnmj) * (dnm - dnmk)); 	// equation 8.113 from Gray's maser book
					const __float128 a1 = pop[i] + oldpops[i][m-j-2];
					const __float128 b1 = oldpops[i][m] + oldpops[i][m-j-1];
					const __float128 a2 = pop[i] + oldpops[i][m-k-2];
					const __float128 b2 = oldpops[i][m] + oldpops[i][m-k-1];
					cascade_summation(err_c, Cjk[j + Ng_order*k],
						(double)( (a1 - b1) * (a2 - b2) )
					); 	// equation 8.113 from Gray's maser book
				}
				Cjk[j + Ng_order*k] += err_c;
				Cjk[k + Ng_order*j] = Cjk[j + Ng_order*k]; 	// matrix Cjk is symmetric
			}
			bk[j] = 0.0;
			err_c = 0.0;
			for (size_t i = mol->levels.size(); i-- > 0; ) {
				//dnm = pop[i] - oldpops[i][m];
				//dnmj = oldpops[i][m-j-1] - oldpops[i][m-j-2];
				//cascade_summation(err_c, bk[j], dnm * (dnm - dnmj)); 	// equation 8.114 from Gray's maser book
				const __float128 a1 = pop[i] + oldpops[i][m-j-2];
				const __float128 b1 = oldpops[i][m] + oldpops[i][m-j-1];
				const __float128 a2 = pop[i];
				const __float128 b2 = oldpops[i][m];
				cascade_summation(err_c, bk[j],
					(double)( (a2 - b2) * (a1 - b1) )
				); 	// equation 8.114 from Gray's maser book
			}
			bk[j] += err_c;
		}

		// solving equation 8.116 from Gray's maser book
		const size_t info = dgefa( Cjk, Ng_order, Ng_order, ipvt ); 	// factorization of the matrix Cjk, degfa is the LINPACK function, see linpack/linpack_d.hpp
		if (info > 0) {
			for (size_t i = mol->levels.size(); i-- > 0; ) {
				for (size_t olp_i = 0; olp_i < m; olp_i++) oldpops[i][olp_i] = oldpops[i][olp_i + 1];
				oldpops[i][m] = pop[i];
				mol->levels[i].pop = max(pop[i], MIN_POP);
				pop_norm += mol->levels[i].pop * mol->levels[i].pop;
			}
			//cerr << "#warning: searching for Ng acceleration coefficients in RT.h have failed ";
			return;
		}
		dgesl( Cjk, Ng_order, Ng_order, ipvt, bk ); 	// obtaining the solution alpha using factorized matrix Cjk, degsl is the LINPACK function, see linpack/linpack_d.hpp

		double sum_b = 0.0;
		err_c = 0.0;
		for (size_t k = 0; k < Ng_order; k++) cascade_summation(err_c, sum_b, bk[k]);
		sum_b = 1.0 - (sum_b + err_c);

		double pops_sum = 0.0, pops_sum_err = 0.0;
		for (size_t i = 0; i < mol->levels.size(); i++ ) {
			mol->levels[i].pop = 0.0;
			err_c = 0.0;
			for (size_t k = 0; k < Ng_order; k++) {
				cascade_summation(err_c, mol->levels[i].pop, bk[k] * oldpops[i][m-k]);
			}
			mol->levels[i].pop += err_c;
			mol->levels[i].pop += sum_b * pop[i]; 	// equation 8.107 from Gray's maser book
			mol->levels[i].pop = max(MIN_POP, mol->levels[i].pop);
			cascade_summation(pops_sum_err, pops_sum, mol->levels[i].pop);
			for (size_t olp_i = 0; olp_i < m; olp_i++) oldpops[i][olp_i] = oldpops[i][olp_i + 1];
			oldpops[i][m] = pop[i];
		}
		pops_sum = partition_function_ratio / (pops_sum + pops_sum_err);
		for (size_t i = 0; i < mol->levels.size(); i++ ) {
			mol->levels[i].pop *= pops_sum;
			pop_norm += mol->levels[i].pop * mol->levels[i].pop;
		}
	}

	bool update_check_pops(const double pop[], size_t & levelWithMaxRPopDiff, double & MaxRPopDiff, const unsigned int & iter, vector <vector <double> > & oldpops_Ng, double & pop_norm, const double & underRelaxFac)	// updates old populations and checks the populations, finds maximum relative difference of pops between iterations; computes norm of vector populations
	{
		double popRDiff;
		bool there_were_bad_levels = false;
		levelWithMaxRPopDiff = 0;
		MaxRPopDiff = -1.0;
		for (size_t i = 0; i < mol->levels.size(); i++) {
			// find maximum relative difference of populations between succesive iterations
			popRDiff = fabs((pop[i] - mol->levels[i].pop) / mol->levels[i].pop);
			if (popRDiff > MaxRPopDiff && fabs(pop[i]) > MIN_POP_FOR_DIFF_CALC && i != 0 && i != mol->levels.size() - 1) { // note: the first and last levels are not taken into account because pop[0]=1-pop[1]-pop[2]-...-pop[n] and the last level population is usually low and may oscillate strongly
				MaxRPopDiff = popRDiff;
				levelWithMaxRPopDiff = i;
			}
			if (pop[i] >= MAX_POP || pop[i] <= 0.0) there_were_bad_levels = true;
		}        
		levelWithMaxRPopDiff += 1; // conversion from 0-based level indexing to 1-based
		pop_norm = 1.e-30;
		if (MaxRPopDiff < 0) there_were_bad_levels = true;
		if (!there_were_bad_levels) {
			if (DoNg && iter > Ng_start && iter % Ng_step == 0) {
				Ng_acceleration(pop, oldpops_Ng, pop_norm);
			} else {
				vector<size_t> levels_to_relax;
				for (size_t i = 0; i < mol->rad_trans.size(); i++) {
					if (mol->rad_trans[i].tau < MAX_TAU_FOR_TRANSITIONS_TO_UNDERELAX) {
						levels_to_relax.push_back(mol->rad_trans[i].up_level);
						levels_to_relax.push_back(mol->rad_trans[i].low_level);
					}
				}
				const double var_under_relax_fac = underRelaxFac / (iter % UNDER_RELAX_FAC_PERIOD + 1.);
				double pops_sum = 0;
				for (size_t i = mol->levels.size(); i-- > 0; ) {
					bool lev_need_to_be_relaxed = false;
					for (size_t j = 0; j < levels_to_relax.size(); j++) {
						if (i == levels_to_relax[j]) {
							lev_need_to_be_relaxed = true;
							break;
						}
					}
					if (lev_need_to_be_relaxed) mol->levels[i].pop = max(pop[i], MIN_POP) * var_under_relax_fac + (1. - var_under_relax_fac) * mol->levels[i].pop;
					else mol->levels[i].pop = max(pop[i], MIN_POP);
					pops_sum += mol->levels[i].pop;
				}
				levels_to_relax.clear();
				pops_sum = partition_function_ratio / pops_sum;
				for (size_t i = mol->levels.size(); i-- > 0; ) {
					mol->levels[i].pop *= pops_sum;
					pop_norm += mol->levels[i].pop * mol->levels[i].pop;
					for (size_t olp_i = 0; olp_i < (Ng_order + 1); olp_i++) oldpops_Ng[i][olp_i] = oldpops_Ng[i][olp_i + 1];
					oldpops_Ng[i][Ng_order + 1] = mol->levels[i].pop;
				}
			}
		}
		return there_were_bad_levels;
	}

public:

	molModel *mol;			// data on molecular levels and transitions, see molModel.h
	string filename_lamda; 		// name of the file with molecular spectroscopic data in LAMDA format
	string filename_pops_in; 	// name of the file with initial level populations
	string filename_pops_out; 	// name of the file to store final populations
	
	virtual int radiative_transfer()		// main function which should be called to perform radiative transfer calculations
	{
		return 0;
	}

	RT() noexcept
	{
		this->mol = new molModel();
		this->lineWidth = 0.0;
		this->maxNumberOfIterations = 1;
		this->initialSolutionSource = 2;
		this->MAX_DpopsDt_EPS = 1.e-6;
		this->beamH = 1.0;
		ifstream fin;
		fin.open("Parameters/RadiativeTransfer.txt", ios::in);
		read_parameters(fin);
		fin.close();
		fin.open("Parameters/Dust_HII_CMB_Jext_Radiation.txt", ios::in);
		this->dust_HII_CMB_Jext_emission = new dust_HII_CMB_Jext_radiation(fin);
		fin.close();
		this->cerr_output_iter_progress = true;
		this->partition_function_ratio = 1.0;
	}
	
	RT(istream & cin)
	{
		this->mol = new molModel();
		this->maxNumberOfIterations = 1;
		this->initialSolutionSource = 2;
		this->MAX_DpopsDt_EPS = 1.e-6;
		this->lineWidth = 0.0;
		this->beamH = 1.0;
		this->read_parameters(cin);
		this->dust_HII_CMB_Jext_emission = new dust_HII_CMB_Jext_radiation(cin);
		this->cerr_output_iter_progress = false;
		this->partition_function_ratio = 1.0;
	}

	RT(const unsigned short & initialSolutionSource, const double & MAX_POPS_EPS, const unsigned long & maxNumberOfIterations, const double & beamH, const double & lineWidth)
	{
		this->initialSolutionSource = initialSolutionSource;
		this->MAX_DpopsDt_EPS = MAX_POPS_EPS;
		this->maxNumberOfIterations = maxNumberOfIterations;
		this->beamH = beamH;
		this->lineWidth = lineWidth;
		this->line_profile_shape = "r";
		this->mol = new molModel();
		this->partition_function_ratio = 1.0;
		this->cerr_output_iter_progress = true;
		ifstream fin;
		fin.open("Parameters/Dust_HII_CMB_Jext_Radiation.txt", ios::in);
		this->dust_HII_CMB_Jext_emission = new dust_HII_CMB_Jext_radiation(fin);
		fin.close();
	}

	virtual ~RT()
	{
		delete dust_HII_CMB_Jext_emission;
		delete mol;
	}
};
