#pragma once
#include "spline.h"

class dust_HII_CMB_Jext_radiation
{
private:
	
	unsigned short read_Jext_from_file; // = 1 - read mean intensity of external radiation from file
	vector<double> nu_J_from_file; // array of frequencies from the file with external radiation mean intensity
	vector<double> J_from_file; // array of occupation numbers from the file with external radiation

	vector<double> time_from_file; // array of time points from file with dependence of dust temperature and dillution factor on time
	vector<double> Td_from_file; // array of dust temperatures from file with dependence of dust temperature and dillution factor on time
	vector<double> Wd_from_file; // array of dust dillution factors from file with dependence of dust temperature and dillution factor on time

	// Dust emission is given by a black body emission multiplied by dilution factor; the mean intensity J = Wd * (1 - exp(tau_nu0 * (nu/nu0)^p)) * planck_function(Td,nu)
	double tau_nu0;				// optical depth at frequency=nu0
	double k_nu0;				// [cm^/g] mass dust absorption coefficient at frequency=nu0
	double nu0;					// [Hz], frequency at which the external dust layer optical depth = tau_nu0
	double p;					// spectral index
	double Td_in;				// [K], temperature of dust within the maser region
	double Td;					// [K], external dust temperature
	double Wd;					// dilution factor for the dust emission
	unsigned short inner_dust_included;	// = 0 - there will be no dust inside the maser region; = 1 - there will be dust inside the maser region
	unsigned short outer_dust_at_LOS;	// = 0 - the dust outside the maser region is not on the line of sight; = 1 - the dust outside the maser region is on the line of sight
	bool use_kabs_from_file;	// = true - take dust absorption coefficients from file; = false - assume that absorption coefficient ~ frequency^p
	string kabs_filename; 		// name of file with dust mass absorption coefficient depenedence on wavelength (the first column is wavelength in microns, the second column is mass absorption coefficient in cm^2/g)
	vector<double> log_wave_kabs_from_file; // array of logarithms of wavelenght from the file with dust absorption coefficient
	vector<double> log_kabs_from_file; // array of logarithms of dust absorptoin coefficients
	akima_spline spline_kabs; // structure for akima spline interpolation of dust absorption coefficients on wavelength
	
	// HII region emission, see e.g. Sobolev et al. 1997, J = WHii * (1-exp(tauHII)) * planck_function(Te,nu), tauHII = (HII_turn_freq/nu)^2
	double WHii;				//dilution factor for the HII region emission
	double HII_turn_freq; 		// [Hz], turnover frequency of HII region emission
	double Te; 					// [K], HII region electron temperature
	unsigned short HII_region_at_LOS_behind_maser;	// = 0 - HII region does not intersect line-of-sight behind the maser region; = 1 - HII region intersects line-of-sight behind the maser region
	unsigned short HII_region_at_LOS_infront_maser;	// = 0 - HII region does not intersect line-of-sight in front of the maser region; = 1 - HII region intersects line-of-sight in front of the maser region

	double Wcross; // dilution factor for the cross-section of HII region and external dust layer

	// Cosmic microwave background temperature, Kelvins; J = planck_function(T_CMB,nu)
	double T_CMB;

	void read_kabs_from_file(const string & filename) 	//read mean intensity of external emission from a file
	{
		ifstream fin;
		istringstream sfin;
		string str;
		double input_wavelength, input_kabs;

		fin.open(filename.c_str(), ios::in);
		if ( !fin.good() ) { 	// check if we found the file
			fin.close();
			throw runtime_error("can't find input file with dust absorption coefficients");
		}

		while (getline(fin,str)) {
			str = trim(str); 	// trim is taken from auxiliary.h
			if (str.size() != 0) {
				sfin.str(str);
				sfin >> input_wavelength >> input_kabs;
				sfin.clear();
				log_wave_kabs_from_file.push_back( log10(input_wavelength * 1.e-4) ); // microns -> cm
				log_kabs_from_file.push_back( log10(input_kabs) ); 
			}
			else break;
		}
		fin.close();
		spline_kabs.akima_init(log_wave_kabs_from_file, log_kabs_from_file); 	// creating spline for dust absorption coefficients from file
	}

	double get_dust_kabs_from_file(const double & freq)
	{
		const double loglam = log10(SPEED_OF_LIGHT / freq);
		const size_t last_id = log_kabs_from_file.size();
		if (loglam < log_wave_kabs_from_file[0]) {
			return pow( 10., log_kabs_from_file[0] + (loglam-log_wave_kabs_from_file[0]) * (log_kabs_from_file[1]-log_kabs_from_file[0])/(log_wave_kabs_from_file[1]-log_wave_kabs_from_file[0]) );
		} else if (loglam > log_wave_kabs_from_file[last_id-1]) {
			return pow( 10., log_kabs_from_file[last_id-2] + (loglam-log_wave_kabs_from_file[last_id-2]) * (log_kabs_from_file[last_id-1]-log_kabs_from_file[last_id-2])/(log_wave_kabs_from_file[last_id-1]-log_wave_kabs_from_file[last_id-2]) );
		} else return pow( 10., spline_kabs.akima_eval(loglam) );
	}

	void read_J_from_file(const string & filename) 	//read mean intensity of external emission from a file
	{
		ifstream fin;
		istringstream sfin;
		string str;
		double input_nu, input_J;

		fin.open(filename.c_str(), ios::in);
		if ( !fin.good() ) { 	// check if we found the file
			fin.close();
			throw runtime_error("can't find input file with occupation number of external emission");
		}

		getline(fin,str);
		getline(fin,str);
		while (getline(fin,str)) {
			str = trim(str); 	// trim is taken from auxiliary.h
			if (str.size() != 0) {
				sfin.str(str);
				sfin >> input_nu >> input_J;
				sfin.clear();
				nu_J_from_file.push_back( input_nu );
				J_from_file.push_back( input_J ); 
			}
			else break;
		}
		fin.close();
	}

	double interpolate_Jext_from_file(const double & freq) 	// interpolate mean intensity that has been read from file for a given frequency
	{
		if (freq <= nu_J_from_file[0]) return J_from_file[0];
		else if (freq >= nu_J_from_file[nu_J_from_file.size()-1]) return J_from_file[J_from_file.size()-1];
		else {
			for (size_t i = 0; i < nu_J_from_file.size()-1; i++ ) {
				if (nu_J_from_file[i] <= freq && freq <= nu_J_from_file[i+1])
					return J_from_file[i] + (J_from_file[i+1] - J_from_file[i]) / (nu_J_from_file[i+1] - nu_J_from_file[i]) * (freq - nu_J_from_file[i]);  
			} 
		}
		return 0.0;
	}

	void read_TdWd_from_file(const string & filename) 	//read mean intensity of external emission from a file
	{
		ifstream fin;
		istringstream sfin;
		string str;
		double input_time, input_Td, input_Wd;

		fin.open(filename.c_str(), ios::in);
		if ( !fin.good() ) { 	// check if we found the file
			fin.close();
			throw runtime_error("can't find input file with dependence of dust temperature on time");
		}

		while (getline(fin,str)) {
			str = trim(str); 	// trim is taken from auxiliary.h
			if (str.size() != 0) {
				sfin.str(str);
				sfin >> input_time >> input_Td >> input_Wd;
				sfin.clear();
				time_from_file.push_back( input_time );
				Td_from_file.push_back( input_Td );
				Wd_from_file.push_back( input_Wd );
			}
			else break;
		}
		fin.close();
	}

	double interpolate_TdWd_from_file(const double & time, const vector <double> & itime_from_file, const vector <double> & Xd_from_file) 	// interpolate dust temperature that has been read from file for a time
	{
		if (time < itime_from_file[0]) return Xd_from_file[0];
		if (time > itime_from_file[itime_from_file.size()-1]) return Xd_from_file[itime_from_file.size()-1];
		for (size_t i = 0; i < itime_from_file.size()-1; i++ ) {
			if (itime_from_file[i] <= time && time <= itime_from_file[i+1])
				return Xd_from_file[i] + (Xd_from_file[i+1] - Xd_from_file[i]) / (itime_from_file[i+1] - itime_from_file[i]) * (time - itime_from_file[i]);
		}
		return 0.0;
	}
	
	template <typename T>
	void read_parameters(T & fin)
	{
		istringstream sfin;
		string str;

		if ( !fin.good() ) 	// check if we found the file or the console input is not corrupted
			throw runtime_error("can't find input Parameters/Dust_HII_CMB_Jext_Radiation.txt file or command line input is corrupted");

		getline(fin, str);
		read_Jext_from_file = readline<unsigned short>(fin); 	// readline function is taken from auxiliary.h
		if (!(read_Jext_from_file == 0 || read_Jext_from_file == 1))
			throw runtime_error("the first parameter in Parameters/Dust_HII_CMB_Jext_Radiation.txt file or command line input should be 0 or 1");

		getline(fin, str);
		getline(fin, str);
		if (read_Jext_from_file == 1) {
			read_J_from_file(trim(str));
		}

		for (short i = 0; i < 4; i++) getline(fin, str);
		sfin.str(trim(str));
		sfin >> Wd >> tau_nu0 >> nu0 >> p >> Td >> inner_dust_included;
		sfin.clear();

		for (short i = 0; i < 2; i++) getline(fin, str);
		sfin.str(trim(str));
		sfin >> Td_in >> k_nu0;
		sfin.clear();

		getline(fin, str);
		use_kabs_from_file = readline<bool>(fin);

		getline(fin, str);
		kabs_filename= readline<string>(fin);

		for (short i = 0; i < 4; i++) getline(fin, str);
		sfin.str(trim(str));
		sfin >> WHii >> HII_turn_freq >> Te >> HII_region_at_LOS;
		sfin.clear();

		getline(fin, str);
		T_CMB = readline<double>(fin);
	}

	void check_sources_at_LOS(const double & iWd, const double & iWhii) // check whether the external radiation sources are at line-of-sight
	{
		HII_region_at_LOS_behind_maser = 0;
		HII_region_at_LOS_infront_maser = 0;
		outer_dust_at_LOS = 1;
		if (iWd   <= 0.5) outer_dust_at_LOS = 0;
		if (iWhii > DBL_EPSILON && HII_region_at_LOS == 1) HII_region_at_LOS_behind_maser = 1;
		if (iWhii > 1. - DBL_EPSILON && HII_region_at_LOS == 1) HII_region_at_LOS_infront_maser = 1;
		if (HII_region_at_LOS == 0 && iWhii > 0.5) {
			HII_region_at_LOS_behind_maser = 1;
			HII_region_at_LOS_infront_maser = 1;
		}
	}

	double tau_HII(const double freq) 	// optical depth of the HII region at a given frequency
	{
		const double freq_ratio = HII_turn_freq / freq;
		return freq_ratio * freq_ratio; 	// see Sobolev & Deguchi 1994, Sobolev et al. 1997
	}

	double var_Td(const double & time) // dependence of external dust temperature on time
	{
		return Td;
		/*return interpolate_TdWd_from_file(time, time_from_file, Td_from_file);
		const double period = 3*3600.;
		if (time > 3 * period) return Td;
		return Td + 3 * sin(2 * PI * time / period);*/
	}

	double var_Wd(const double & time) // dependence of dust dillution factor on time
	{
		return Wd;
		/*return interpolate_TdWd_from_file(time, time_from_file, Wd_from_file);
		const double period = 3*3600.;
		if (time > 3 * period) return Wd;
		return Wd + 0.1 * sin(2 * PI * time / period);*/
	}

	double var_Td_in(const double & time) // dependence of internal dust temperature on time
	{
		return Td_in;
	}

	double get_Wcross(const double & iWd, const double & iWhii)  //calculates Wcross
	{
		if (HII_region_at_LOS == 0) return min(iWd, iWhii);
		if (iWd + iWhii <= (1./8.)) return 0;
		const double sqrtWmin = sqrt(min(iWd, iWhii));
		const double sqrtWmax = sqrt(max(iWd, iWhii));
		return 0.5 * sqrtWmin * (sqrtWmax + sqrtWmin - 0.5 + max(0., sqrtWmin - 0.5));
	}

public:

	unsigned short HII_region_at_LOS;	// = 0 - HII region is not centered on the line of sight; = 1 - HII region is centered on the line of sight

	double tau_HII_infront(const double freq) 	// optical depth of the HII region at a given frequency on line-of-sight in front of the maser region
	{
		return HII_region_at_LOS_infront_maser * tau_HII(freq); 	// see Sobolev & Deguchi 1994, Sobolev et al. 1997
	}

	double tau_dust(const double & freq) 	// optical depth of external dust at a given frequency
	{
		if (use_kabs_from_file) return tau_nu0 * get_dust_kabs_from_file(freq) / get_dust_kabs_from_file(nu0);
		else return tau_nu0 * pow( freq / nu0, p); // see e.g. Sobolev et al. 1997
	}

	double tau_dust_LOS(const double & freq) 	// optical depth of external dust at a given frequency on the line of sight
	{
		return outer_dust_at_LOS * tau_dust(freq);
	}

	double tau_dust_LOS(const double & freq, const double & time) 	// optical depth of external dust at a given frequency on the line of sight depending on time
	{
		check_sources_at_LOS(var_Wd(time), WHii);
		return outer_dust_at_LOS * tau_dust(freq);
	}

	double tau_dust_in(const double & freq, const double & lineWidth) 	// optical depth of the dust inside the maser region at a given frequency and for a given line width
	{
		if (use_kabs_from_file) {
			return inner_dust_included * get_dust_kabs_from_file(freq) * AMU * modelPhysPars::dust_gas_mass_ratio * modelPhysPars::mean_mol_weight * lineWidth * modelPhysPars::max_NH2dV;
		} else {
			//return inner_dust_included * 2.6e-25 * lineWidth * modelPhysPars::max_NH2dV * pow( freq / nu0, p); // based on Sherwood et al. 1980;
			return inner_dust_included * k_nu0 * AMU * modelPhysPars::dust_gas_mass_ratio * modelPhysPars::mean_mol_weight * lineWidth * modelPhysPars::max_NH2dV * pow( freq / nu0, p);
		}
	}

	double outer_dust_source_function(const double & freq) 	// source function for the external dust emission
	{
		return planck_function(Td, freq);
	}

	double outer_dust_source_function(const double & freq, const double & time) 	// source function for the external dust emission
	{
		return planck_function(var_Td(time), freq);
	}

	double inner_dust_source_function(const double & freq) 	// source function for the dust emission inside the maser region
	{
		return inner_dust_included * planck_function(Td_in, freq);
	}

	double inner_dust_source_function(const double & freq, const double & time) 	// source function for the dust emission inside the maser region
	{
		return inner_dust_included * planck_function(var_Td_in(time), freq);
	}
	
	double compute_Jext_dust_CMB_file(const double & freq)		// returns mean intensity taken from file or sum of dust and cosmic microwave background (CMB) mean intensities
	{
		if (read_Jext_from_file == 0) {
			const double J_CMB = Wcross * exp( -(tau_dust(freq) + tau_HII(freq)) ) + \
								(Wd - Wcross) * exp(-tau_dust(freq)) + \
								(WHii - Wcross) * exp(-tau_HII(freq)) + \
								(1 - (Wd + WHii) + Wcross); // cosmic microwave background emission taking into account an absorption by the external dust layer and HII region

			const double J_dust = Wd *  oneMinusExp(tau_dust(freq))  * outer_dust_source_function(freq); // Note that this a general formula for dust emission (e.g. Sobolev et al. 1997), while eq. 4 in van der Walt 2014 gives optically thin case

			return J_dust + J_CMB * planck_function(T_CMB,freq);
		} else {
			// !!!!!! NOTE !!!!!! that this code has been used for calculations with the mean intensity provided by CLOUDY code;
			// CLOUDY outputs photon occupation number which is converted into the mean intensity below
			return interpolate_Jext_from_file(freq) * (2. * PLANK_CONSTANT * freq) * (freq / SPEED_OF_LIGHT) * (freq / SPEED_OF_LIGHT); // photon occupation number into mean intensity;
		}
	}

	double compute_Jext_dust_CMB_file(const double & freq, const double & time)		// returns mean intensity taken from file or sum of dust and cosmic microwave background (CMB) mean intensities
	{
		if (read_Jext_from_file == 0) {
			const double var_Wcross = get_Wcross(var_Wd(time), WHii);
			const double J_CMB = var_Wcross * exp( -(tau_dust(freq) + tau_HII(freq)) ) + \
								(var_Wd(time) - var_Wcross) * exp(-tau_dust(freq)) + \
								(WHii - var_Wcross) * exp(-tau_HII(freq)) + \
								(1 - (var_Wd(time) + WHii) + var_Wcross); // cosmic microwave background emission taking into account an absorption by the external dust layer and HII region

			const double J_dust = var_Wd(time) *  oneMinusExp(tau_dust(freq))  * outer_dust_source_function(freq, time); // Note that this a general formula for dust emission (e.g. Sobolev et al. 1997), while eq. 4 in van der Walt 2014 gives optically thin case

			return J_dust + J_CMB * planck_function(T_CMB,freq);
		} else {
			// !!!!!! NOTE !!!!!! that this code has been used for calculations with the mean intensity provided by CLOUDY code;
			// CLOUDY outputs photon occupation number which is converted into the mean intensity below
			return interpolate_Jext_from_file(freq) * (2. * PLANK_CONSTANT * freq) * (freq / SPEED_OF_LIGHT) * (freq / SPEED_OF_LIGHT); // photon occupation number into mean intensity;
		}
	}

	double compute_JextHII(const double & freq)		// mean intensity of the HII region emission, needs to be separated from other types of external emission because of beaming
	{
		if (read_Jext_from_file == 0)
			return (Wcross * exp(-tau_dust(freq)) + (WHii - Wcross)) * oneMinusExp(tau_HII(freq)) * planck_function(Te,freq);
		else
			return 0.0;
	}

	double compute_JextHII(const double & freq, const double & time)		// mean intensity of the HII region emission, needs to be separated from other types of external emission because of beaming
	{
		if (read_Jext_from_file == 0) {
			const double var_Wcross = get_Wcross(var_Wd(time), WHii);
			return (var_Wcross * exp(-tau_dust(freq)) + (WHii - var_Wcross)) * oneMinusExp(tau_HII(freq)) * planck_function(Te,freq);
		}
		else
			return 0.0;
	}

	double external_dust_layer_emission(const double & freq) //emission of exteranl dust layer on LOS
	{
		return oneMinusExp(tau_dust_LOS(freq)) * outer_dust_source_function(freq);
	}

	double external_dust_layer_emission(const double & freq, const double & time) //emission of exteranl dust layer on LOS
	{
		return oneMinusExp(tau_dust_LOS(freq, time)) * outer_dust_source_function(freq, time);
	}

	double continuum_behind_maser_region(const double & freq) 	// computes continuum behind the maser region
	{
		const double THII = exp(- HII_region_at_LOS_behind_maser * tau_HII(freq)) * planck_function(T_CMB,freq) + oneMinusExp(HII_region_at_LOS_behind_maser * tau_HII(freq)) * planck_function(Te,freq);  // the HII region background emission and CMB emission absorbed by the HII region. This evaluates to non-attenuated CMB emission if the HII region is not at the line-of-sight
		const double Td1 = exp(-tau_dust_LOS(freq)) * THII + external_dust_layer_emission(freq); //the contribution of external dust behind the maser region in the case if the external dust is at the line-of-sight
		return Td1;
	}

	double continuum_behind_maser_region(const double & freq, const double & time) 	// computes continuum behind the maser region
	{
		const double THII = exp(- HII_region_at_LOS_behind_maser * tau_HII(freq)) * planck_function(T_CMB,freq) + oneMinusExp(HII_region_at_LOS_behind_maser * tau_HII(freq)) * planck_function(Te,freq);  // the HII region background emission and CMB emission absorbed by the HII region. This evaluates to non-attenuated CMB emission if the HII region is not at the line-of-sight
		const double Td1 = exp(-tau_dust_LOS(freq, time)) * THII + external_dust_layer_emission(freq, time); //the contribution of external dust behind the maser region in the case if the external dust is at the line-of-sight
		return Td1;
	}

	template <typename T>
	dust_HII_CMB_Jext_radiation(T & fin)
	{
		T_CMB = 2.728;
		Te = 8000.0;
		WHii = 0;
		HII_turn_freq = 1.e8;
		Td = T_CMB;
		Td_in = Td;
		Wd = 0;
		p = 2.0;
		nu0 = 300.e9;
		tau_nu0 = 0.0;
		k_nu0 = 0.0;
		read_Jext_from_file = 0;
		inner_dust_included = 0;
		HII_region_at_LOS = 1;
		use_kabs_from_file = false;
		read_parameters(fin);
		if (Wd > 0.99999) {
			Wd = 1.0e00;
		}
		check_sources_at_LOS(Wd, WHii);
		Wcross = get_Wcross(Wd, WHii);
		if (use_kabs_from_file) read_kabs_from_file(kabs_filename);
		//read_TdWd_from_file("TdWd.txt");
	}

	~dust_HII_CMB_Jext_radiation()
	{
		nu_J_from_file.clear();
		J_from_file.clear();
		time_from_file.clear();
		Td_from_file.clear();
		Wd_from_file.clear();
		log_wave_kabs_from_file.clear();
		log_kabs_from_file.clear();
		spline_kabs.clear();
	}
};
