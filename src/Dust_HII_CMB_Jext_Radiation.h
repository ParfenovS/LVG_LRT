#pragma once

class dust_HII_CMB_Jext_radiation
{
private:
	
	unsigned short read_Jext_from_file; // = 1 - read mean intensity of external radiation from file
	vector<double> nu_from_file; // array of frequencies from the file with external radiation
	vector<double> J_from_file; // array of occupation numbers from the file with external radiation

	vector<double> time_from_file_small; // array of time points from file with dependence of dust temperature on time
	vector<double> Td_from_file_small; // array of dust temperatures from file with dependence of dust temperature on time
	vector<double> time_from_file_large; // array of time points from file with dependence of dust temperature on time
	vector<double> Td_from_file_large; // array of dust temperatures from file with dependence of dust temperature on time

	// Dust emission is given by modified black body emission multiplied by Wd; the mean intensity J = Wd * tau0 * (nu/nu0)^p * planck_function(Td,nu) (see e.g. Sobolev et al. 1997, van der Walt 2014)
	double tau_nu0;				// optical depth at frequency=nu0 from the modified black body function
	double nu0;					// [Hz], nu0 from the modified black body function
	double p;					// spectral index
	double Wd;					// dillution factor for the dust emission
	double Td;					// [K], dust temperature
	unsigned short inner_dust_included;	// = 0 - there will be no dust inside the maser region; = 1 - there will be dust inside the maser region
	unsigned short outer_dust_at_LOS;	// = 0 - the dust outside the maser region is not on the line of sight; = 1 - the dust outside the maser region is on the line of sight
	
	// HII region emission, see e.g. Sobolev et al. 1997, J = WHii * (1-exp(tauHII)) * planck_function(Te,nu), tauHII = (HII_turn_freq/nu)^2
	double WHii;				//dillution factor for the HII region emission
	double HII_turn_freq; 		// [Hz], turnover frequency of HII region emission
	double Te; 					// [K], HII region electron temperature
	unsigned short HII_region_included;	// = 0 - there will be no background HII region emission; = 1 - there will be the background HII region emission

	// Cosmic microwave background temperature, Kelvins; J = planck_function(T_CMB,nu)
	double T_CMB;

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
				nu_from_file.push_back( input_nu );
				J_from_file.push_back( input_J ); 
			}
			else break;
		}
		fin.close();
	}

	double interpolate_Jext_from_file(const double & freq) 	// interpolate mean intensity that has been read from file for a given frequency
	{
		if (freq <= nu_from_file[0]) return J_from_file[0];
		else if (freq >= nu_from_file[nu_from_file.size()-1]) return J_from_file[J_from_file.size()-1];
		else {
			for (size_t i = 0; i < nu_from_file.size()-1; i++ ) {
				if (nu_from_file[i] <= freq && freq <= nu_from_file[i+1])
					return J_from_file[i] + (J_from_file[i+1] - J_from_file[i]) / (nu_from_file[i+1] - nu_from_file[i]) * (freq - nu_from_file[i]);  
			} 
		}
		return 0.0;
	}

	void read_Td_from_file(const string & filename) 	//read mean intensity of external emission from a file
	{
		ifstream fin;
		istringstream sfin;
		string str;
		double input_time, input_Td;

		fin.open(filename.c_str(), ios::in);
		if ( !fin.good() ) { 	// check if we found the file
			fin.close();
			throw runtime_error("can't find input file with dependence of dust temperature on time");
		}

		while (getline(fin,str)) {
			str = trim(str); 	// trim is taken from auxiliary.h
			if (str.size() != 0) {
				sfin.str(str);
				sfin >> input_time >> input_Td;
				sfin.clear();
				time_from_file_small.push_back( input_time );
				Td_from_file_small.push_back( input_Td ); 
			}
			else break;
		}
		fin.close();
	}

	double interpolate_Td_from_file(const double & time, const vector <double> & time_from_file, const vector <double> & Td_from_file) 	// interpolate dust temperature that has been read from file for a time
	{
		if (time < time_from_file[0] || time > time_from_file[time_from_file.size()-1]) return 0.0;
		for (size_t i = 0; i < time_from_file.size()-1; i++ ) {
			if (time_from_file[i] <= time && time <= time_from_file[i+1])
				return Td_from_file[i] + (Td_from_file[i+1] - Td_from_file[i]) / (time_from_file[i+1] - time_from_file[i]) * (time - time_from_file[i]);  
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

		for (short i = 0; i < 4; i++) getline(fin, str);
		sfin.str(trim(str));
		sfin >> WHii >> HII_turn_freq >> Te >> HII_region_at_LOS;
		sfin.clear();

		getline(fin, str);
		T_CMB = readline<double>(fin);
	}

	double tau_HII(const double freq) 	// optical depth of the HII region at a given frequency
	{
		return ( HII_turn_freq / freq) * ( HII_turn_freq / freq); 	// see Sobolev & Deguchi 1994, Sobolev et al. 1997
	}

	double var_Td(const double & time) // dependence of dust temperature on time
	{
		//double curTd = 150;
		//if (time > 3600) curTd = 120;
		//
		double curTd = Td;
		double maxTd = 170.;
		if ( time <= 600 ) curTd = curTd + (maxTd - curTd) / 600. * time;
		else curTd = curTd + (maxTd - curTd) * exp( - 0.5*pow((time - 600)/18000., 2.0) );
		return curTd;
	}

	double var_Td(const double & time, const vector <double> & time_from_file, const vector <double> & Td_from_file) // dependence of dust temperature on time
	{
		return Td + interpolate_Td_from_file(time, time_from_file, Td_from_file);	
	}

	double var_Td_small(const double & time) // dependence of dust temperature on time
	{
		double curTd = Td;
		double maxTd = 170.;
		if ( time <= 600 ) curTd = curTd + (maxTd - curTd) / 600. * time;
		else curTd = curTd + (maxTd - curTd) * exp( - 0.5*pow((time - 600)/18000., 2.0) );
		return curTd;
	}

	double var_Td_large(const double & time) // dependence of dust temperature on time
	{
		double curTd = Td;
		double maxTd = 130.;
		if ( time <= 600 ) curTd = curTd + (maxTd - curTd) / 600. * time;
		else curTd = curTd + (maxTd - curTd) * exp( - 0.5*pow((time - 600)/18000., 2.0) );
		return curTd;
	}

public:

	unsigned short HII_region_at_LOS;	// = 0 - HII region is not on the line of sight; = 1 - HII region is on the line of sight
	
	double tau_dust(const double & freq) 	// optical depth of external dust at a given frequency
	{
		return tau_nu0 * pow( freq / nu0, p); // see e.g. Sobolev et al. 1997
		////const double norm = 1 + pow( (4*nu0) / (nu0*0.25), 1.5);
		//const double norm = pow( nu0 / (4*nu0), 2.0) + pow( nu0 / (nu0*0.25), 1.5);
		//return tau_nu0 * ( pow( freq / (4*nu0), 2.0) + pow( freq / (nu0*0.25), 1.5) ) / norm;
	}

	double tau_dust_LOS(const double & freq) 	// optical depth of external dust at a given frequency on the line of sight
	{
		return outer_dust_at_LOS * tau_dust(freq);
	}

	double tau_dust_in(const double & freq, const double & lineWidth) 	// optical depth of the dust inside the maser region at a given frequency and for a given line width
	{
		//return inner_dust_included * 2.6e-25 * lineWidth / modelPhysPars::abundance * modelPhysPars::NdV * pow( freq / nu0, p); // based on Sherwood et al. 1980;
		return inner_dust_included * (0.899 * 1.5e-26) * lineWidth / modelPhysPars::abundance * modelPhysPars::NdV * pow( freq / nu0, p); // based on Ossenkopf & Henning et al. 1994;
	}

	double outer_dust_source_function(const double & freq) 	// source function for the external dust emission
	{
		return planck_function(Td, freq);
		/*double time = 0;
		double k_ratio = pow( freq / (nu0*0.25), 1.5) / pow( freq / (4*nu0), 2.0);
		double j = planck_function(var_Td_small(time), freq) + planck_function(var_Td_large(time), freq) * k_ratio;
		double k = 1. + k_ratio;
		return j / k;*/
	}

	double outer_dust_source_function(const double & freq, const double & time) 	// source function for the external dust emission
	{
		return planck_function(var_Td(time), freq);
		/*double k_ratio = pow( freq / (nu0*0.25), 1.5) / pow( freq / (4*nu0), 2.0);
		double j = planck_function(var_Td_small(time), freq) + planck_function(var_Td_large(time), freq) * k_ratio;
		double k = 1. + k_ratio;
		return j / k;*/
	}

	double inner_dust_source_function(const double & freq) 	// source function for the dust emission inside the maser region
	{
		return inner_dust_included * planck_function(Td, freq);
	}

	double inner_dust_source_function(const double & freq, const double & time) 	// source function for the dust emission inside the maser region
	{
		return inner_dust_included * planck_function(var_Td(time), freq);
	}
	
	double compute_Jext_dust_CMB_file(const double & freq)		// returns mean intensity taken from file or sum of dust and cosmic microwave background (CMB) mean intensities
	{
		if (read_Jext_from_file == 0) {
			const double J_dust = Wd * ( 1. - exp(-tau_dust(freq)) ) * outer_dust_source_function(freq); // Note that this a general formula for dust emission (e.g. Sobolev et al. 1997), while eq. 4 in van der Walt 2014 gives optically thin case

			const double J_CMB = Wd * exp(-tau_dust(freq)) * planck_function(T_CMB,freq) + (1. - Wd) * planck_function(T_CMB,freq); // cosmic microwave background emission taking into account an absorption by the external dust layer

			return J_dust + J_CMB;
		} else {
			// !!!!!! NOTE !!!!!! that this code has been used for calculations with the mean intensity provided by CLOUDY code;
			// CLOUDY outputs photon occupation number which is converted into the mean intensity below
			return interpolate_Jext_from_file(freq) * 2.*PLANK_CONSTANT*pow(freq/SPEED_OF_LIGHT,2.)*freq; // photon occupation number into mean intensity;
		}
	}

	double compute_Jext_dust_CMB_file(const double & freq, const double & time)		// returns mean intensity taken from file or sum of dust and cosmic microwave background (CMB) mean intensities
	{
		if (read_Jext_from_file == 0) {
			const double J_dust = Wd * ( 1. - exp(-tau_dust(freq)) ) * outer_dust_source_function(freq, time); // Note that this a general formula for dust emission (e.g. Sobolev et al. 1997), while eq. 4 in van der Walt 2014 gives optically thin case

			const double J_CMB = Wd * exp(-tau_dust(freq)) * planck_function(T_CMB,freq) + (1. - Wd) * planck_function(T_CMB,freq); // cosmic microwave background emission taking into account an absorption by the external dust layer

			return J_dust + J_CMB;
		} else {
			// !!!!!! NOTE !!!!!! that this code has been used for calculations with the mean intensity provided by CLOUDY code;
			// CLOUDY outputs photon occupation number which is converted into the mean intensity below
			return interpolate_Jext_from_file(freq) * 2.*PLANK_CONSTANT*pow(freq/SPEED_OF_LIGHT,2.)*freq; // photon occupation number into mean intensity;
		}
	}

	double compute_JextHII(const double & freq)		// mean intensity of the HII region emission, needs to be separated from other types of external emission because of beaming
	{
		if (read_Jext_from_file == 0)
			return WHii * ( 1. - exp(-tau_HII(freq)) ) * planck_function(Te,freq);
		else
			return 0.0;
	}

	double continuum(const double & freq) 	// computes continuum part in the equation for Tb which is similar to the one given in Appendix A of Sobolev et al. 1997
	{
		return outer_dust_source_function(freq) * (1. - exp(-tau_dust_LOS(freq))) + 
			   planck_function(T_CMB,freq) * exp(-tau_dust_LOS(freq)) + HII_region_included*( 1. - exp(-tau_HII(freq)) ) * planck_function(Te,freq);
	}

	double continuum(const double & freq, const double & time) 	// computes continuum part in the equation for Tb which is similar to the one given in Appendix A of Sobolev et al. 1997
	{
		return outer_dust_source_function(freq, time) * (1. - exp(-tau_dust_LOS(freq))) + 
			   planck_function(T_CMB,freq) * exp(-tau_dust_LOS(freq)) + HII_region_included*( 1. - exp(-tau_HII(freq)) ) * planck_function(Te,freq);
	}

	template <typename T>
	dust_HII_CMB_Jext_radiation(T & fin)
	{
		T_CMB = 2.728;
		Te = 8000.0;
		WHii = 0.0;
		HII_turn_freq = 1.e8;
		Td = T_CMB;
		Wd = 0.0;
		p = 2.0;
		nu0 = 300.e9;
		tau_nu0 = 1.0;
		read_Jext_from_file = 0;
		inner_dust_included = 0;
		HII_region_included = 1;
		HII_region_at_LOS = 1;
		outer_dust_at_LOS = 1;
		read_parameters(fin);
		if (Wd   <= 0.5) outer_dust_at_LOS = 0;
		if (WHii <= 100. * DBL_EPSILON || HII_region_at_LOS == 0) HII_region_included = 0;
		//read_Td_from_file("Td.txt");
	}

	~dust_HII_CMB_Jext_radiation()
	{
		nu_from_file.clear();
		J_from_file.clear();
		time_from_file_small.clear();
		Td_from_file_small.clear();
		time_from_file_large.clear();
		Td_from_file_large.clear();
	}
};
