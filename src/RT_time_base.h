#pragma once
#include "RT.h"
#include "physconsts.h"
#include <fstream>

class RT_time_base : public RT		// base class for integration of kinetic equations for level populations over time;
{
protected:

	void prepare_results_for_output(beta_LVG & LVG_beta, const double & time)
	{
		// preparing quantities for output
		double dummy_S = 0.0;
		double dummy_beta = 0.0;
		double dummy_betaS = 0.0;
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			compute_tau(i); 		// computing final optical depths that can be used for output
			compute_J_S_beta(i, LVG_beta, dummy_S, dummy_beta, dummy_betaS); 	// computing final mean intensities that can be used for output
			compute_Tex(i); 		// computing excitation temperature that can be used for output
			compute_brightness_temperature(i, time); 	// computing brightness temperature and intensity of the emission
		}
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			// the output optical depth is an optical depth along the line of sight, thus one need to take into account beaming
			mol->rad_trans[i].tau *= beamH;
		}
	}

	void write_rad_trans_data_into_binary_file(beta_LVG & LVG_beta, const double & time, ofstream & binTbrfile, ofstream & binTexfile, ofstream & bintaufile) {
		binTbrfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
		binTexfile.write(reinterpret_cast<const char*>(&time), sizeof(double));
		bintaufile.write(reinterpret_cast<const char*>(&time), sizeof(double));
		prepare_results_for_output(LVG_beta, time);
		double temp_var = 0.0;
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			temp_var = mol->rad_trans[i].Tbr;
			binTbrfile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
			temp_var = mol->rad_trans[i].Tex;
			binTexfile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
			temp_var = mol->rad_trans[i].tau;
			bintaufile.write(reinterpret_cast<const char*>(&temp_var), sizeof(double));
		}
	}

	void clear_mem_close_files(double A[], double pop[], double dpop_dt[], vector <vector <double> > & oldpops_Ng, vector <vector <double> > & oldpops_time, ofstream & binpopfile, ofstream & binTbrfile, ofstream & binTexfile, ofstream & bintaufile)
	{
		delete[] A; delete[] pop; delete[] dpop_dt;
		oldpops_Ng.clear(); oldpops_time.clear();
		if (cerr_output_iter_progress) {
			binpopfile.close(); binTbrfile.close(); binTexfile.close(); bintaufile.close();
		}
	}

	void clear_mem_close_files(double A[], double pop[], double Jac[], double dpop_dt[], vector <vector <double> > & oldpops_Ng, vector <vector <double> > & oldpops_time, ofstream & binpopfile, ofstream & binTbrfile, ofstream & binTexfile, ofstream & bintaufile)
	{
		delete[] A; delete[] pop; delete[] Jac; delete[] dpop_dt;
		oldpops_Ng.clear(); oldpops_time.clear();
		if (cerr_output_iter_progress) {
			binpopfile.close(); binTbrfile.close(); binTexfile.close(); bintaufile.close();
		}
	}

	void update_external_emission(const double & time)
	{
		// set the external emission mean intensity, optical depth, absorption and emission coefficients
		for (size_t i = 0; i < mol->rad_trans.size(); i++) {
			mol->rad_trans[i].JExt = dust_HII_CMB_Jext_emission->compute_Jext_dust_CMB_file(mol->rad_trans[i].nu, time); //external emission from dust, cosmic microwave background or file
			mol->rad_trans[i].emiss_dust = mol->rad_trans[i].kabs_dust * dust_HII_CMB_Jext_emission->inner_dust_source_function(mol->rad_trans[i].nu, time); //absorption coefficient of the dust inside the maser region
		}
	}

	void compute_brightness_temperature(const size_t & i, const double & time)		// computes brightness temperature and intensity for radiative transition i
	{ 	// it is similar to the last Equation for Tbr in Appendix A of Sobolev et al. 1997
		const double & nu = mol->rad_trans[i].nu;
		mol->rad_trans[i].Tbr = exp(-dust_HII_CMB_Jext_emission->tau_dust_LOS(nu)) * (
			(1. - exp(-mol->rad_trans[i].tau*beamH)) * compute_source_function(i) +
			(exp(-mol->rad_trans[i].tau*beamH) - exp(-mol->rad_trans[i].taud_in*beamH)) * dust_HII_CMB_Jext_emission->continuum(nu, time) -
			(1. - exp(-mol->rad_trans[i].taud_in*beamH)) * dust_HII_CMB_Jext_emission->inner_dust_source_function(nu, time)
		) * (pow(SPEED_OF_LIGHT/nu, 2.0) / (2. * BOLTZMANN_CONSTANT));
	}

public:

	RT_time_base() noexcept : RT()
	{}

	RT_time_base(istream & cin) : RT(cin)
	{}

	RT_time_base(const unsigned short & initialSolutionSource, const double & MAX_POPS_EPS, const unsigned long & maxNumberOfIterations, const double & beamH, const double & lineWidth) : RT(initialSolutionSource, MAX_POPS_EPS, maxNumberOfIterations, beamH, lineWidth)
	{}
};
