#pragma once
#include "spline.h"
#include "hiddenParameters.h"
#include <functional>
#include "cubature/cubature.h"

typedef struct {
	double sig;
	double tau;
} fdata_struct;

int integrand_for_cubature(unsigned ndim, const double *x, void * fdata_in, unsigned fdim, double * fval) // function under the integral in equation A1 from Langer & Watson
{
	fdata_struct fdata = *(fdata_struct *)fdata_in;
	const double sig = fdata.sig;
	const double tau = fdata.tau;
	//double integ_fac = ( 1. + sig * y *y );
	const double integ_fac = fma(sig * x[0], x[0], 1.0);

	fval[0] = integ_fac * exp( - tau/integ_fac ); 	// = (1+sigma*y*y) * exp[ -tau/(1+sigma*y*y) ]
	return 0;
}

class beta_LVG // class to compute escape probability, beta; see equations from Castor 1970 and Langer & Watson 1984
{
private:

	double tauCutOff;	    // tau limit below which beta=(1-exp(-tau))/tau will be expanded as a Taylor series
	double sig;        		// equation A6 from Langer & Watson
	double oPlusSigDiv3; 	// = (1 + sig / 3)
	double tau_min; 		// minimum optical depth for which there will be precomputed values of beta in the case of sig != 0
	double tau_max; 		// maximum optical depth for which there will be precomputed values of beta in the case of sig != 0

	akima_spline spline_neg; // structure for akima spline interpolation of beta for negative taus
	akima_spline spline_pos; // structure for akima spline interpolation of beta for positive taus

	double integrate_F(const double & tau) 		// the integral in equation A1 from Langer & Watson
	{
		fdata_struct fdata = {sig, tau};
		const double x_min[1] = {0.};
        const double x_max[1] = {1.};
        double integ_res[1];
        double integ_err[1];
        int integration_success = 0;
        integration_success = hcubature(1, integrand_for_cubature, &fdata, 1, x_min, x_max, 100000000, 1.e-60, BETA_ACCURACY, ERROR_INDIVIDUAL, integ_res, integ_err);
		if (integration_success != 0) throw runtime_error("integration for escape probability has failed");
		return integ_res[0];
	}

	void precalculate_beta() 	// precalculates the values of beta for different values of the optical depth and for a given value of sig;
	{ // beta is computed with equation A2 from Langer & Watson
		const double tau1 = tau_min;
		const double tau2 = 0.0;
		const double dtau1 = (tau2-tau1) / double(BETA_NUM_SPLINE_POINTS);
		const double dtau2 = (tau_max-tau2) / double(BETA_NUM_SPLINE_POINTS);
		double tau, F;
		vector <double> tau_interp_array; 		// array of optical depths for which there will be precomputed values of beta; these values will be used to create the spline for interpolation
		vector <double> log_beta_interp_array; 	// array of precomputed values of log(beta); these values will be used to create the spline for interpolation

		tau = tau1;
		for (size_t i = 0; i < BETA_NUM_SPLINE_POINTS; i++) { 	// computing beta for different values of optical depth which are < 0
			tau_interp_array.push_back(tau);
			F = integrate_F(tau);
			log_beta_interp_array.push_back( log( (oPlusSigDiv3 - F) / tau ) ); 	// see equation A2 from Langer & Watson 1984, integral F is computed numerically
			tau += dtau1;
		}
		tau_interp_array.push_back(0.0e0);
		log_beta_interp_array.push_back(0.0e0);
		spline_neg.akima_init(tau_interp_array, log_beta_interp_array); 	// creating spline for negative optical depths

		tau_interp_array.clear();
		log_beta_interp_array.clear();

		tau_interp_array.push_back(0.0e0);
		log_beta_interp_array.push_back(0.0e0);
		tau = tau2 + dtau2;
		for (size_t i = 0; i < BETA_NUM_SPLINE_POINTS; i++) { 	// computing beta for different values of optical depth which are > 0
			tau_interp_array.push_back(tau);
			F = integrate_F(tau);
			log_beta_interp_array.push_back( log( (oPlusSigDiv3 - F) / tau ) ); 	// see equation A2 from Langer & Watson 1984, integral F is computed numerically
			tau += dtau2;
		}
		spline_pos.akima_init(tau_interp_array, log_beta_interp_array); 	// creating spline for positive optical depths

		tau_interp_array.clear();
		log_beta_interp_array.clear();
	}

	double beta_noBeam(const double & tau)			// LVG escape probability without beaming, see Appendix A in Langer & Watson 1984, Castor 1970
	{
		// without beaming, beta = (1 - exp(-tau)) / tau
		if (fabs(tau) < tauCutOff) return 1. - tau * (1. - tau*(1. / 3.)) * 0.5; //Taylor expansion in the case of small absolute value of tau
		return (1. - exp(-tau)) / tau;
	}

	double beta_Beam(const double & tau)			// LVG escape probability with beaming, see Appendix A in Langer & Watson 1984, Castor 1970
	{
		if (tau > tau_max) return 1. / tau * oPlusSigDiv3; 	// see equation A.6 from Castor 1970	
		if (tau < 0.0) return exp(spline_neg.akima_eval(tau));
		else return exp(spline_pos.akima_eval(tau));
	}

	double tauDbetaDtau_noBeam(const double & tau)			// = (tau * derivative of beta on tau); without beaming
	{
		if (fabs(tau) < tauCutOff) return (tau * (1. / 3.) - 0.5) * tau;
		const double expTau = exp(-tau);
		return expTau - (1. - expTau)/tau;
	}

	double tauDbetaDtau_Beam(const double & tau)			// = (tau * derivative of beta on tau); with beaming
	{
		if (tau > tau_max) return  - 1. / tau * oPlusSigDiv3; // see equation A.6 from Castor 1970	
		if (tau < 0.0) return tau * exp(spline_neg.akima_eval(tau)) * spline_neg.akima_eval_deriv(tau);
		else return tau * exp(spline_pos.akima_eval(tau)) * spline_pos.akima_eval_deriv(tau);
	}

public:

	std::function<double(const double &)> beta; 			// function that returns LVG escape probability with or without beaming
	std::function<double(const double &)> tauDbetaDtau; 	// function that returns (tau * derivative of beta on tau) with or without beaming
	//double (*beta)(const double &);
	//double (*tauDbetaDtau)(const double &);

	double betaHII_LOS(const double & taui, const double & beamH)			// LVG escape probability for the HII region backgroung radiation in the case if the HII region is on the line of sight, see Appendix A in Sobolev et al. 1997
	{
		const double tau = taui * beamH;
		return beta_noBeam(tau);
	}

	double tauDbetaHIIDtau_LOS(const double & taui, const double & beamH)			// (tau *dbeta/dtau) for the HII region backgroung radiation in the case if the HII region is on the line of sight, see Appendix A in Sobolev et al. 1997
	{
		const double tau = taui * beamH;
		return tauDbetaDtau_noBeam(tau);
	}

	double betaHII_pump(const double & tau, const double & beamH)			// LVG escape probability for the HII region backgroung radiation in the case if the HII region is not on the line of sight; in a similar manner as in Appendix A from Sobolev et al. 1997
	{
		return beta_noBeam(tau);
	}

	double tauDbetaHIIDtau_pump(const double & tau, const double & beamH)			// (tau *dbeta/dtau) of the LVG escape probability for the HII region backgroung radiation  in the case if the HII region is not on the line of sight; in a similar manner as in Appendix A from Sobolev et al. 1997
	{
		return tauDbetaDtau_noBeam(tau);
	}

	beta_LVG(const double & beamH)
	{
		sig = 1.0e00 / beamH - 1.0e00;					// equation A6 from Langer & Watson 1984
		oPlusSigDiv3 = 1.e0 + sig / 3.e0;
		tauCutOff = pow((24.*DBL_EPSILON), 0.25); 	// taken from LIME code (Brinch & Hogerheijde 2010)
		tau_min = MIN_TAU; 							// MIN_TAU is defined in hiddenParameters.h
		tau_max = 100. * max(1., 1. + sig); 		// this optical depth limit is consistent with equation A.6 from Castor 1970
		if (fabs(sig) >= (100.*DBL_EPSILON)) { 		// the case with beaming;
			precalculate_beta();
			beta = [this](const double & tau) -> double {
				return this->beta_Beam(tau);
			};
			tauDbetaDtau = [this](const double & tau) -> double {
				return this->tauDbetaDtau_Beam(tau);
			};
			//beta = &beta_LVG::beta_Beam;
			//tauDbetaDtau = &beta_LVG::tauDbetaDtau_Beam;
		} else { 									//no beaming
			beta = [this](const double & tau) -> double {
				return this->beta_noBeam(tau);
			};
			tauDbetaDtau = [this](const double & tau) -> double {
				return this->tauDbetaDtau_noBeam(tau);
			};
			//beta = &beta_LVG::beta_noBeam;
			//tauDbetaDtau = &beta_LVG::tauDbetaDtau_noBeam;
		}
	}

	~beta_LVG()
	{
		spline_neg.clear();
		spline_pos.clear();
	}
};
