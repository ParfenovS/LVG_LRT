#pragma once
#include "spline.h"
#include "hiddenParameters.h"
#include <functional>

class beta_LVG // class to compute escape probability, beta; see equations from Castor 1970 and Langer & Watson 1984
{
private:

	double tauCutOff;	    // tau limit below which beta=(1-exp(-tau))/tau will be expanded as a Taylor series
	double sig;        		// equation A6 from Langer & Watson
	double oPlusSigDiv3; 	// = (1 + sig / 3)
	double tau_min; 		// minimum optical depth for which there will be precomputed values of beta in the case of sig != 0
	double tau_max; 		// maximum optical depth for which there will be precomputed values of beta in the case of sig != 0
	
	// the following objects will be used for interpolation of precomputed beta values on tau with splines
	tk::spline beta_spline_neg; 	// for negative optical depths
	tk::spline beta_spline_pos; 	// for positive optical depths
	
	double integrate_F(const double & tau) 		// the integral in equation A1 from Langer & Watson
	{
		unsigned int n;
		double F, F_old, y, step, sum;

		auto integrand = [&] (const double & y) -> double {
			double integ_fac = ( 1. + sig * y *y );
			return integ_fac * exp( - tau/integ_fac ); 	// = (1+sigma*y*y) * exp[ -tau/(1+sigma*y*y) ]
		}; // function under the integral in equation A1 from Langer & Watson

		// the following is an integration with trapezoidal rule with adaptive step
		sum = ( integrand(0.0) + integrand(1.0) ) * 0.5;
		F_old = sum;
		n = 2;
		do {
			step = 1.0 / n;
			for (unsigned int i = 1; i <= n/2; i++) {
				y = (2. * i - 1.) * step;
				sum += integrand(y);
			}
			F = sum * step;

			if ( fabs(F - F_old) < BETA_ACCURACY * fabs(F_old) ) return F; 	// BETA_ACCURACY is defined in hiddenParameters.h

			n = 2 * n;
			F_old = F; 
		} while (true);

		return 0.0;
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
		beta_spline_neg.set_points(tau_interp_array, log_beta_interp_array); 	// creating spline for negative optical depths

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
		beta_spline_pos.set_points(tau_interp_array, log_beta_interp_array); 	// creating spline for positive optical depths

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
		if (tau < 0.0) return exp(beta_spline_neg(tau));
		else return exp(beta_spline_pos(tau));
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
		if (tau < 0.0) return tau * exp(beta_spline_neg(tau)) * beta_spline_neg.deriv(1, tau);
		else return tau * exp(beta_spline_pos(tau)) * beta_spline_pos.deriv(1, tau);
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
	}
};

/*
double beta_LVG::tauCutOff;
double beta_LVG::sig;
double beta_LVG::oPlusSigDiv3;
double beta_LVG::tau_min;
double beta_LVG::tau_max;

tk::spline beta_LVG::beta_spline_neg; // interpolation for negative optical depths
tk::spline beta_LVG::beta_spline_pos; // interpolation for positive optical depths
*/
