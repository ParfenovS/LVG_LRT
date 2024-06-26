#pragma once
#include "spline.h"
#include "hiddenParameters.h"
#include <functional>
#include "cubature/cubature.h"
//#include"gamma.h"

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
	//double epsb;			// = 1. / beamH;
	double tau_min; 		// minimum optical depth for which there will be precomputed values of beta in the case of sig != 0
	double tau_max; 		// maximum optical depth for which there will be precomputed values of beta in the case of sig != 0

	akima_spline spline_neg; // structure for akima spline interpolation of beta for negative taus
	akima_spline spline_pos; // structure for akima spline interpolation of beta for positive taus
	
	/*double F1(const double & lam, const double & phi, const double & x) // equation A12 from Langer & Watson
	{
		const double phi2 = phi * phi;
		const double phi3 = phi2 * phi;
		const double d = (1. - phi) / (lam * phi); // equation A13
		const double d2 = d * d;
		const double d3 = d2 * d;
		const double a = (1 - x * x) * sig * lam / ((1. + sig) * (1 + sig * x * x)); // equation A9

		return pow(1. - phi, 3.) / (2. * phi * lam) * exp(lam / (1. - phi)) * (
			lowgamma(1, -a) +
			0.5 * d * (1 + 5 * phi) * lowgamma(2, -a) + 
			(1./ 8.) * d2 * (3 + 10 * phi + 35 * phi2) * lowgamma(3, -a) +
			(5./ 16.) * d3 * (1 + 3 * phi + 7 * phi2 + 21 * phi3) * lowgamma(4, -a) + 
			(5./ 128.) * d2 * d2 * (7 + 20 * phi + 42 * phi2 + 84 * phi3 + 231 * phi2 * phi2) * lowgamma(5, -a) + 
			(7./ 256.) * d3 * d2 * (9 + 25 * phi + 50 * phi2 + 90 * phi3 + 165 * phi2 * phi2 + 429 * phi3 * phi2) * lowgamma(6, -a)
		);
	}

	double F3(const double & lam, const double & x) // equation A16 from Langer & Watson
	{
		const double sig2 = sig * sig;
		const double sig3 = sig2 * sig;
		const double x3 = x * x * x;
		const double lam2 = lam * lam;

		return exp(lam) * (
			1 - x - (1 - x3) / 3. * sig * (lam - 1) + (1 - x3*x*x) * sig2 * 0.1 * lam2 -
			(1 - x3 * x3 * x) * sig3 * (1. / 14.) * (lam / 3. + 1) * lam2 +
			(1 - x3 * x3 * x3) * sig2 * sig2 * (1. / 9.) * (lam2 / 24. + lam / 3. + 0.5) * lam2 -
			(1 - x3 * x3 * x3 * x * x) * sig3 * sig2 * (1. / 22.) * (lam2 * lam / 60. + lam2 * 0.25 + lam + 1) * lam2
		);
	}

	double G5(const double & lam, const double & phi, const double & x) // equation A21 from Langer & Watson
	{
		const double oneminx = 1 / (1 - phi * x * x);
		const double oneminx2 = oneminx * oneminx;
		const double sphi = sqrt(phi);
		const double Lt = log( (1 + sphi) / (1 - sphi) * (1 - sphi * x) / (1 + sphi * x)) / sphi; // equation A22 from Langer & Watson
		const double lam2 = lam * lam;
		const double inveps = 1. / epsb;
		const double inveps2 = inveps * inveps;

		return (1 - x + lam * 0.25 * Lt + // note, there is no factor of 1/2 in comparison to eq. A21 from Langer & Watson
			lam2 * (inveps - x * oneminx + Lt * 0.5) * (1. / 12.) +
			lam2 * lam * (inveps2 - x * oneminx2 + 1.5 * inveps - 1.5 * x * oneminx + Lt * (3./4.)) * (1. / 96.) +
			lam2 * lam2 * (inveps2 * inveps * (1./3.) - x * oneminx2 * oneminx * (1./3.) + inveps2 * (5./12.) - x * oneminx2 * (5./12.) + inveps * (5./8.) - x * oneminx  * (5./8.) + Lt * (5./16.)) * (1. / 240.) + 
			lam2 * lam2 * lam * (inveps2 * inveps2 * (1./8.) - x * oneminx2 * oneminx2 * (1./8.) + inveps2 * inveps * (7./48.) - x * oneminx2 * oneminx * (7./48.) + inveps2 * (35./192.) - x * oneminx2 * (35./192.) + inveps * (35./128.) - x * oneminx * (35./128.) + Lt * (35./256.)) * (1. / 720.) // note there is a difference in power of one of the members in comparison to eq. A21 from Langer & Watson
		);
	}*/

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
		/*if (tau < 0.0) { //return exp(spline_neg.akima_eval(tau));
			if (sig < 0.0) {
				const double lam = -tau; // equation A.8 from Langer & Watson
				const double ilam = 1.0 / lam; // equation A.8 from Langer & Watson
				const double phi = -sig; // equation A.11 from Langer & Watson
				double xtilde = 0.0;
				if (lam > 3.5) xtilde = 0.0;
				else if (lam < 3.5 * epsb) xtilde = 1.0;
				else xtilde = sqrt((1 - lam / 3.5) / phi); // equation A.23 from Langer & Watson

				if (-1 < sig && sig <= -0.4) {
					return G5(lam, phi, 0) - G5(lam, phi, xtilde) - ilam * (1 - xtilde - phi * (1 - xtilde*xtilde*xtilde) * (1./3.)) + ilam * F1(lam, phi, xtilde);
				} else {
					if (lam < 1 + 42 *(1 - 4 * phi + 4.17 * phi * phi)) return G5(lam, phi, xtilde) - ilam * (1 - xtilde - phi * (1 - xtilde*xtilde*xtilde) * (1./3.)) + ilam * F3(lam, xtilde);
					return G5(lam, phi, xtilde) - ilam * (1 - xtilde - phi * (1 - xtilde*xtilde*xtilde) * (1./3.)) + ilam * F1(lam, phi, xtilde);
				}
			} else return exp(spline_neg.akima_eval(tau));
		}*/
		if (tau < 0.0) return exp(spline_neg.akima_eval(tau));
		else return exp(spline_pos.akima_eval(tau));
	}

	double DbetaDtau_noBeam(const double & tau)			// = derivative of beta on tau; without beaming
	{
		if (fabs(tau) < tauCutOff) return tau * (1. / 3.) - 0.5;
		const double expTau = exp(-tau);
		const double inv_tau = 1 / tau;
		return (expTau - (1. - expTau) * inv_tau ) * inv_tau;
	}

	double DbetaDtau_Beam(const double & tau)			// = derivative of beta on tau; with beaming
	{
		if (tau > tau_max) return  - 1. / (tau * tau) * oPlusSigDiv3; // see equation A.6 from Castor 1970	
		if (tau < 0.0) return exp(spline_neg.akima_eval(tau)) * spline_neg.akima_eval_deriv(tau);
		else return exp(spline_pos.akima_eval(tau)) * spline_pos.akima_eval_deriv(tau);
	}

public:

	std::function<double(const double &)> beta; 			// function that returns LVG escape probability with or without beaming
	std::function<double(const double &)> DbetaDtau; 	// function that returns derivative of beta on tau with or without beaming
	//double (*beta)(const double &);
	//double (*DbetaDtau)(const double &);

	double betaHII_LOS(const double & taui, const double & beamH)			// LVG escape probability for the HII region backgroung radiation in the case if the HII region is on the line of sight, see Appendix A in Sobolev et al. 1997
	{
		const double tau = taui * beamH;
		return beta_noBeam(tau);
	}

	double DbetaHIIDtau_LOS(const double & taui, const double & beamH)			// dbeta/dtau for the HII region backgroung radiation in the case if the HII region is on the line of sight, see Appendix A in Sobolev et al. 1997
	{
		const double tau = taui * beamH;
		return DbetaDtau_noBeam(tau);
	}

	double betaHII_pump(const double & tau, const double & beamH)			// LVG escape probability for the HII region backgroung radiation in the case if the HII region is not on the line of sight; in a similar manner as in Appendix A from Sobolev et al. 1997
	{
		return beta_noBeam(tau);
	}

	double DbetaHIIDtau_pump(const double & tau, const double & beamH)			// dbeta/dtau of the LVG escape probability for the HII region backgroung radiation  in the case if the HII region is not on the line of sight; in a similar manner as in Appendix A from Sobolev et al. 1997
	{
		return DbetaDtau_noBeam(tau);
	}

	beta_LVG(const double & beamH)
	{
		//epsb = 1.0e00 / beamH;
		sig = 1.0e00 / beamH - 1.0e00;					// equation A6 from Langer & Watson 1984
		oPlusSigDiv3 = 1.e0 + sig / 3.e0;
		tauCutOff = pow((24.*DBL_EPSILON), 0.25); 	// taken from LIME code (Brinch & Hogerheijde 2010)
		tau_min = max(MIN_TAU, -705 * (sig + 1));	// MIN_TAU is defined in hiddenParameters.h, the limit is set to avoid an overflow in integrand_for_cubature
		tau_max = 100. * max(1., 1. + sig); 		// this optical depth limit is consistent with equation A.6 from Castor 1970
		if constexpr (ESCAPE_PROBABILITY_METHOD == 0) {
			if (fabs(sig) >= (100.*DBL_EPSILON)) { 		// the case with beaming;
				precalculate_beta();
				beta = [this](const double & tau) -> double {
					return this->beta_Beam(tau);
				};
				DbetaDtau = [this](const double & tau) -> double {
					return this->DbetaDtau_Beam(tau);
				};
				//beta = &beta_LVG::beta_Beam;
				//DbetaDtau = &beta_LVG::DbetaDtau_Beam;
			} else { 									//no beaming
				beta = [this](const double & tau) -> double {
					return this->beta_noBeam(tau);
				};
				DbetaDtau = [this](const double & tau) -> double {
					return this->DbetaDtau_noBeam(tau);
				};
				//beta = &beta_LVG::beta_noBeam;
				//DbetaDtau = &beta_LVG::DbetaDtau_noBeam;
			}
		} else {
			if constexpr (ESCAPE_PROBABILITY_METHOD == 1) {
				//     Uniform sphere formula from Osterbrock (Astrophysics of
				//     Gaseous Nebulae and Active Galactic Nuclei) Appendix 2
				//     with power law approximations for large and small tau
				beta = [](const double & tau) -> double {
					const double taur = tau * 0.5;
					if( fabs(taur) < 0.02) return 1.0 - 0.75 * taur + (taur * taur) / 2.5 
												- pow(taur, 3.) / 6.0 + pow(taur, 4.) / 17.5;
					else if (fabs(taur) > 5.e2) return 0.75 / taur;
					else return 0.75 / taur * (1. - 1. / (2 * (taur * taur)) +
								(1. / taur + 1. / (2*(taur * taur))) * exp(-2.*taur));
				};
				DbetaDtau = [](const double & tau) -> double {
					const double taur = tau * 0.5; // note that Dbeta/Dtau = Dtaur/Dtau * Dbeta/Dtaur = 0.5 * Dbeta/Dtaur
					if( fabs(taur) < 0.02) return 0.5 * ( 0.2285714285714286 * pow(taur, 3.) - 0.5 * (taur * taur) + 0.8 * taur - 0.75 );
					else if (fabs(taur) > 5.e2) return 0.5 * ( - 0.75 / (taur * taur) );
					else return 0.5 * ( -(exp(-2 * taur) * ((6 * (taur*taur) - 9) * exp(2*taur) + 12 * (taur*taur) + 18 * taur + 9)) / (8. * pow(taur, 4.)) );
				};
			} else if constexpr (ESCAPE_PROBABILITY_METHOD == 2) {
				//     Expanding sphere = Large Velocity Gradient (LVG) or Sobolev case.
				//     Formula from De Jong, Boland and Dalgarno (1980, A&A 91, 68)
				//     corrected by factor 2 in order to match beta(tau=0)=1
				beta = [](const double & tau) -> double {
					const double taur = tau * 0.5;
					if (fabs(taur) < 1.e-8) return 1.0;
					else if (fabs(taur) < 6.96756) return 2.0 * (1.0 - exp(-2.34 * taur)) / (4.68 * taur);
					else return 2.0 / (taur * 4.0 * (sqrt(log(taur / sqrt(PI)))));
				};
				DbetaDtau = [](const double & tau) -> double {
					const double taur = tau * 0.5;
					if (fabs(taur) < 1.e-8) return 0.0;
					else if (fabs(taur) < 6.96756) return 0.5 * ( (exp(-2.34 * taur)) / taur - (0.4273504273504274 * (1. - exp(-2.34 * taur))) / (taur * taur) );
					else return 0.5 * ( -0.5 / ((taur*taur) * sqrt(log(taur / sqrt(PI)))) - 0.25 / ((taur*taur) * pow(log(taur / sqrt(PI)), (1.5))) );
				};
			} else if constexpr (ESCAPE_PROBABILITY_METHOD == 3) {
				//     Slab geometry (e.g., shocks): de Jong, Dalgarno & Chu 1975,
				//     ApJ 199, 69 (again with power law approximations)
				beta = [](const double & tau) -> double {
					if (fabs(3.0 * tau) < 1.e-6) return 1.0 - 1.5 * (tau - tau*tau);
					else if (fabs(3.0 * tau) > 50.0) return 1. / (3.0 * tau);
					else return (1. - exp(-3.0 * tau)) / (3.0 * tau);
				};
				DbetaDtau = [](const double & tau) -> double {
					if (fabs(3.0 * tau) < 1.e-6) return -1.5 * (1 - 2 * tau);
					else if (fabs(3.0 * tau) > 50.0) return - 1. / (3.0 * tau * tau);
					else return (exp(-3.0 * tau)) / tau - (1 - exp(-3.0 * tau)) / (3.0 * tau*tau);
				};
			}
		}
	}

	~beta_LVG()
	{
		spline_neg.clear();
		spline_pos.clear();
	}
};
