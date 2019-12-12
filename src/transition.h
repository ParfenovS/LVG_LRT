#pragma once
#include <vector>

class radiative_transition		// stores the data for a given radiative transition
{
private:

	struct lineOverlapping 		// this structure is used to store the data on lines which are overlpapped with a given transition; the overlapping is local only
	{
		size_t id; 				// id of radiative transition which is blended with a given transition
		double fac; 			// the fraction of the blended line which overlaps with a given transition

		lineOverlapping(const size_t & blend_id, const double & vel_fac) {
			id = blend_id;
			fac = vel_fac;
		}

		~lineOverlapping() {}
	};

public:

	size_t trans_id;
	size_t up_level, low_level;
	double A, Bul, Blu;
	double nu;
	double tau;
	double J, JExt, JExtHII;
	double taud_in, emiss_dust, kabs_dust;
	double Tex, Tbr;
	std::vector <lineOverlapping> blends; //the definition of lineOverlapping structure is in lineOverlapping.h file

	radiative_transition() noexcept
	{
		trans_id = 0; 			// begins from 1;
		A = 0.0; 				// spontaneous Einstein coeff. A [1/s]; can be set in molModel.h with set_trans function;
		Bul = 0.0; 				// Einstein coeff. B from upper to lower level; can be set in molModel.h with set_trans function;
		Blu = 0.0; 				// Einstein coeff. B from lower to upper level; can be set in molModel.h with set_trans function;
		nu = 0.0;				// transition frequency [Hz]; can be set in molModel.h with set_trans function;
		tau = 0.0;				// optical depth at the transition
		Tex = 0.0; 				// excitation temperature, [K]
		J = 0.0;				// total mean intensity [erg*s/cm3]; can be set in molModel.h with set_trans function;
		JExt = 0.0;				// mean intensity [erg*s/cm3] of the external radiation from dust, cosmic microwave background or file
		JExtHII = 0.0;			// mean intensity of the external radiation from HII region, needs to be separated from JExt because of the maser beaming (see equation for Jav in Appendix A of Sobolev et al. 1997)
		taud_in = 0.0;			// optical depth of the dust inside the maser region
		emiss_dust = 0.0;		// emission coefficient of the dust inside the maser region
		kabs_dust = 0.0;		// absorption coefficient of the dust inside the maser region
		Tbr = 0.0;				// brightness temperature for a given transition, [K] 
		up_level = 0;			// index of upper level of the transition, zero-based
		low_level = 0;			// index of lower level of the transition, zero-based
	}

	~radiative_transition()
	{
		blends.clear();
	}

	void add_overlapped_line(const size_t & blend_id, const double & vel_fac) 	// adds an element in the vector of lines that are overlapped with a given transition
	{
		blends.push_back(lineOverlapping(blend_id, vel_fac));
	}
};
