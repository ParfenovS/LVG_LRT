#pragma once
#include "level.h"
#include "transition.h"
#include "physconsts.h"
#include "modelPhysPars.h"
#include <variant>

class molModel				// stores molecular spectroscopic data
{
private:

	vector<vector <size_t> > C_counter;
	
	void compute_C(const size_t & id_up, const size_t & id_low, const size_t & num_of_coll_partners) 				// compute collisional rate of transition from lower (j) to upper (i) level; i and j begin from 0;
	{
	    if (C_counter[id_up][id_low] == num_of_coll_partners) {	// check if one already assigned Cji coef. and if one took into account all collisional agents;
			if (levels[id_low].g != 0.0 && modelPhysPars::Tks != 0.0) {			// check if one have all quantities to compute Cji;
				const double dE = levels[id_low].E - levels[id_up].E;
				coll_trans[id_low][id_up] = (levels[id_up].g / (double)levels[id_low].g) * coll_trans[id_up][id_low] * exp(dE / modelPhysPars::Tks * (SPEED_OF_LIGHT*PLANK_CONSTANT/BOLTZMANN_CONSTANT));
			} else coll_trans[id_low][id_up] = 0.0;
			// computing collisional part of diagonal elements of matrix A from the statistical equilibrium equation A*X=B
			// Cii = sum{k=1,Nlevel}(Cik)
			coll_trans[id_up][id_up] 	+= coll_trans[id_up][id_low];
			coll_trans[id_low][id_low]  += coll_trans[id_low][id_up];
		}
	}

public:
	
	size_t idspec; 							// index of a given molecular species
	size_t num_of_coll_partners;			// number of collisional partners
	vector<level> levels;					// array of levels; see level class definition in level.h;
	vector<vector <double> > coll_trans;	// 2d-array of collisional transitions;
	vector<radiative_transition> rad_trans; // array of radiative transitions; see rad_transition class definition in transition.h;
		
	// assign a value (par_value) to a given parameter (par_name) of a given transition
	void set_trans(const size_t & trans_id_1, const size_t & id_up_1, const size_t & id_low_1, const std::string_view par_name, const std::variant<size_t, double> par_value)
	{
		// converting to zero-based indexes
		const size_t trans_id	= trans_id_1	- 1;	// transition id
		const size_t id_up		= id_up_1		- 1;	// upper level id
		const size_t id_low		= id_low_1		- 1;	// lower level id

		if (levels.size() < 2) throw runtime_error("Input levels data first");
		if ( par_name != "C") {
			while (rad_trans.size() < trans_id_1) rad_trans.push_back(radiative_transition());
			rad_trans[trans_id].up_level = id_up;
			rad_trans[trans_id].low_level = id_low;
		}
/*
Compiler produces multiple copies of a template function depending on the types (T) of par_value. Copies differ from each other only by the type given by T. 
For the following if-statement, this may lead to assignment of variables with different types.
get<> allows to avoid compilation errors/warnings related to difference of types of operands during the assignment 
*/
		if      (par_name == "id")  rad_trans[trans_id].trans_id = get<size_t>(par_value);		// setting up transition id
		else if (par_name == "A")   rad_trans[trans_id].A = get<double>(par_value);				// setting up Einstein A coefficient
		else if (par_name == "Bul") rad_trans[trans_id].Bul = get<double>(par_value);			// setting up Einstein B coefficient
		else if (par_name == "nu") {															// setting up radiative transition frequency
			rad_trans[trans_id].nu = get<double>(par_value);
			//rad_trans[trans_id].nu = (levels[id_up].E - levels[id_low].E) * SPEED_OF_LIGHT;
			if (fabs(levels[id_up].E - levels[id_low].E) <= 1.e-20) { 	// treating the special case when the levels energy are (aslmost) equal however there is a radiative transition between them with non-neglible frequency
				if (rad_trans[trans_id].nu < 1.0) levels[id_up].E += rad_trans[trans_id].nu*1.e9 / SPEED_OF_LIGHT;	// modify levels energy according to radiative transition frequency in case of frequencies given in GHz
				else levels[id_up].E += rad_trans[trans_id].nu / SPEED_OF_LIGHT; 									// modify levels energy according to radiative transition frequency in case of frequencies given in Hz
			}
		}
		else if (par_name == "J")   rad_trans[trans_id].J = get<double>(par_value);				// setting up mean intensity of emission for a given radiative transition
		else if (par_name == "C") 																// setting up collisional rate for a given radiative transition
		{
			if (coll_trans.size() < levels.size()) {			// allocate array of collisional rates if it is hasn't been allocated already
				for (size_t k = 0; k < levels.size(); k++) {
					coll_trans.push_back( vector <double>() );
					C_counter.push_back( vector <size_t>() );
					for (size_t j = 0; j < levels.size(); j++) {
						coll_trans[k].push_back( 0.0 );
						C_counter[k].push_back( 0 );
					}
				}
			}
			if (C_counter[id_up][id_low] > num_of_coll_partners) throw runtime_error("error: C_counter > fraction_H2.size() in set_trans function in molModel.h");
			coll_trans[id_up][id_low] += modelPhysPars::fraction_H2[modelPhysPars::collPartCounter + C_counter[id_up][id_low]]*get<double>(par_value)*modelPhysPars::H2dens; 	// sum Cij for all collisional agents
			C_counter[id_up][id_low] += 1;
			compute_C(id_up, id_low, num_of_coll_partners);
		}
		else throw runtime_error("unknown par_name input parameter in set_trans function in molModel.h");
	}
	
	// assign a value (par_value) to a given parameter (par_name) of a given level
	void set_level(const size_t & id, const std::string_view par_name, const std::variant<size_t, double, string> par_value)
	{
		const size_t i = id - 1; 	// converting to zero-based index

		if (levels.size() < id) levels.push_back(level());

		if      (par_name == "id")  levels[i].id	 = get<size_t>(par_value);
		else if (par_name == "E")   levels[i].E		 = get<double>(par_value);
		else if (par_name == "g")   levels[i].g		 = static_cast<int>(round(get<double>(par_value)));
		else if (par_name == "pop") levels[i].pop	 = get<double>(par_value);
		//else if (par_name == "Q")   levels[i].Q_num	 = get<string>(par_value);
		else throw runtime_error("There is no such a field in the level class");
	}

	void check_data() 	// check whether frequencies are in Hz and Einstein B coefficients were set
	{
		// check frequencies
		bool freqsAreInGHz = false;
		for (size_t i = 0; i < rad_trans.size(); i++) {
			if (rad_trans[i].nu < 50.0) {
				freqsAreInGHz = true;
				break;
			}
		}
		if (freqsAreInGHz) {
			for (size_t i = 0; i < rad_trans.size(); i++) rad_trans[i].nu *= 1.e9; 	// converting from GHz to Hz
		}

		// check Einstein B coeffs.
		if (fabs(rad_trans[0].Bul) == 0.0 && fabs(rad_trans[rad_trans.size()-1].Bul) == 0.0) {
			for (size_t i = 0; i < rad_trans.size(); i++) {
				rad_trans[i].Bul = rad_trans[i].A * (SPEED_OF_LIGHT / rad_trans[i].nu) * (SPEED_OF_LIGHT / rad_trans[i].nu) / (PLANK_CONSTANT * rad_trans[i].nu) * 0.5;
				rad_trans[i].Blu = levels[rad_trans[i].up_level].g * rad_trans[i].Bul / (double)levels[rad_trans[i].low_level].g;
			}
		}

		C_counter.clear();
	}

	void compute_T(vector <vector <double> > & T, vector <vector <short> > & isItCollisionalDominated) // compute net population flow rates (see e.g. Sobolev & Deguchi 1994) and return it in T; should be called before pop_flow (declared in MonteCarlo.h);
	{
		double TC;
		if (T.size() < levels.size()) {
			T.resize(levels.size());
			isItCollisionalDominated.resize(levels.size());

			for (size_t k = 0; k < levels.size(); k++) {
				T[k].resize(levels.size()); 
				isItCollisionalDominated[k].resize(levels.size()); 
			}
			for (size_t i = 0; i < levels.size(); i++) {
				for (size_t j = 0; j < levels.size(); j++) {
					T[i][j] = 0.0;
					isItCollisionalDominated[i][j] = 0;
				}
			}
		}

		for (size_t i = 0; i < rad_trans.size(); i++) {
			const size_t & up = rad_trans[i].up_level;
			const size_t & low = rad_trans[i].low_level;
			T[up][low] = levels[up].pop*(rad_trans[i].A + rad_trans[i].Bul*rad_trans[i].J) - levels[low].pop*rad_trans[i].Blu*rad_trans[i].J;
			T[low][up] = - T[up][low];
		}
		
		for (size_t i = 0; i < levels.size(); i++) {
			for (size_t j = 0; j < levels.size(); j++) {
				TC = levels[i].pop*coll_trans[i][j] - levels[j].pop*coll_trans[j][i];
				T[i][j] += TC;
				if (fabs(TC) > fabs(0.9*T[i][j])) isItCollisionalDominated[i][j] = 1; //check whether the transition is collision dominated
			}
			T[i][i] = 0.0;
		}
		
	}

	molModel() noexcept
	{
		idspec = 0;
		num_of_coll_partners = 0;
	}

	molModel(const size_t & idspec_in) noexcept
	{
		idspec = idspec_in;
		num_of_coll_partners = 0;
	}
	
	~molModel()
	{
		levels.clear();
		coll_trans.clear();
		rad_trans.clear();
		C_counter.clear();
	}
};
