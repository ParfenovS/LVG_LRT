#pragma once
#include "molModel.h"
#include "auxiliary.h"

class molDataRead			// reads and sets spectroscopic data
{
private:

	ifstream fin;
	istringstream sfin;
	
	double interp_C(double *T, double *C, const size_t & num_of_coll_temps) 	//collisional coefficient interpolation on temperature Temp
	{
		const double & Tk = modelPhysPars::Tks;
		if (Tk <= T[0]) return C[0];
		if (Tk >= T[num_of_coll_temps-1]) return C[num_of_coll_temps-1];
		for (size_t i = 1; i < num_of_coll_temps; i++) {
			if (T[i-1] <= Tk && Tk <= T[i]) return C[i-1] + (Tk - T[i-1]) * (C[i] - C[i-1]) / (T[i] - T[i-1]);  
		}
		return 0.0;
	}

/*
** reading file with molecular data in LAMDA format
*/
	void read_LAMDA(const string & filename, molModel *mod)
	{
		string str,Q_num;
		size_t id, id_up, id_low;
		double E, g, A, nu, T_C;
		size_t num_of_levels_lamda,num_of_rad_trans,num_of_coll_trans,num_of_coll_temps,num_of_coll_partners;

		fin.open(filename.c_str(), ios::in);

		if ( !fin.good() ) { 	// check if we found the file
			fin.close();
			throw runtime_error("can't find LAMDA file");
		}

		for (short i = 0; i < 5; i++) getline(fin,str);

		num_of_levels_lamda = readline<size_t>(fin); 	// readline function is taken from auxiliary.h
		
		getline(fin,str);

// reading levels data				
		for (size_t i = 0; i < num_of_levels_lamda; i++) {
			fin >> id >> E >> g;
			getline(fin,Q_num);
			mod->set_level(id, "id", id   );
			mod->set_level(id, "E",  E    ); 
			mod->set_level(id, "g",  g    ); 
			mod->set_level(id, "Q",  Q_num);
		}
	
		getline(fin,str);

		num_of_rad_trans = readline<size_t>(fin);
		
		getline(fin,str);

// radiative transitions
		for (size_t i = 0; i < num_of_rad_trans; i++) {
			fin >> id >> id_up >> id_low >> A >> nu;
			getline(fin,str);
			
			mod->set_trans(id, id_up, id_low, "id", id);
			mod->set_trans(id, id_up, id_low, "A",  A );
			mod->set_trans(id, id_up, id_low, "nu", nu);
		}

/*
** collisional coefficients
*/
		getline(fin,str);
		num_of_coll_partners = readline<size_t>(fin);
		for (size_t k = 0; k < num_of_coll_partners; k++) {
			for (short i = 0; i < 3; i++) getline(fin, str);

			num_of_coll_trans = readline<size_t>(fin);

			getline(fin, str);

			num_of_coll_temps = readline<size_t>(fin);

			getline(fin, str);

			double* T = new double[num_of_coll_temps];
			double* C = new double[num_of_coll_temps];

			getline(fin, str);
			sfin.clear();
			sfin.str(str);
			for (size_t i = 0; i < num_of_coll_temps; i++) { 	// temperatures for collisions with a given collisional agent;
				sfin >> T_C; T[i] = T_C;
			}

			getline(fin, str);

			for (size_t i = 0; i < num_of_coll_trans; i++) { 	// collisional transitions for a given collisional agent;
				getline(fin, str);
				sfin.clear();
				sfin.str(str);
				sfin >> id >> id_up >> id_low;
				for (size_t j = 0; j < num_of_coll_temps; j++) {
					sfin >> T_C; C[j] = T_C;
				}
				
				mod->set_trans(id, id_up, id_low, "C", interp_C(T, C, num_of_coll_temps));
			}
			delete[] T;
			delete[] C;
		}
				
		fin.close();
	}

/*
** reading file with levels population
*/
	void read_pops(const string & filename, molModel *mod)
	{
		string str;

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

/*
** reading file with the emission mean intensity
*/
	void read_J(const string & filename, molModel *mod)
	{
		size_t id_up, id_low;
		double J;
		string str;

		fin.open(filename.c_str(), ios::in);
		sfin.clear();

		if (!fin.good())	// check if we found the file
		{
			fin.close();
			throw runtime_error("can't find file with the mean intensity");
		}
		size_t id = 1;
		while (!fin.eof())
		{
			getline(fin, str);
			str = trim(str);
			if (str.size() == 0) break;
			sfin.str(str);
			sfin >> id_up >> id_low >> J;
			sfin.clear();
			mod->set_trans(id, id_up, id_low, "J", J);
			mod->set_trans(id, id_low, id_up, "J", J);
			id += 1;
		}

		fin.close();
	}

public:

	molDataRead() noexcept
	{
	}

	~molDataRead()
	{
	}

	void read_data(const string & filename_lamda, const string & filename_pops, molModel *mod)
	{
		read_LAMDA(filename_lamda, mod);
		read_pops(filename_pops, mod);
		mod->check_data(); // check whether frequencies were given in Hz; computes Einstein B coefficients if needed
	}

	void read_data(const string & filename_lamda, const string & filename_pops, const string & filename_J, molModel *mod)
	{
		read_LAMDA(filename_lamda, mod);
		read_pops(filename_pops, mod);
		mod->check_data(); // check whether frequencies were given in Hz; computes Einstein B coefficients if needed
		read_J(filename_J, mod);
	}
};
