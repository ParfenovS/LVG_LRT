#pragma once
#include <vector>
#include <random>
#include <float.h>
#include <iostream>

using namespace std;

class MonteCarloSearchCycles			// search for the maser pumping cycles (see e.g. Sobolev 1986; Sobolev & Deguchi 1994)
{
private:

	struct cycle_class					// stores levels on a path from the i_lev to the k_lev levels
	{
		size_t cycle_counter = 1;
		bool A_is_not_positive_for_final_level = false;
		vector<size_t> lms;
		vector<bool> CollDom;
	};

	vector<cycle_class> cycles;			// stores all the found paths from the i_lev to the k_lev levels
	vector<vector<double> > P;
	vector<vector<double> > P0;

	size_t i_lev;						// initial level initialized in pop_flow()
	size_t k_lev;						// final level initialized in pop_flow()
	unsigned int seed;					// random seed
	size_t num_of_tries;				// number of tries to search cycles
	bool close_levels;					// = true - search for elementary cycles; = false - search for cycles with recursion
	size_t MAX_STEPS;					// maximum allowed number of steps in a path; calculated in compute_P();

	void compute_P()					// compute fractions of population flow from all to all levels;
	{
		if (P.size() < V.size()) {
			P.resize(V.size());
			for (size_t k = 0; k < V.size(); k++) P[k].resize(V.size());
			for (size_t i = 0; i < V.size(); i++) {
				for (size_t j = 0; j < V.size(); j++) {
					P[i][j] = 0.0;
				}
			}
		}

		for (size_t i = 0; i < V.size(); i++) {
			double sumV = 0.0;
			//double y, t, c = 0.0;
//#pragma omp parallel for reduction(+:sumV, c)			//parallelized Kahan summation of population flows V
			for (size_t q = 0; q < V.size(); q++) {
				/*if (V[i][q] > 0.0) {
					y = V[i][q] - c;
					t = sumV + y;
					c = (t - sumV) - y;
					sumV = t;	
				}*/
				if (V[i][q] > 0.0) sumV += V[i][q];
			}
			//sumV = sumV - c;

			sumV = 1. / sumV;
			for (size_t j = 0; j < V.size(); j++) {
				if (V[i][j] <= 0.) P[i][j] = 0.0;
				else P[i][j] = V[i][j] * sumV;
			}
		}
		P[i_lev][k_lev] = 0.0;
		P[i_lev][i_lev] = 0.0;
		P0 = P;
		MAX_STEPS = 10 * P.size();
	}

	double compute_A(const size_t &lm)			// compute population flow from lm level to another levels
	{
		double A = 0.0;
		for (size_t q = 0; q < P.size(); q++) A += P[lm][q];
		return A;
	}

	void close_level(const size_t &lm)			// closing level lm, i.e. zeroing fraction of population flow (Pqlm) from all levels to lm level; 
	{
		for (size_t k = 0; k < P.size(); k++) P[k][lm] = 0.0;
	}

	size_t compare(const double &AR, const size_t &lm)		// searching for the next level (lm1) to transition to such that sum(q=1,lm1-1)Plmq < AR < sum(q=1,lm1)Plmq; the transition is from lm to lm1 levels;
	{
		double sumP1 = 0.0; 
		double sumP2 = 0.0;
		//const double AR = A * R;
		for (size_t q = 0; q < P.size(); q++) {
			if (P[lm][q] > 0.0) {
				sumP2 += P[lm][q];
				if (sumP1 <= AR && AR <= sumP2) {		
					return q;
				}
				sumP1 = sumP2;
			}
		}
		return i_lev;
	}

	void search_for_cycles(mt19937 &gen, uniform_real_distribution<double> &dist)
	{
		size_t lm, lm1, steps;
		double A, R;
		compute_P();
		
		if (compute_A(i_lev) <= 0.0) throw runtime_error("A<=0 for the initial level, it's impossible to find a cycle");

		for (size_t c = 0; c < num_of_tries; c++) {
			steps = 0;
			cycle_class cycle;
			lm = i_lev;
			P = P0;
			A = compute_A(lm);
			cycle.lms.push_back(lm);

			while (true) {
				steps += 1;
				R = dist(gen);							// generating random R
				lm1 = compare(A*R, lm); 				// searching for the next level (lm1) such that sum(q=1,lm1-1)Plmq < AR < sum(q=1,lm1)Plmq; the transition is from lm to lm1 levels
				double Ai = compute_A(lm1);
				if (Ai <= 0.0 && lm1 != k_lev) {
					if (steps >= MAX_STEPS) {
						steps = 0;
						cycle.lms.clear();
						cycle.CollDom.clear();
						lm = i_lev;
						cycle.lms.push_back(lm);
						P = P0;
						A = compute_A(lm);
					}
					continue;
				}
				if (close_levels) close_level(lm);		// closing the current level
				cycle.CollDom.push_back(isItCollisionalDominated[lm][lm1]);
				cycle.lms.push_back(lm1);
				lm = lm1;
				if (close_levels) A = compute_A(lm);
				else A = Ai;
				if (lm == k_lev) {
					if (A <= 0.0) cycle.A_is_not_positive_for_final_level = true;
					steps = 0;
					break;
				}
			}

			if (cycles.size() == 0) {
				cycles.push_back(cycle);
				cycle.CollDom.clear();
				cycle.lms.clear();
				continue;
			}
			size_t same_cycle_id = 0;
			bool same_cycle = false;
			for (size_t n = 0; n < cycles.size(); n++) {		// check whether the cycle has been found at previous attempts
				if (cycle.lms.size() == cycles[n].lms.size()) {
					same_cycle_id = n;
					same_cycle = true;
					for (size_t m = 0; m < cycles[n].lms.size(); m++) {
						if (cycles[n].lms[m] != cycle.lms[m]) {
							same_cycle = false;
							break;
						}
					}
					if (same_cycle) break;
				}
			}
			if (!same_cycle) cycles.push_back(cycle);
			else cycles[same_cycle_id].cycle_counter += 1;
			cycle.CollDom.clear();
			cycle.lms.clear();
		}

		// calculate the relative efficiency and power of cycles
		vector <double> fE(cycles.size()); // relative efficiency
		vector <double> fW(cycles.size()); // power

		double sum_V = 0.0;
		for (size_t q = 0; q < V.size(); q++) {
			if (V[i_lev][q] > 0.0) sum_V += V[i_lev][q];
		}

		for (size_t c = 0; c < cycles.size(); c++)
		{
			fE[c] = 1;
			for (size_t m = 1; m < cycles[c].lms.size(); m++)
			{
				fE[c] *= P0[cycles[c].lms[m - 1]][cycles[c].lms[m]];
			}
			fW[c] = fE[c] * sum_V;
		}

		// sorting by efficiency
		vector <size_t> sort_indexes = argsort(fE);

		// output results
		for (size_t ic = 0; ic < cycles.size(); ic++)
		{
			size_t c = sort_indexes[ic];
			if (fabs(fE[c]) < 1.e-30) continue;
			cout << cycles[c].lms.size() << " " << cycles[c].cycle_counter << " " << fE[c] << " " << fW[c];
			if (cycles[c].A_is_not_positive_for_final_level) cout << " # warning A<=0 for the final level in this route";
			cout << "\n";

			for (size_t m = 0; m < cycles[c].lms.size(); m++)
			{
				if (cycles[c].CollDom[m] && cycles[c].lms[m] != k_lev) cout << "*";
				cout << cycles[c].lms[m] + 1 << "\t";
			}
			cout << "\n";
		}
	}

public:

	vector<vector <double> > V;							// net population flow rates; should be initialized from outside the class; allocated with compute_V function in molModel.h 
	vector<vector <bool> >isItCollisionalDominated;

	void pop_flow(const size_t &i, const size_t &k)		// main function that performs the search of paths or cycles from i to k levels; i and k begin from 1
	{
		this->i_lev = i - 1; 	// convert to zero-based indexes
		this->k_lev = k - 1; 	// convert to zero-based indexes
		mt19937 gen(seed);
		uniform_real_distribution<double> dist(0, 1.0f);
		search_for_cycles(gen, dist);
	}
	
	MonteCarloSearchCycles() noexcept
	{
		seed = static_cast<unsigned int>(time(NULL));
		close_levels = true;
		num_of_tries = 500;
		MAX_STEPS = 1;
		k_lev = 0;
		i_lev = 0;
	}

	MonteCarloSearchCycles(const int &seed, const size_t &num_of_tries, const int &close_levels)
	{
		if (seed < 0) this->seed = static_cast<unsigned int>(time(NULL));
		else this->seed =seed;
		if (num_of_tries<=0) this->num_of_tries = 500;
		else this->num_of_tries = num_of_tries;
		if (close_levels == 0) this->close_levels = false;
		else this->close_levels = true;
		MAX_STEPS = 1;
		k_lev = 0;
		i_lev = 0;
	}

	~MonteCarloSearchCycles()
	{
		V.clear();
		P.clear();
		P0.clear();
		isItCollisionalDominated.clear();
		cycles.clear();
	}
};

