#pragma once
#include <string>

struct level			// stores information on a given rotational level
{
	size_t id;
	double E;
	int g;
	double pop;
	std::string Q_num;

	level() noexcept
	{
		id = 0; 		// level id, begins from 1;
		E = -1.0;		// energy [cm-1]; can be set in molModel.h with set_level function;
		g = -1;		// statistical weight
		pop = -1.0;		// population [1/cm3]
		Q_num = "";		// string with all quantum numbers for a given level
	}
};
