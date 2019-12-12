#pragma once
#include "math.h"
#include <float.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "physconsts.h"

using namespace std;

inline double planck_function(const double & T, const double & nu)
{
    constexpr double inv_squared_c = 2. / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    const double exp_factor = (PLANK_CONSTANT * nu) / (BOLTZMANN_CONSTANT * T);
    if (fabs(exp_factor) > 1.e-4) return PLANK_CONSTANT * (nu*nu*nu*inv_squared_c) / ( exp(exp_factor) - 1.0 );
    else return (nu*nu*inv_squared_c) * BOLTZMANN_CONSTANT * T;
}

template <typename T, typename Tfin>
inline T readline(Tfin & fin)       // auxiliary function to read lines with empty spaces after symbols and before \n
{
    T a;
    string str;
    istringstream sfin;
    getline(fin, str);
    sfin.str(str);
    sfin >> a;
    sfin.clear();
    return a;
}

inline string trim(string & str)    // remove left and right spaces in a string
{
    str.erase(0, str.find_first_not_of(' '));       // prefixing spaces
    str.erase(str.find_last_not_of(' ') + 1);       // surfixing spaces
    return str;
}

template <typename T>
inline void two_sum(T & err, T & sum, const T & input)  // Knuth's Two-Sum algorithm
{
    const T a = sum;
    const T b = input;

    const T s = a + b;
    const T a0 = s - b;
    const T b0 = s - a0;
    const T da = a - a0;
    const T db = b - b0;
    err = da + db;
    sum = s;
}

template <typename T>
inline T cascade_summation(T & err, T & sum, const T & input)   // cascaded summation algorithm of Rump, Ogita, and Oishi 
{
    T err_i = 0.0;
    two_sum(err_i, sum, input);
    err += err_i;
    return sum;
}
